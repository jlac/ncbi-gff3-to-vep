#!/usr/bin/env python3
"""Repair NCBI/VEuPathDB GFF3 transcript models and coordinate-sort for Ensembl VEP.

Python >=3.8; no Python package dependencies. For .gz output requires bgzip;
--index additionally requires tabix. Uses memory proportional to GFF feature count.
"""
import argparse
import collections
import gzip
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from urllib.parse import quote, unquote

TRANSCRIPT_TYPES = {'mRNA', 'transcript', 'ncRNA', 'lnc_RNA', 'lncRNA',
                    'rRNA', 'tRNA', 'snRNA', 'snoRNA', 'miRNA', 'scRNA',
                    'antisense_RNA', 'pseudogenic_transcript', 'primary_transcript',
                    'SRP_RNA', 'RNase_P_RNA', 'RNase_MRP_RNA', 'tmRNA'}
GENE_TYPES = {'gene', 'protein_coding_gene', 'ncRNA_gene', 'pseudogene'}
KEEP_TYPES = TRANSCRIPT_TYPES | GENE_TYPES | {'CDS', 'exon',
                                 'five_prime_UTR', 'three_prime_UTR', 'start_codon',
                                 'stop_codon'}
RANK = {**{x: 0 for x in GENE_TYPES}, **{x: 1 for x in TRANSCRIPT_TYPES},
        'exon': 2, 'five_prime_UTR': 3, 'three_prime_UTR': 3,
        'start_codon': 4, 'CDS': 4, 'stop_codon': 4}


def attrs_parse(s):
    ans = {}
    if s == '.':
        return ans
    for field in s.split(';'):
        if not field:
            continue
        if '=' not in field:
            continue
        k, v = field.split('=', 1)
        ans[unquote(k)] = unquote(v)
    return ans


def attrs_format(d):
    safe = ':@|^*+-._'  # commas are preserved in Parent, handled specially
    parts = []
    for k, v in d.items():
        value = ','.join(quote(p, safe=safe) for p in v.split(',')) if k == 'Parent' else quote(v, safe=safe)
        parts.append(quote(k, safe=safe) + '=' + value)
    return ';'.join(parts) if parts else '.'


def parents(f):
    return [p for p in f.attrs.get('Parent', '').split(',') if p]


@dataclass
class Feature:
    seqid: str
    source: str
    type: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    attrs: dict
    ordinal: int

    def line(self):
        return '\t'.join([self.seqid, self.source, self.type, str(self.start),
                          str(self.end), self.score, self.strand, self.phase,
                          attrs_format(self.attrs)]) + '\n'


def load(path, keep_other, warnings):
    features = []
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8') as fh:
        for lineno, raw in enumerate(fh, 1):
            s = raw.rstrip('\r\n')
            if s == '##FASTA':
                warnings.append('Skipped embedded FASTA section')
                break
            if not s or s.startswith('#'):
                continue
            parts = s.split('\t')
            if len(parts) != 9:
                raise ValueError(f'Line {lineno}: expected 9 tab-separated fields, found {len(parts)}')
            seqid, source, typ, start, end, score, strand, phase, attr = parts
            try:
                start, end = int(start), int(end)
            except ValueError:
                raise ValueError(f'Line {lineno}: invalid coordinates') from None
            if start < 1 or end < start:
                raise ValueError(f'Line {lineno}: invalid interval {start}-{end}')
            if not keep_other and typ not in KEEP_TYPES:
                continue
            features.append(Feature(seqid, source, typ, start, end, score, strand,
                                    phase, attrs_parse(attr), lineno))
    return features


def unique_id(proposed, known):
    result = proposed
    i = 2
    while result in known:
        result = f'{proposed}_{i}'
        i += 1
    known.add(result)
    return result


def create_feature(template, typ, start, end, attrs, ordinal, phase='.'):
    return Feature(template.seqid, template.source, typ, start, end, '.',
                   template.strand, phase, attrs, ordinal)


def convert(features, warnings, skip_pseudogenes=False):
    all_ids = collections.defaultdict(list)
    for f in features:
        if f.attrs.get('ID'):
            all_ids[f.attrs['ID']].append(f)
    duplicated = [k for k, v in all_ids.items() if len(v) > 1]
    if duplicated:
        warnings.append(f'{len(duplicated)} duplicate IDs in input; disambiguation not attempted')
    by_id = {k: v[0] for k, v in all_ids.items()}
    known = set(by_id)
    stats = collections.Counter()
    next_ordinal = max((f.ordinal for f in features), default=0) + 1
    added = []

    # NCBI bacterial GFF commonly gives CDS Parent=gene directly, with no mRNA.
    # PlasmoDB protein_coding_gene is also recognized as a gene feature.
    # One synthetic transcript per parent gene covers its directly attached CDSs.
    gene_cdss = collections.defaultdict(list)
    for f in features:
        if f.type != 'CDS':
            continue
        for parent in parents(f):
            p = by_id.get(parent)
            if p and p.type in GENE_TYPES:
                gene_cdss[parent].append(f)

    for gene_id, cdss in gene_cdss.items():
        gene = by_id[gene_id]
        if gene.type == 'pseudogene' and skip_pseudogenes:
            stats['skipped_pseudogene_groups'] += 1
            continue
        # Skip pseudogene CDSs by default? No: explicit opt-in skipping only.
        tr_id = unique_id(f'{gene_id}.vep_transcript', known)
        start = min(f.start for f in cdss)
        end = max(f.end for f in cdss)
        # Use CDS bounds for the synthetic coding transcript, not necessarily the
        # whole gene span (which may include promoter regions).
        attrs = {'ID': tr_id, 'Parent': gene_id}
        for tag in ('Name', 'gene', 'gene_name', 'product', 'locus_tag', 'Dbxref',
                    'ebi_biotype', 'biotype', 'gene_biotype'):
            if gene.attrs.get(tag):
                attrs[tag] = gene.attrs[tag]
        if gene.attrs.get('pseudo') == 'true' or gene.type == 'pseudogene':
            attrs['pseudo'] = 'true'
        tr = create_feature(gene, 'mRNA', start, end, attrs, next_ordinal)
        next_ordinal += 1
        added.append(tr)
        by_id[tr_id] = tr
        stats['created_transcripts'] += 1
        for cds in cdss:
            existing = parents(cds)
            updated = [tr_id if x == gene_id else x for x in existing]
            cds.attrs['Parent'] = ','.join(dict.fromkeys(updated))
            stats['reparented_cds_links'] += 1

    # Recover exons where a transcript has CDS but no exon features. Keep
    # existing exon coordinates (including UTRs) untouched when present.
    exons_by_tx = collections.defaultdict(list)
    cdss_by_tx = collections.defaultdict(list)
    for f in features:
        if f.type not in {'CDS', 'exon'}:
            continue
        for p in parents(f):
            tx = by_id.get(p)
            if tx and tx.type in TRANSCRIPT_TYPES:
                (exons_by_tx if f.type == 'exon' else cdss_by_tx)[p].append(f)

    for tx_id, cdss in cdss_by_tx.items():
        if exons_by_tx.get(tx_id):
            continue
        tx = by_id[tx_id]
        # Deduplicate exact intervals, but do not coalesce adjacent distinct
        # CDS intervals: exon structure/phase should not be guessed.
        spans = sorted({(f.start, f.end) for f in cdss})
        for num, (start, end) in enumerate(spans, 1):
            exon_id = unique_id(f'{tx_id}.vep_exon{num}', known)
            attrs = {'ID': exon_id, 'Parent': tx_id}
            added.append(create_feature(tx, 'exon', start, end, attrs, next_ordinal))
            next_ordinal += 1
            stats['created_exons'] += 1

    result = features + added
    # Lexicographic seqid order matches standard VEP / GNU sort example;
    # tie break by parent-before-child priority, then end, then input order.
    # For equal starts, parent before child is also preferred for GFF consumers.
    result.sort(key=lambda f: (f.seqid, f.start,
                               RANK.get(f.type, 9), f.end, f.ordinal))

    prev = None
    for f in result:
        key = (f.seqid, f.start)
        if prev is not None and key < prev:
            raise AssertionError('Sorting invariant failed')
        prev = key
    stats['output_features'] = len(result)
    stats['input_features_retained'] = len(features)
    stats['output_cds'] = sum(f.type == 'CDS' for f in result)
    orphans = sum(1 for f in result for p in parents(f) if p not in by_id)
    if orphans:
        warnings.append(f'{orphans} Parent links reference IDs absent from retained output')
    # Existing CDS without exons and other nonstandard relationships are flagged.
    unresolved = sum(1 for f in result if f.type == 'CDS' and not any(
        by_id.get(p) and by_id[p].type in TRANSCRIPT_TYPES for p in parents(f)))
    if unresolved:
        warnings.append(f'{unresolved} CDS features remain without a transcript Parent')
    return result, stats


def write_output(features, output, index):
    compressed = output.endswith('.gz')
    if index and not compressed:
        raise ValueError('--index requires an output filename ending in .gz')
    if compressed:
        if not shutil.which('bgzip'):
            raise RuntimeError('bgzip is required for .gz output (install htslib/samtools)')
        with open(output, 'wb') as outfile:
            proc = subprocess.Popen(['bgzip', '-c'], stdin=subprocess.PIPE,
                                    stdout=outfile, text=True)
            try:
                for f in features:
                    proc.stdin.write(f.line())
                proc.stdin.close()
                if proc.wait() != 0:
                    raise RuntimeError('bgzip failed')
            except BaseException:
                proc.kill()
                proc.wait()
                raise
    else:
        with open(output, 'w', encoding='utf-8') as fh:
            for f in features:
                fh.write(f.line())
    if index:
        if not shutil.which('tabix'):
            raise RuntimeError('tabix is required for --index (output .gz was created)')
        subprocess.run(['tabix', '-f', '-p', 'gff', output], check=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('input', help='NCBI/VEuPathDB GFF3 or GFF3.gz input')
    ap.add_argument('output', help='sorted .gff3 or BGZF-compressed .gff3.gz output')
    ap.add_argument('--index', action='store_true', help='create .tbi index (requires .gz output and tabix)')
    ap.add_argument('--keep-other', action='store_true', help='retain other feature types (e.g. region, match)')
    ap.add_argument('--skip-pseudogenes', action='store_true', help='do not synthesize mRNA under pseudogenes')
    args = ap.parse_args()
    if os.path.abspath(args.input) == os.path.abspath(args.output):
        ap.error('input and output must be different paths')
    warnings = []
    try:
        features = load(args.input, args.keep_other, warnings)
        if not features:
            raise ValueError('No GFF3 features retained from input')
        fixed, counts = convert(features, warnings, args.skip_pseudogenes)
        write_output(fixed, args.output, args.index)
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as exc:
        ap.exit(1, f'ERROR: {exc}\n')
    print('Output:', args.output, file=sys.stderr)
    if args.index:
        print('Index:', args.output + '.tbi', file=sys.stderr)
    for k, v in sorted(counts.items()):
        print(f'{k}: {v}', file=sys.stderr)
    for w in warnings:
        print('WARNING:', w, file=sys.stderr)


if __name__ == '__main__':
    main()
