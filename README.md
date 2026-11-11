# NCBI GFF3 to VEP

Convert and coordinate-sort NCBI GFF3 annotations for use with [Ensembl Variant Effect Predictor (VEP)](https://www.ensembl.org/info/docs/tools/vep/index.html).

## Features

- Repairs common NCBI GFF3 gene → CDS structures by generating transcript features where necessary.
- Generates exon features from CDS intervals for transcripts lacking exons.
- Reconnects CDS `Parent` attributes to generated transcripts.
- Sorts output by sequence ID and coordinate, prioritizing parents at equal starting positions.
- Accepts uncompressed or gzip-compressed input GFF3.
- Writes plain-text GFF3 or BGZF-compressed GFF3 and optionally creates a tabix index.
- Reports missing parent links and unresolved CDS-to-transcript relationships.

## Requirements

- Python 3.8+ (standard library only)
- `bgzip` on `PATH` for `.gz` output
- `tabix` on `PATH` when using `--index`

## Usage

```bash
python convert_ncbi_gff_for_vep.py genomic.gff.gz genomic.vep.gff3.gz --index
```

For plain-text output without indexing:

```bash
python convert_ncbi_gff_for_vep.py genomic.gff genomic.vep.gff3
```

Other flags:

```text
--keep-other         Retain other GFF3 feature types
--skip-pseudogenes   Do not generate synthetic mRNA under pseudogenes
```

Then annotate variants with VEP:

```bash
vep --input_file variants.vcf.gz \
    --output_file variants.vep.vcf.gz \
    --vcf --compress_output bgzip --offline \
    --gff genomic.vep.gff3.gz \
    --fasta reference.fasta \
    --symbol --protein --biotype
```

**Important:** The GFF3 sequence IDs must match the FASTA/VCF contig names. Synthesized exons are based on annotated CDS intervals and do not infer missing UTRs. Inspect warnings and validate the output on your specific annotation, particularly for complex eukaryotic transcript models. Sorting is lexicographic by sequence ID.

## License

No license has been selected yet. Add one if you plan to distribute under an open-source license.
