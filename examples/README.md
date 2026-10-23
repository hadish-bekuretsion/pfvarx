# Examples

This directory contains a worked, runnable example that requires only the
Python core (no external bioinformatics binaries).

## Run it

```bash
bash examples/run_demo.sh
```

The script:

1. prints the version and release status;
2. generates a synthetic reference and BAM with three variants at allele
   fractions 0.50, 0.20, and 0.08;
3. runs the statistical caller and writes a VCF plus a 46-column evidence table;
4. runs unified QC;
5. attempts report export; inspect any export error, as optional exporter dependencies may be unavailable for a platform.

Outputs are written under `examples/output/` by default; pass a different
directory as the first argument.

## What to look for

- `examples/output/results/DEMO/DEMO.evidence.tsv` — the scientific record.
  The 0.50 and 0.20 variants should be `ALT_DETECTED_EXPLORATORY`; the 0.08
  variant falls in the minority bin.
- `examples/output/results/DEMO/DEMO.calls.vcf` — `QUAL` is `.` for every
  record because the policy is provisional.
- `examples/output/results/run_qc.json` — the run-level QC summary.

## Real data

For real FASTQ/BAM/VCF data, follow `docs/USAGE.md`. You must supply a pinned
reference assembly and record its SHA-256; see
`config/resources/reference_manifest.yaml`.
