# Installation and first run

PfVarX uses **Pixi** as the primary environment and task manager. The committed `pixi.lock` declares `linux-64`, `linux-aarch64`, and Intel macOS (`osx-64`). The 2026-10-23 audit installed both environments and ran software tests and a synthetic FASTQ workflow on Linux ARM64; it did not validate biological performance or every declared platform. Native macOS ARM64 and Windows are not declared targets. Windows users may use WSL2 where appropriate. No separate user-managed Conda environment is required.

## 1. Install Pixi environments

From the project directory:

```bash
curl -fsSL https://pixi.sh/install.sh | bash
git clone https://github.com/hadish-bekuretsion/PfVarX_Pipeline.git
cd PfVarX_Pipeline
pixi install --locked -e default
pixi run -e default status
pixi run -e default test
```

Install the `workflow` environment for FASTQ preprocessing, reference indexing, and Snakemake execution:

```bash
pixi install --locked -e workflow
pixi run -e workflow snakemake --version
```

`default` contains the Python CLI, caller, objectives, and reporting dependencies. `workflow` adds Snakemake and the declared external bioinformatics tools, including BWA-MEM2, SAMtools, fastp, Picard, and fgbio. fgbio is invoked only in the configured FASTQ UMI branch; its presence does not mean that UMI performance is validated. The optional cluster executor plugins are installed, but Slurm/Kubernetes deployments are not verified by the audit.

## 2. Stage 1 — create input folders only

Copy the example configuration and edit the data locations, reference identity, sample-sheet path, and Stage 2 output destinations:

```bash
cp config/config.yaml config/my-run.yaml
```

The example uses these fields:

```yaml
reference:
  assembly_id: Pf3D7_v68
  fasta: data/reference/Pf3D7.fa
  sha256: REPLACE_WITH_64_HEX_SHA256

inputs:
  root: data
  raw_dir: data/raw
  prepared_dir: data/prepared
  metadata_dir: null
  samplesheet: data/prepared/samples.tsv

outputs:
  root: results/my-run
  logs_dir: null
  work_dir: null
```

Paths in configuration are interpreted from the directory where you invoke the command. Choose a unique `outputs.root` for each analysis. `logs_dir` and `work_dir` are optional independent locations; if unset they default to `<outputs.root>/logs` and `<outputs.root>/intermediate`.

Run the Stage 1 initializer:

```bash
pixi run -e default pfvarx init --config config/my-run.yaml
```

It creates configured `raw_dir` and `prepared_dir`, any explicitly configured `metadata_dir`, and parent directories for the declared sample sheet, FASTA, and configured optional reference/annotation/metadata resources. It does not create output, results, log, or run-work directories and does not create or download the input files. Input and output path overlap is rejected.

Expected input layout with the example values:

```text
data/
├── raw/                 # place reads, BAM/CRAM, or VCF/BCF/gVCF here
├── prepared/
│   └── samples.tsv      # create manually
└── reference/
    └── Pf3D7.fa         # supply manually
```

If using an external metadata table, set `inputs.metadata_dir` (for example `data/metadata`) or set `downstream.metadata_table` to its path. `init` then creates the configured folder or file parent. The packaged resistance catalogue is used by default; no separate copy is needed unless you configure an override.

## 3. Stage 1 — prepare files manually

PfVarX does not bundle reference genomes. Obtain the FASTA for the exact assembly used for the input alignments/calls; copy it to `reference.fasta`, then generate its checksum and FASTA index:

```bash
sha256sum data/reference/Pf3D7.fa
pixi run -e workflow samtools faidx data/reference/Pf3D7.fa
```

Paste the checksum into `reference.sha256`. It must be the exact 64-character hexadecimal digest, not `null` or an unresolved placeholder. `pfvarx validate`, `pfvarx run`, direct `pfvarx call`/`import`/`annotate`, and the Snakemake validation stage recompute the FASTA digest and stop before analysis if it is a placeholder or does not match. This guards the reference identity; `inputs.require_checksums` does not currently verify checksums for sample-sheet read/alignment/variant files.

For the **paired FASTQ route only**, prepare BWA-MEM2 index files in the reference location:

```bash
pixi run -e workflow bwa-mem2 index data/reference/Pf3D7.fa
```

The current pipeline does not generate or consume a Picard/GATK sequence dictionary (`.dict`); it is not required unless another tool you use requires it. Existing BAM/CRAM and VCF/BCF inputs do not go through FASTQ alignment. Keep their reference assembly and contig names consistent with `reference.fasta`.

Create a tab-delimited `samples.tsv` with at least `sample_id` and `input_type`. Accepted input types are `fastq`, `bam`, `cram`, `vcf`, `bcf`, `gvcf`, and `auto`. Paired reads require both `fastq_1` and `fastq_2`. Paths in the sheet are resolved relative to the sheet's directory unless absolute paths or the CLI `--base-dir` are used.

```tsv
sample_id	input_type	assay	fastq_1	fastq_2	bam	vcf	metadata_json
PF001	fastq	wgs	../raw/PF001_R1.fastq.gz	../raw/PF001_R2.fastq.gz			
PF002	bam	wgs			../raw/PF002.bam		
PF003	vcf	wgs				../raw/PF003.vcf.gz	{"country":"Ethiopia","collection_date":"2023-11-04"}
```

For `drug_resistance`, provide a GFF3 gene model matching the FASTA and set `annotation.gff`. For `population_structure`, provide metadata including `country`; for `spatial_temporal`, include both `country` and `collection_date`. A separate metadata table is keyed by `sample_id`; alternatively, small metadata values can be included in `metadata_json`. `relatedness` and `phylogenetics` require at least two samples.

## 4. Stage 2 — validate and execute

Set the selected objective list in the configuration, e.g.:

```yaml
objectives:
  - diversity
  - mixed_infection
```

Before execution, run preflight validation. Without `--out`, validation prints the route/objective plan and does not create result directories:

```bash
pixi run -e default pfvarx validate \
  --config config/my-run.yaml \
  --samplesheet data/prepared/samples.tsv
```

For paired FASTQ data, use the Snakemake workflow because `pfvarx run` does not align raw reads:

```bash
pixi run -e workflow snakemake \
  --snakefile workflow/Snakefile \
  --configfile config/my-run.yaml \
  --profile workflow/profiles/local --dry-run

pixi run -e workflow snakemake \
  --snakefile workflow/Snakefile \
  --configfile config/my-run.yaml \
  --profile workflow/profiles/local
```

For BAM/CRAM or pre-called VCF/BCF/gVCF, the default CLI supports direct processing/import:

```bash
pixi run -e default pfvarx run \
  --config config/my-run.yaml \
  --samplesheet data/prepared/samples.tsv \
  --objectives diversity,relatedness
```

After preflight, Stage 2 creates the necessary output, log, and work directories automatically. Final results are written under `outputs.root`; stage/sample logs go to `outputs.logs_dir`; FASTQ trimming/alignment intermediates go to `outputs.work_dir`. If the latter paths are unset, they are created under the results root. Snakemake also creates its per-run validation records, evidence, reports, and manifest as targets execute.

`outputs.overwrite` defaults to `false`. The CLI refuses known file collisions; setting it to `true` permits replacement of those corresponding files but does not clean the tree. Use a new `outputs.root` for a new run.

## 5. Outputs and interpretation

`run_manifest.json` is the main record of objective status, artifacts, warnings, and errors. The direct CLI also writes `run_summary.json`; Snakemake writes evidence/VCFs under `<outputs.root>/<sample>/`, while direct BAM/CRAM or VCF import uses `<outputs.root>/variants/<sample>/`. Selected analyses are under `<outputs.root>/objectives/<objective>/`, with reports under `<outputs.root>/reports/`. Inspect logs at the configured log location.

`ALT_DETECTED_EXPLORATORY` is not a validated call. `ALT_UNCERTAIN`, `NO_CALL`, `failed`, and `unavailable` must not be read as evidence of reference genotype or absence of resistance. Imported variants do not acquire PfVarX-computed likelihoods. Software tests and synthetic runs check implementation/dataflow, not analytical validity.

## 6. Optional container path

`container/Dockerfile` and `container/pfvarx.def` are optional container recipes. Routine environment installation and task execution should use Pixi. The Dockerfile's build-time compatibility specification is separate from the Pixi workflow and is not the primary installation route.

## 7. Troubleshooting

| Symptom | Likely action |
|---|---|
| Pixi reports a stale lock | Review `pixi.toml`, run `pixi lock`, inspect the lockfile diff, and retain the updated `pixi.lock` with the manifest change. |
| `pfvarx` or a bioinformatics executable is missing | Install and invoke the correct Pixi environment (`default` for the core CLI; `workflow` for Snakemake and external tools). |
| `pfvarx init` rejects an input path | Move the input location outside `outputs.root`, `outputs.logs_dir`, or `outputs.work_dir`. |
| Reference checksum/FASTA validation fails | Confirm the exact reference path, assembly ID, checksum, and generated `.fai`; use the same assembly used for BAM/VCF inputs. |
| BWA-MEM2 reports missing index data | Run `pixi run -e workflow bwa-mem2 index <reference.fa>` before the FASTQ workflow. |
| Sample files appear missing | Paths are relative to the sample-sheet directory unless absolute or overridden with `--base-dir`. |
| FASTQ is not aligned by the CLI | Use the Snakemake command; the in-process `pfvarx run` does not align FASTQ. |
| A multi-sample objective refuses to run | Provide at least two valid samples or remove the objective. |
| Existing outputs block a run | Select a fresh `outputs.root`; only enable overwrite when intentional, noting that old files are not deleted. |

See [`TROUBLESHOOTING.md`](TROUBLESHOOTING.md), [`PARAMETERS.md`](PARAMETERS.md), [`OBJECTIVES.md`](OBJECTIVES.md), and [`VALIDATION.md`](VALIDATION.md) for further detail and limitations.
