# PfVarX Quickstart

This is a concise walkthrough of the two-stage user workflow. The [README](../README.md) contains the complete configuration, input-format, and output reference.

> **Research use only.** PfVarX is not clinically validated and must not inform diagnosis or treatment decisions.

## 0. Install Pixi environments

From the project directory:

```bash
curl -fsSL https://pixi.sh/install.sh | bash
git clone https://github.com/hadish-bekuretsion/PfVarX_Pipeline.git
cd PfVarX_Pipeline
pixi install --locked -e default
pixi install --locked -e workflow   # needed for FASTQ, Snakemake, and reference indexing
```

No separate Conda environment is required. The `default` environment provides the core CLI; the `workflow` environment adds Snakemake and external bioinformatics tools.

## 1. Stage 1 — configure and create input folders

Copy the example and set the reference, sample-sheet, and input folder paths. Keep `outputs.root` separate from all input paths.

```bash
cp config/config.yaml config/my-run.yaml
```

For example, configure `reference.fasta: data/reference/Pf3D7.fa`, `inputs.samplesheet: data/prepared/samples.tsv`, `inputs.raw_dir: data/raw`, and `outputs.root: results/my-run`. Set `outputs.logs_dir` or `outputs.work_dir` for independent locations; if null, each defaults beneath `outputs.root`.

Run input initialization:

```bash
pixi run -e default pfvarx init --config config/my-run.yaml
```

This creates only configured input directories and parents for declared input resource paths. It does not create results, logs, or work folders and does not populate the input files.

## 2. Stage 1 — place and prepare input files

Copy the matching reference FASTA to the configured path and generate its checksum and `.fai` index:

```bash
sha256sum data/reference/Pf3D7.fa
pixi run -e workflow samtools faidx data/reference/Pf3D7.fa
```

Set `reference.sha256` to the resulting 64-character digest. For paired FASTQ only, create the BWA-MEM2 index:

```bash
pixi run -e workflow bwa-mem2 index data/reference/Pf3D7.fa
```

Create the configured sample sheet as tab-separated text. Paths in it are relative to the sample-sheet directory unless absolute.

```tsv
sample_id	input_type	assay	fastq_1	fastq_2	bam	vcf	metadata_json
PF001	fastq	wgs	../raw/PF001_R1.fastq.gz	../raw/PF001_R2.fastq.gz			
PF002	bam	wgs			../raw/PF002.bam		
```

At least `sample_id` and `input_type` are required. FASTQ requires paired `fastq_1`/`fastq_2`; BAM/CRAM reuses the alignment; VCF/BCF/gVCF imports existing variants. The pipeline does not generate a `.dict` sequence dictionary. A matching GFF3 is required only for `drug_resistance`; metadata fields are needed for `population_structure` and `spatial_temporal`.

## 3. Stage 2 — select objectives and preflight

Set objectives in `config/my-run.yaml` (for example `objectives: [diversity, mixed_infection]`) and run preflight without `--out`:

```bash
pixi run -e default pfvarx validate \
  --config config/my-run.yaml \
  --samplesheet data/prepared/samples.tsv
```

The command validates/routs the inputs and previews the directories Stage 2 will create; it does not create the results, log, or work directories. Use `pixi run -e default pfvarx objectives` to list the available objectives and prerequisites.

## 4. Stage 2 — execute

For paired FASTQ, use the workflow environment because the core CLI does not align raw reads:

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

For BAM/CRAM or VCF/BCF/gVCF, the direct CLI is supported:

```bash
pixi run -e default pfvarx run \
  --config config/my-run.yaml \
  --samplesheet data/prepared/samples.tsv \
  --objectives diversity
```

During Stage 2, the pipeline creates final result, log, and work directories from `outputs.root`, `outputs.logs_dir`, and `outputs.work_dir`. The latter two default to `outputs.root/logs` and `outputs.root/intermediate`.

## 5. Review the run

Start with `run_manifest.json` for objective statuses and artifact paths, then review each evidence TSV/VCF and the relevant stage logs. `completed` is software completion, not scientific validation. `unavailable`, `ALT_UNCERTAIN`, and `NO_CALL` are not evidence of a reference genotype or absence of a biological signal.

See the [README](../README.md), [installation guide](INSTALL.md), [objective guide](OBJECTIVES.md), and [validation record](VALIDATION.md) for more detail.
