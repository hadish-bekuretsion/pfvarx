# Usage

> New to PfVarX? Start with [`QUICKSTART.md`](QUICKSTART.md) and
> [`OBJECTIVES.md`](OBJECTIVES.md). This page is the detailed reference.

This guide walks through a complete PfVarX run: preparing a samplesheet,
configuring the run, validating and routing inputs, executing the pipeline,
and interpreting the outputs. Every command shown here is implemented in
`src/pfvarx/cli.py` and can be inspected with `pixi run -e default pfvarx
<command> --help`. CLI examples below omit that prefix for readability; run
them through `pixi run -e default pfvarx ...`. Use the `workflow` environment
for Snakemake and its bioinformatics executables.

---

## 1. Command overview

```
pfvarx init         Create input directories only (Stage 1)
pfvarx objectives    List the available analysis objectives
pfvarx validate      Validate inputs, route samples and preview the plan
pfvarx run           Execute the selected objectives
pfvarx call          Call variants for one BAM/CRAM
pfvarx calibrate     Fit error rate and over-dispersion
pfvarx benchmark     Run benchmark and ablation study
pfvarx report        Export reports or launch the dashboard
pfvarx parameters    Print or export the parameter dictionary
pfvarx status        Print the release-status matrix
pfvarx demo          Generate synthetic data and run an end-to-end call
pfvarx import        Import a pre-called VCF/BCF/gVCF (no realignment)
pfvarx annotate      Annotate an evidence table
pfvarx downstream    Diversity and mixed-infection summaries
pfvarx relatedness   Pairwise relatedness and NJ tree
pfvarx qc            Unified QC across a results tree
```

Configuration and validation options vary by command; consult that command's
`--help` output. Options such as `--json` and `--check-paths` are not accepted
by every subcommand. `init` creates only configured input folders and resource
parents; results, logs, and work directories are created during Stage 2.

---

## 2. The samplesheet

The samplesheet is a TSV (or CSV) with a fixed header. One row per sample:

```
sample_id   input_type   assay   fastq_1   fastq_2   bam   vcf
```

| Column | Meaning |
|---|---|
| `sample_id` | Unique identifier; used for all output paths |
| `input_type` | `fastq`, `bam`, `cram`, `vcf`, `bcf`, `gvcf`, or `auto` |
| `assay` | `wgs`, `tbs_no_umi`, or `tbs_umi` |
| `fastq_1`, `fastq_2` | Paired-end FASTQ paths (may be `.gz`) |
| `bam` | BAM/CRAM path for reuse |
| `vcf` | VCF/BCF/gVCF path for import |

Only the columns relevant to a sample's `input_type` need to be populated.
`auto` triggers content-based detection. A reference FASTA placed in a sample
slot is rejected: a reference is not sample evidence.

Example:

```
sample_id	input_type	assay	fastq_1	fastq_2	bam	vcf
S1	fastq	wgs	data/S1_R1.fastq.gz	data/S1_R2.fastq.gz		
S2	bam	wgs			data/S2.bam	
S3	vcf	wgs				calls/S3.vcf.gz
```

---

## 3. Configuration

A run configuration is a YAML file validated by the Pydantic schema in
`config/schema.py`. Start from `config/config.yaml` and set the required fields:

```yaml
schema_version: 1
project: PfVarX
profile: wgs

run:
  mode: wgs
  intent: research
  policy_status: PROVISIONAL_DEVELOPMENT
  run_id: RUN_001

reference:
  assembly_id: Pf3D7_01_v3          # pin a real assembly for VALIDATED
  fasta: reference/Pf3D7.fa
  sha256: <64-hex-digest>
  callable_mask_bed: reference/mask.bed

inputs:
  samplesheet: config/samples.tsv
  input_type: auto
  accept_existing_bam: true
  require_checksums: false

outputs:
  root: results
  retain_raw_evidence: true
  atomic_writes: true
  overwrite: false
```

Precedence is `defaults < profile < run config < CLI`. The `--profile` flag
merges a mode profile (e.g. `config/profiles/wgs.yaml`) beneath the run config,
so shared scientific defaults live in the profile and per-run overrides live in
the run config.

Three profiles ship with the platform:

| Profile | Mode | Key differences |
|---|---|---|
| `config/profiles/wgs.yaml` | Whole-genome | `wgs_min_independent_depth: 10`, duplicate marking on |
| `config/profiles/tbs_no_umi.yaml` | Targeted, no UMI | higher depth thresholds, no duplicate marking |
| `config/profiles/tbs_umi.yaml` | Targeted, UMI | `umi_design` required, `umi_min_family_size: 2` |

All profile values are `PROVISIONAL_DEVELOPMENT`. They are engineering starting
points, not validated cutoffs. See `docs/PARAMETERS.md`.

---

## 4. Validate and route

Always validate before running. This catches mislabelled files, missing
references, and samplesheet errors before any compute is spent:

```bash
pfvarx validate \
    --config config/config.yaml \
    --samplesheet config/samples.tsv \
    --objectives drug_resistance,diversity \
    --check-paths \
    --json
```

`validate` also resolves the selected objectives, checks their prerequisites and
prints the directories that will be created, so a missing metadata table or GFF
is caught here rather than halfway through a run. See
[`OBJECTIVES.md`](OBJECTIVES.md) for the objective list and dependency rules.

The command prints a routing plan. Each sample receives one of four routes:

- **`alignment`** — raw FASTQ will be trimmed and aligned;
- **`alignment_reuse`** — an existing BAM/CRAM is reused as-is (never realigned);
- **`variant_import`** — a VCF/BCF/gVCF is imported and never re-derived from reads;
- **`rejected`** — validation failed; the reason is reported.

Write the plan to a file for audit:

```bash
pfvarx validate --config config/config.yaml \
    --samplesheet config/samples.tsv --out results/routing_plan.json
```

---

## 5. Run the pipeline

### 5.1 Via the CLI (single sample or small batches)

```bash
pfvarx run --config config/config.yaml --samplesheet config/samples.tsv \
    --objectives drug_resistance,diversity --out results
```

`--objectives` overrides `objectives:` in the configuration file. Only the
stages required by the selected objectives (plus those implied by each sample's
input route) are executed, and only their directories are created. The run
follows a fixed ten-step sequence — read config, validate, determine start
stage, prepare inputs, build the plan, create directories, execute, record
logs/status, generate results, summarise — and writes `run_manifest.json` as
the single source of truth for reporting.

Exit codes: `0` completed; `1` partial, unavailable, or failed run; `2` CLI
usage or missing input/evidence file; `3` configuration or other PfVarX
application error. Missing `--evidence-dir` files and known output collisions
are checked before the CLI creates run outputs.

If evidence tables have already been produced (for example by the Snakemake
workflow, or a previous `pfvarx call`), add `--evidence-dir RESULTS` to skip
evidence production and execute only the selected objectives against the
existing tables. This is how the Snakemake tier hands off to the same
objective-driven executor, and it is safe to use interactively:

```bash
pfvarx run --config config/config.yaml --samplesheet config/samples.tsv \
    --objectives diversity --evidence-dir results --out results
```

### 5.2 Via Snakemake (recommended for batches)

```bash
snakemake --profile workflow/profiles/local
```

The Snakefile reads `config/config.yaml`, validates the samplesheet, plans
routes, and only schedules the stages each route requires. It also honours the
`objectives:` list: the workflow resolves your objectives, produces evidence
for every sample (FASTQ → align, BAM/CRAM → reuse, VCF → import), then runs
**only** the selected analyses through `pfvarx run --evidence-dir`. An
unselected analysis is never scheduled, and no directory is created for it.
Re-running resumes from the last completed step; `--rerun-incomplete` is set in
every profile.

Choose a profile for your environment:

```bash
snakemake --profile workflow/profiles/local    # workstation, 8 GB
snakemake --profile workflow/profiles/server   # shared server
snakemake --profile workflow/profiles/hpc      # Slurm
snakemake --profile workflow/profiles/cloud    # Kubernetes
```

---

## 6. Call a single BAM

For a targeted region or a quick check, call one BAM directly:

```bash
pfvarx call \
    --config config/config.yaml \
    --bam data/S1.bam \
    --sample-id S1 \
    --fasta reference/Pf3D7.fa \
    --chrom Pf3D7_01_v3 --start 1 --end 200000 \
    --out results/S1
```

Outputs:

- `results/S1/S1.calls.vcf` — VCF with `PVX_DECISION` and `PVX_POLICY_ID` INFO;
- `results/S1/S1.evidence.tsv` — the 46-column evidence table;
- a JSON summary printed to stdout (or `--json`).

---

## 7. Import a pre-called VCF

When a collaborator supplies calls, import them without touching reads:

```bash
pfvarx import \
    --config config/config.yaml \
    --vcf calls/S3.vcf.gz \
    --sample-id S3 \
    --out results/S3
```

Imported calls are labelled `model = "imported"` and
`calibration_profile = "IMPORTED"`. They can never be `ALT_PASS_VALIDATED`,
because no evidence was re-derived. They are restricted to
`ALT_DETECTED_EXPLORATORY`, `REFERENCE_SUPPORTED`, or `NO_CALL`.

---

## 8. Calibrate

The calibration utility estimates error-rate/over-dispersion parameters from
control data. Provide a TSV with `depth` and `alt_count` columns:

```bash
pfvarx calibrate --config config/config.yaml --controls controls.tsv --out calibration.json
```

The caller treats configured `statistics.error_rate` and
`statistics.beta_binomial_rho` as fixed null parameters; it does not propagate
their estimation uncertainty or model context-specific error rates. A calibration
output, pinned reference, or `run.policy_status: VALIDATED` does not by itself
establish assay validation or a verifiable limit-of-detection (LoD). The current
caller consumes no assay-validation/LoD provenance and cannot emit a validated
`ALT_PASS_VALIDATED` call. Treat results as exploratory/uncertain, not validated.

---

## 9. Benchmark and ablation

```bash
pfvarx benchmark --config config/config.yaml \
    --replicates 3 --coverage 60 --seed 11 \
    --out benchmarks/results
```

This runs a synthetic simulation benchmark and can report SNP/INDEL precision,
recall, F1, false positives, and Wilson confidence intervals against simulation
truth. Such tests/demo results show software behaviour under specified synthetic
conditions; they do not establish biological validation or empirical assay LoD.
Historical benchmark/ablation labels are not the caller's inference path: caller
inference is an exact one-sided ALT-count tail with a complete-scan Bonferroni
correction. Use `--no-ablation` for the benchmark without ablation. See
`benchmarks/README.md` for recorded results and limitations.

---

## 10. Annotate, downstream, relatedness, QC, report

```bash
# Functional annotation (internal engine; SnpEff selection is not implemented)
pfvarx annotate --config config/config.yaml \
    --evidence results/S1/S1.evidence.tsv --out results/S1/S1.annotated.tsv

# Diversity and mixed-infection summaries
pfvarx downstream --config config/config.yaml \
    --evidence results/S1/S1.evidence.tsv --out results/S1/downstream.json

# Pairwise relatedness and NJ tree across a results tree
pfvarx relatedness --config config/config.yaml --results results --out results/relatedness

# Unified QC across a results tree
pfvarx qc --config config/config.yaml --results results --out results/run_qc.json

# Reports
pfvarx report workbook --config config/config.yaml --results results --out results/report
pfvarx report pdf      --config config/config.yaml --results results --out results/report
pfvarx report dashboard --config config/config.yaml --results results
```

---

## 11. Demo (no external data required)

The fastest way to confirm the platform works end to end:

```bash
pfvarx demo --config config/profiles/wgs.yaml --out /tmp/pfvarx_demo --coverage 60 --seed 11
```

This generates a synthetic reference, a synthetic BAM with three variants at
allele fractions 0.50, 0.20, and 0.08, runs the caller, and writes a VCF plus an
evidence table. Inspect the counts and decisions rather than expecting specific
calls: the complete-scan Bonferroni bound can leave even simulated variants
`ALT_UNCERTAIN`, depending on coverage and the configured small-indel universe.
This demo verifies software flow, not call sensitivity or biological accuracy.

---

## 12. Output layout

```
results/
├── routing_plan.json
├── <sample>/
│   ├── <sample>.calls.vcf
│   ├── <sample>.evidence.tsv
│   ├── <sample>.annotated.tsv        # if annotation enabled
│   └── downstream.json               # if downstream enabled
├── relatedness/
│   ├── relatedness_matrix.tsv
│   └── tree.nwk
├── run_qc.json
└── report/
    ├── report.pdf
    ├── calls.xlsx
    └── figures/*.png
```

All writes are atomic (write-to-temp then rename) when
`outputs.atomic_writes: true`, so an interrupted run never leaves a truncated
artefact. Re-running with identical inputs and configuration produces
byte-identical evidence tables.

---

## 13. Interpreting the evidence table

The evidence table has 46 columns. The most important are:

| Column | Meaning |
|---|---|
| `ref_count`, `alt_count`, `other_count` | Base counts after filtering |
| `depth`, `informative_depth` | Raw and informative depth |
| `alt_fraction` | `alt_count / informative_depth` |
| `af_ci_low`, `af_ci_high` | Jeffreys credible interval |
| `model`, `model_version` | Count model used |
| `log_likelihood_h0`, `log_likelihood_h1` | Log-likelihoods |
| `lrt_statistic` | `2*(logL1 - logL0)` |
| `lrt_pvalue`, `lrt_pvalue_adjusted` | Exact one-sided null count-tail p-value and Bonferroni family-adjusted value (legacy field names; not LRT p-values or BH q-values) |
| `posterior_h1` | Exploratory Bayesian quantity when enabled; not assay-validated probability of biological truth |
| `decision` | One of the five decision categories |
| `filters` | Semicolon-separated filter names |
| `bin` | Minority-fraction bin |
| `sequencing_mode`, `umi_status` | Assay context |
| `policy_status`, `policy_id` | Provenance |

Decision categories:

- `NO_CALL` — a hard gate failed (low depth, outside callable region, reference mismatch);
- `REFERENCE_SUPPORTED` — evidence consistent with the reference allele;
- `ALT_DETECTED_EXPLORATORY` — alternate supported but not validated;
- `ALT_PASS_VALIDATED` — reserved category; the current caller cannot emit it without verifiable assay-validation/LoD provenance;
- `ALT_UNCERTAIN` — evidence present but flagged for review.

Minority bins: `f < 0.01` below threshold; `0.01–0.05`
`LOW_FRACTION_EXPLORATORY`; `0.05–0.20` `MINORITY_CANDIDATE`; `0.20–0.80`
`INTERMEDIATE_FRACTION`; `≥ 0.80` `PREDOMINANT_ALT`.

---

## 14. Interpreting the VCF

PfVarX emits a standards-compliant VCF. Two INFO fields carry platform
semantics:

- `PVX_DECISION` — the decision category;
- `PVX_POLICY_ID` — the configuration/parameter provenance identifier.

`QUAL` is not a validated probability or quality score. The caller currently
has no verifiable assay-validation/LoD provenance and cannot produce a validated
PASS; do not interpret a configured status or posterior as biological validation.
The raw likelihood ratio is not a posterior. See `docs/METHODS.md` for details.

---

## 15. Common workflows

**Reuse existing BAMs (no realignment).** Put BAM paths in the `bam` column and
set `input_type: bam`. The router selects `alignment_reuse`; the alignment
rules are never scheduled.

**Combine cohorts.** Run each sample into a shared `results/` root, then run
`pfvarx relatedness --results results` to compute the pairwise matrix and tree.

**Calibration and validation status.** `pfvarx calibrate` can estimate null
parameters from controls, and configuration supports a `VALIDATED` status, but
neither substitutes for documented assay validation and LoD provenance. The
current caller does not consume that provenance and will not emit
`ALT_PASS_VALIDATED`; do not label tests, demo outputs, or configured status as
biologically validated.

See `docs/TROUBLESHOOTING.md` for runtime issues and `docs/DEPLOYMENT.md` for
production deployment.
