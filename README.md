# PfVarX

**Research use only. No component of this pipeline is clinically validated.** Do not use PfVarX for diagnosis, treatment selection, or patient-facing decisions. See the [current audit report](docs/AUDIT_REPORT_2026-10-23.md) for verified software changes, test evidence, scientific limitations, and deployment scope.

PfVarX is a Pixi-managed pipeline for exploratory detection and characterization of SNPs and small insertions/deletions (up to 50 bp) in haploid *Plasmodium falciparum* whole-genome sequencing data. The caller accounts for mixed infections, but its defaults are provisional and uncalibrated. The pipeline supports three input routes:

| Input in the sample sheet | Processing route | Use this run command |
|---|---|---|
| Paired FASTQ/FASTQ.GZ | Trim, align, optionally mark duplicates, call, then analyze | Snakemake in the `workflow` Pixi environment |
| BAM or CRAM | Reuse the supplied alignment and call | `pfvarx run` in the `default` environment, or Snakemake |
| VCF, BCF, or gVCF | Validate/import existing calls; no realignment or PfVarX likelihood calling | `pfvarx run` in the `default` environment, or Snakemake |

The FASTQ workflow is for paired-end short reads using the fastp/BWA-MEM2 toolchain (Illumina-compatible input). Configured modes are WGS, targeted/no-UMI, and targeted/UMI; the UMI route requires an assay-matched read structure and is not scientifically validated. Single-end and long-read technologies (including Nanopore/PacBio workflows) are not established. Imported calls remain exploratory; CRAM, BCF, and gVCF were not independently end-to-end smoke-tested in the audit, so verify those routes against your actual files before use.

## Workflow at a glance

PfVarX has two separate stages. Stage 1 creates **input directories only** and is followed by manual placement and preparation of input files. Stage 2 validates the inputs, then creates the configured output, log, and work directories and runs the selected analysis.

```text
Stage 1: copy/edit config → `pfvarx init` → place and prepare input files
Stage 2: select objectives → validate → run pipeline → inspect results/logs
```

`pfvarx init` does not create `results`, `logs`, or `work` directories. The pipeline creates those only after Stage 2 preflight succeeds.

## 1. Requirements and Pixi installation

Use [Pixi](https://pixi.sh) for dependency installation and commands. Do not create or activate a separate Conda environment. The committed `pixi.lock` covers Linux x86-64, Linux ARM64, and Intel macOS. The audit verified native Linux ARM64 installation and execution; the other declared platforms were not natively tested in that audit. Windows is not a declared platform; use WSL2 if appropriate for your environment.

Pixi manages the project's Python 3.11 runtime and resolves packages from the declared conda-forge/Bioconda channels plus the Python package index where specified.

The local Snakemake profile defaults to four cores and one job and reserves 8 GB for heavy rules. Adjust the profile/configuration for your system and available memory before running a larger dataset.

Install Pixi and clone the project, then enter the project directory:

```bash
curl -fsSL https://pixi.sh/install.sh | bash
git clone https://github.com/hadish-bekuretsion/PfVarX_Pipeline.git
cd PfVarX_Pipeline
pixi install --locked -e default
pixi run -e default status
pixi run -e default test
```

For paired-FASTQ processing and reference indexing, install the workflow environment too:

```bash
pixi install --locked -e workflow
pixi run -e workflow snakemake --version
```

The `workflow` environment contains Snakemake and the declared bioinformatics executables, including `fastp`, `bwa-mem2`, `samtools`, `picard`, and `fgbio`. SnpEff is installed as a tool but is not integrated into PfVarX's per-call annotation path. PfVarX's configurable annotation engine accepts `internal` or `auto` only. Optional DEploid, hmmIBD, and IQ-TREE wrappers are not provisioned by this environment.

## 2. Stage 1 — create input directories and prepare files

### 2.1 Copy and edit a configuration file

Copy the example configuration, then edit it for your data and run. Run subsequent commands from the repository directory so relative paths in the configuration resolve predictably.

```bash
cp config/config.yaml config/my-run.yaml
```

At minimum, edit these locations and identifiers in `config/my-run.yaml`:

```yaml
reference:
  assembly_id: Pf3D7_v68
  fasta: data/reference/Pf3D7.fa
  sha256: REPLACE_WITH_64_HEX_SHA256

inputs:
  root: data
  raw_dir: data/raw
  prepared_dir: data/prepared
  metadata_dir: null             # set a path only when you will use a separate metadata file
  samplesheet: data/prepared/samples.tsv

outputs:
  root: results/my-run           # final results; choose a distinct directory for each run
  logs_dir: null                 # defaults to <outputs.root>/logs
  work_dir: null                 # defaults to <outputs.root>/intermediate
```

`reference.sha256` must be the actual 64-character SHA-256 digest of the supplied FASTA. You can initially leave the example's `REQUIRED_SHA256` placeholder while creating folders, but `validate`, `run`, and direct `call`/`import`/`annotate` commands reject placeholders and verify the digest against the exact FASTA before reference-dependent analysis. A mismatch fails before output directories are created. The `inputs.*_dir` fields tell `init` where to create input folders. They do not rewrite file paths in the sample sheet; each data-file path must still be given in the sample sheet or its appropriate configuration field.

### 2.2 Create only the input directories

```bash
pixi run -e default pfvarx init --config config/my-run.yaml
```

The command creates the configured raw and prepared directories, any explicitly configured metadata directory, and parent directories for declared reference, sample-sheet, optional BED, GFF3, and metadata paths. Existing directories are left intact. If an input location overlaps the configured results, logs, or work location, initialization stops instead of creating it. No sequencing files are copied, downloaded, indexed, or modified by `init`.

With the example paths, the Stage 1 layout is:

```text
data/
├── raw/                       # user places FASTQ, BAM/CRAM, or VCF/BCF/gVCF files here
├── prepared/
│   └── samples.tsv            # user creates this tab-separated sample sheet
├── reference/
│   ├── Pf3D7.fa               # user supplies the matching reference FASTA
│   ├── Pf3D7.fa.fai           # user generates with samtools faidx
│   └── bwa-mem2 index files   # user generates for the FASTQ route only
└── metadata/                  # optional; configure inputs.metadata_dir to create it
```

If `annotation.gff`, `reference.callable_mask_bed`, `preprocessing.primer_bed`, `downstream.metadata_table`, or a custom resistance catalogue is configured at another path, `init` creates the parent folder for that declared file as well. It does not create any file contents.

### 2.3 Prepare the reference and required resources

PfVarX does not bundle a reference genome. Supply a FASTA for the exact assembly used to produce any existing BAM/CRAM or VCF/BCF inputs. Keep the assembly identifier, reference sequence, and annotation consistent.

Generate and record the FASTA checksum and index:

```bash
sha256sum data/reference/Pf3D7.fa
pixi run -e workflow samtools faidx data/reference/Pf3D7.fa
```

Copy the resulting SHA-256 digest into `reference.sha256`. The `.fai` index is needed for reference access. For paired FASTQ data, build the BWA-MEM2 index once in the same reference location:

```bash
pixi run -e workflow bwa-mem2 index data/reference/Pf3D7.fa
```

The current workflow does **not** generate or consume a Picard/GATK `.dict` file, so a sequence dictionary is not a PfVarX input requirement. Create one only for another tool that independently requires it.

Additional files are required only for the relevant configuration/objective:

| File/resource | When needed | Format and location |
|---|---|---|
| GFF3 gene model | `drug_resistance` objective | GFF3 matching the exact FASTA assembly; set `annotation.gff` |
| Callable-region BED | Optional | BED matching the reference contigs; set `reference.callable_mask_bed` |
| Primer BED | Optional targeted preprocessing | BED appropriate to the assay; set `preprocessing.primer_bed` |
| Sample metadata table | `population_structure` or `spatial_temporal` | TSV/CSV keyed by `sample_id`; set `downstream.metadata_table` or pass `--metadata` |
| Resistance catalogue | Optional override | PfVarX-compatible catalogue YAML; set `downstream.resistance_catalogue`; otherwise the packaged catalogue is used |

### 2.4 Prepare the sample sheet

Create a tab-separated file at the configured `inputs.samplesheet` path. It must contain one row per sample and at least `sample_id` and `input_type`. File paths in the sheet are resolved relative to the **sample-sheet directory** unless they are absolute (or `--base-dir` is supplied to the CLI).

```tsv
sample_id	input_type	assay	fastq_1	fastq_2	bam	vcf	metadata_json
PF001	fastq	wgs	../raw/PF001_R1.fastq.gz	../raw/PF001_R2.fastq.gz			
PF002	bam	wgs			../raw/PF002.bam		
PF003	vcf	wgs				../raw/PF003.vcf.gz	{"country":"Ethiopia","collection_date":"2023-11-04"}
```

| Column | Required | Meaning |
|---|---|---|
| `sample_id` | Yes | Unique sample label used in output names and directories. |
| `input_type` | Yes | One of `fastq`, `bam`, `cram`, `vcf`, `bcf`, or `gvcf`; `auto` is also accepted. |
| `fastq_1`, `fastq_2` | FASTQ rows | Paired-end reads. Both files are required for this route. |
| `bam` | BAM/CRAM rows | Existing alignment. The reference must match the alignment header. |
| `vcf` | VCF/BCF/gVCF rows | Pre-called variants; this route imports, but does not realign or recalculate PfVarX evidence. |
| `assay` | Optional | Informational assay label, commonly `wgs` or `tbs`; `run.mode` controls PfVarX's configured processing mode. |
| `metadata_json` | Optional | Inline per-sample JSON metadata; useful for small tables. |

Paths may be absolute or relative to the sheet. Use consistent paired-read names such as `<sample>_R1.fastq.gz` and `<sample>_R2.fastq.gz` to make the sheet easy to review; the parser does not require that particular naming pattern. Sample IDs should be unique and should not contain path separators.

For a separate metadata TSV, use one row per sample and include the columns required by the selected objective. For example:

```tsv
sample_id	country	collection_date	population
PF001	Ethiopia	2023-11-04	Amhara
PF002	Ethiopia	2023-12-19	Oromia
```

`population_structure` requires metadata including `country`; `spatial_temporal` requires `country` and `collection_date`. `relatedness` and `phylogenetics` require at least two samples. The remaining objectives are evaluated per sample, subject to their other prerequisites. Run `pfvarx objectives` for the current registry.

## 3. Configuration reference

The YAML loader merges built-in defaults, an optional mode profile, the run configuration, and then command-line overrides. Unknown schema keys are rejected. Paths are interpreted from the directory where the command runs; sample-sheet data paths are the exception and are relative to the sheet directory by default.

| Section / field | Purpose |
|---|---|
| `profile` | Optional loader directive: `wgs`, `tbs_umi`, or `tbs_no_umi`; profile values are starting points, not assay validation. |
| `run.mode` | Processing mode: `wgs`, `tbs_no_umi`, or `tbs_umi`. UMI mode requires a verified `preprocessing.umi_design`. |
| `run.intent`, `run.policy_status` | Record intended use and validation status. Keep `policy_status: PROVISIONAL_DEVELOPMENT` unless the assay and thresholds have been independently validated. |
| `run.run_id`, `run.description` | Optional run label and free-text description. |
| `reference.assembly_id`, `reference.fasta`, `reference.sha256` | Identify, locate, and checksum the reference assembly. |
| `reference.callable_mask_bed` | Optional callable-region mask. |
| `inputs.root`, `inputs.raw_dir`, `inputs.prepared_dir`, `inputs.metadata_dir` | Locations used by `pfvarx init` to create input folders. These are not automatic file-discovery rules. |
| `inputs.samplesheet` | Configured sample-sheet path; the run CLI also requires `--samplesheet`. |
| `inputs.input_type`, `accept_existing_bam`, `require_checksums` | Input-routing and reuse preferences; `require_checksums` is not currently enforced for read/alignment/variant files. The row's `input_type` and explicit file paths determine routing. |
| `preprocessing.*` | FASTQ quality trimming, adapter detection, primer, mapping/base-quality, overlap, duplicate-marking, and UMI fields; not all are wired to the workflow (see the complete field table below). |
| `callability.*` | Mode-specific depth/callability values, some of which are used as the caller's depth gate; extended review thresholds are not all wired (see below). Defaults are provisional. |
| `candidate.*` | Minimum alternate observations/fraction, maximum small-indel size, homopolymer review length, and strand-bias threshold. |
| `statistics.*` | Count model, empirical error rate/profile, over-dispersion, alternative-hypothesis fraction, null calibration, bootstrap count, legacy `*_fdr_q` decision cutoffs, and optional Bayesian settings. The caller's implemented correction and limitations are described in the [methods](docs/METHODS.md) and [audit report](docs/AUDIT_REPORT_2026-10-23.md). |
| `reporting.*` | Fraction bins and accepted PASS/uncertainty preferences; the current policy always fails closed for validated PASS and emits uncertainty regardless of dormant flags. |
| `benchmark.*` | Declared target precision/recall and confidence values; the current synthetic benchmark runner does not consume these targets and they are not evidence of achieved performance. |
| `resources.*` | Threads per sample and per-heavy-rule memory are used where workflow rules request them; `concurrent_samples` does not override the shipped local profile's fixed job count. |
| `outputs.root` | Final results root. Choose a unique path for each run. |
| `outputs.logs_dir`, `outputs.work_dir` | Optional independent log and temporary/intermediate roots. When null, defaults are `<outputs.root>/logs` and `<outputs.root>/intermediate`. |
| `outputs.retain_raw_evidence`, `atomic_writes`, `overwrite` | `overwrite` controls known file-collision refusal; the retention and atomic-write preferences are accepted but do not currently control runtime behavior. `overwrite: true` does not clean the output tree. |
| `annotation.enabled`, `annotation.engine`, `annotation.gff` | Annotation switch, implemented engine (`internal` or `auto`), and matching GFF3. Drug-resistance analysis requires the GFF3. |
| `downstream.enabled`, `resistance_catalogue`, `metadata_table` | Legacy objective list (prefer top-level `objectives`), optional catalogue override, and metadata table. |
| `objectives` | List of objectives to execute. This is the preferred selection interface. |

### Field-by-field reference

The table below covers the fields accepted by the current strict configuration schema. “Accepted” means the YAML validator accepts and checks the value; it does **not** imply that every field is consumed by every analysis route. In particular, dormant settings are identified so users do not assume that changing them alters analysis behavior. Numeric defaults and accepted ranges are available from `pfvarx parameters` and [`docs/PARAMETERS.md`](docs/PARAMETERS.md).

| Key | Purpose and implementation status |
|---|---|
| `schema_version` | Configuration schema version; currently `1`. |
| `project` | Human-readable project label included in run metadata. |
| `profile` | Loader directive selecting the `wgs`, `tbs_umi`, or `tbs_no_umi` profile before run-specific overrides. |
| `run.mode` | Required assay-processing mode: `wgs`, `tbs_no_umi`, or `tbs_umi`. UMI mode requires a `preprocessing.umi_design`. |
| `run.intent` | `research`, `surveillance`, or `clinical_research`; the latter requires a validated status in schema but does not establish clinical validity. |
| `run.policy_status` | `PROVISIONAL_DEVELOPMENT` or `VALIDATED`. Defaults to provisional; a label alone never validates thresholds or enables a validated call. |
| `run.run_id` | Optional 1–64 character run label, restricted to letters, digits, dot, underscore, and hyphen; begins with a letter/digit. |
| `run.description` | Optional free-text run description. |
| `reference.assembly_id` | Assembly identifier corresponding to the supplied FASTA and any BAM/VCF inputs. |
| `reference.fasta` | Reference FASTA path. The user supplies the file and matching indexes. |
| `reference.sha256` | Required 64-hex SHA-256 for the exact FASTA. Stage 1 accepts a `REQUIRED_*` setup placeholder; `validate`, `run`, and direct `call`, `import`, and `annotate` recompute and compare the actual file digest, rejecting placeholders and mismatches before reference-dependent analysis. |
| `reference.callable_mask_bed` | Optional BED mask for callable regions. |
| `reference.snpeff_db`, `reference.snpeff_jar` | Legacy SnpEff metadata fields; SnpEff is not integrated into PfVarX annotation and these fields do not activate it. |
| `inputs.samplesheet` | Required path to the TSV/CSV sample sheet; individual file paths are resolved relative to the sheet by default. |
| `inputs.root` | Base input directory for Stage 1 defaults. |
| `inputs.raw_dir` | Directory Stage 1 creates for raw read/alignment/variant inputs; users place files there manually. |
| `inputs.prepared_dir` | Directory Stage 1 creates for the prepared sample sheet and other manually prepared files. |
| `inputs.metadata_dir` | Optional metadata directory; Stage 1 creates it only when explicitly configured. |
| `inputs.input_type` | Default input-type hint: `auto`, `fastq`, `bam`, `cram`, `vcf`, `bcf`, or `gvcf`; sample-row paths and type determine the route. |
| `inputs.accept_existing_bam` | Whether existing BAM/CRAM input may be reused; controls input routing. |
| `inputs.require_checksums` | Accepted by the schema but not currently enforced for sample-sheet read/alignment/variant files. Reference FASTA checksum comparison is enforced separately through `reference.sha256`. |
| `preprocessing.fastp_min_phred` | Minimum Phred quality used by the FASTQ trimming rule. |
| `preprocessing.fastp_min_length` | Minimum retained read length passed to fastp. |
| `preprocessing.fastp_max_n_fraction` | Maximum-N fraction is declared in the schema/registry but is not currently forwarded to the fastp workflow rule. |
| `preprocessing.adapter_autodetect` | Enables fastp paired-end adapter autodetection in the workflow. |
| `preprocessing.primer_bed` | Optional primer BED path; Stage 1 creates its parent, but the current workflow does not use it for primer trimming/clipping. |
| `preprocessing.mapping_quality_min` | Minimum mapping quality for read evidence. |
| `preprocessing.base_quality_min` | Minimum base quality for read evidence. |
| `preprocessing.count_overlapping_mates_once` | Prevents overlapping paired-end mate observations from being counted twice where implemented by the evidence extractor. |
| `preprocessing.wgs_mark_duplicates` | Enables Picard duplicate marking in the WGS preprocessing route. |
| `preprocessing.duplicate_exclusion_from_evidence` | Excludes marked duplicates from read-level evidence when enabled. |
| `preprocessing.umi_design` | UMI read-structure description used by the fgbio branch; it must match a verified assay layout. |
| `preprocessing.umi_min_family_size` | Minimum fgbio family size for UMI consensus processing. |
| `callability.wgs_min_independent_depth` | WGS minimum depth threshold returned for the WGS mode. |
| `callability.targeted_min_independent_depth` | Targeted, non-UMI minimum depth threshold. |
| `callability.umi_min_consensus_molecules_majority` | UMI-mode majority-depth threshold. |
| `callability.umi_min_consensus_molecules_exploratory` | Accepted setting for an exploratory UMI molecule threshold; not currently used as a distinct decision cutoff. |
| `callability.wgs_high_depth_multiplier_for_review` | Declared high-depth review parameter; not currently wired into the caller's decision path. |
| `callability.wgs_max_depth_floor` | Declared maximum-depth parameter; not currently wired into the caller's decision path. |
| `callability.wgs_min_callable_fraction` | Declared minimum callable-fraction parameter; not currently wired into the caller's decision path. |
| `candidate.snp_min_alt_observations` | Minimum alternate observations required to retain SNP candidates. |
| `candidate.indel_min_alt_observations` | Minimum alternate observations required to retain small-indel candidates. |
| `candidate.exploratory_min_alt_fraction` | Minimum alternate fraction for candidate discovery; this is a discovery threshold, not an LoD. |
| `candidate.max_small_indel_bp` | Maximum supported small-indel length; schema permits 1–50 bp. |
| `candidate.homopolymer_review_length` | Homopolymer length at/above which indel evidence is flagged for review. |
| `candidate.strand_bias_alpha` | Threshold for strand-bias review. |
| `statistics.count_model` | `binomial`, `beta_binomial`, or `auto`; selects the implemented count-tail null model. |
| `statistics.empirical_error_profile` | Accepted profile-path field, but the runtime does not currently load this path; use only explicit supported error parameters and do not infer calibration from this field. |
| `statistics.error_rate` | Null alternate-read probability used by the count model; must come from suitable data for scientific interpretation. |
| `statistics.beta_binomial_rho` | Optional beta-binomial over-dispersion parameter; its estimation uncertainty is not propagated. |
| `statistics.min_alt_fraction_h1` | Lower alternate-fraction boundary for the composite alternative. |
| `statistics.lrt_null_calibration` | Selects whether the caller requests auxiliary bootstrap diagnostics; final per-site calling currently uses exact count-tail p-values, not the legacy LRT field. |
| `statistics.bootstrap_replicates` | Accepted bootstrap-count setting, but the current caller does not pass this value into the bootstrap implementation; it does not change the exact count-tail decision. |
| `statistics.exploratory_fdr_q` | Exploratory adjusted-p-value cutoff; legacy “FDR” name notwithstanding, the caller applies Bonferroni family-wise adjustment, not BH FDR. |
| `statistics.sensitivity_fdr_q` | Sensitivity adjusted-p-value cutoff; also a Bonferroni-adjusted p-value threshold, not a BH FDR level. |
| `statistics.bayesian_enabled` | Enables exploratory posterior output when a threshold is explicitly configured. |
| `statistics.bayesian_prior_h1` | Prior probability used by optional exploratory Bayesian output. |
| `statistics.posterior_threshold_if_validated` | Explicit posterior threshold required when Bayesian mode is enabled; does not itself validate a call. |
| `statistics.low_depth_fallback` | Accepted as `conservative` or `fail`; it is not currently a separate active caller branch. |
| `reporting.minority_candidate_fraction` | Lower boundary of the reporting fraction categories. |
| `reporting.intermediate_fraction` | Boundary used for the intermediate allele-fraction category. |
| `reporting.predominant_fraction` | Boundary for the predominant-alt category; boundaries must be non-decreasing. |
| `reporting.pass_requires_assay_validated_lod` | Schema setting for a validated-LoD PASS gate; the current policy fails closed for validated PASS regardless of this flag. |
| `reporting.uncertainty_always_visible` | Accepted reporting preference; uncertainty is currently emitted by the implementation and this flag does not suppress it. |
| `benchmark.confidence_level` | Declared benchmark confidence target; the current synthetic benchmark runner does not consume this setting. |
| `benchmark.plan_snp_precision`, `benchmark.plan_snp_recall` | Declared SNP performance targets; not consumed by the current runner and not achieved-performance claims. |
| `benchmark.plan_indel_precision`, `benchmark.plan_indel_recall` | Declared indel performance targets; not consumed by the current runner and not achieved-performance claims. |
| `benchmark.plan_lod_sensitivity` | Declared LoD sensitivity target; not consumed by the current runner and not an empirically established LoD. |
| `resources.threads_per_sample` | Threads passed to processing tools where the workflow uses this setting. |
| `resources.concurrent_samples` | Accepted concurrency preference; the shipped local Snakemake profile currently fixes its own job count. |
| `resources.mem_mb_per_heavy_rule` | Memory request used for heavy Snakemake rules. |
| `outputs.root` | Final-results directory. |
| `outputs.logs_dir` | Optional separate log directory; defaults to `<outputs.root>/logs`. |
| `outputs.work_dir` | Optional separate temporary/intermediate directory; defaults to `<outputs.root>/intermediate`. |
| `outputs.retain_raw_evidence` | Accepted preference, but currently does not control evidence retention. |
| `outputs.atomic_writes` | Accepted preference, but currently does not switch the writer's atomic-write behavior. |
| `outputs.overwrite` | Controls known file-collision refusal. `true` does not clean existing directories. |
| `annotation.enabled` | Master annotation preference. |
| `annotation.engine` | Implemented choices are `internal` or `auto`; SnpEff is not a supported engine. |
| `annotation.snpeff_db` | Legacy field; does not activate SnpEff annotation. |
| `annotation.gff` | GFF3 gene model used for annotation; required for the drug-resistance objective. |
| `downstream.enabled` | Legacy objective list. Prefer top-level `objectives`; schema requires both selections to agree if both are supplied. |
| `downstream.resistance_catalogue` | Optional custom resistance catalogue; otherwise the packaged catalogue is used. |
| `downstream.metadata_table` | Optional sample metadata file for objectives that require it. |
| `objectives` | Preferred list of selected objectives; only requested objectives are executed. |

The full parameter registry, defaults, accepted ranges, and provisional scientific rationale are in [`docs/PARAMETERS.md`](docs/PARAMETERS.md). To print the actual registry from the installed code:

```bash
pixi run -e default pfvarx parameters
pixi run -e default pfvarx parameters --mode wgs
pixi run -e default pfvarx objectives
```

## 4. Stage 2 — configure objectives, validate, and run

### 4.1 Select objectives

Set only the analyses you want in `config/my-run.yaml`:

```yaml
objectives:
  - diversity
  - mixed_infection
```

Available objectives are `drug_resistance`, `mixed_infection`, `coi`, `diversity`, `population_structure`, `relatedness`, `phylogenetics`, and `spatial_temporal`. A CLI `--objectives` list overrides the configuration. Missing prerequisites are reported; `unavailable` does not mean that no biological signal was found.

For example, use one of these objective lists in the config (replace the prior `objectives` list rather than appending every analysis):

Per-sample summaries:

```yaml
objectives: [diversity, mixed_infection, coi]
```

Resistance review (also configure `annotation.gff` for the matching reference):

```yaml
objectives: [drug_resistance]
```

Cohort analyses (provide at least two samples):

```yaml
objectives: [relatedness, phylogenetics]
```

Geographic/time summaries (configure `downstream.metadata_table` or include
suitable metadata in each row's `metadata_json`):

```yaml
objectives: [population_structure, spatial_temporal]
```

Only the objective list should be used in a given run configuration. The examples illustrate prerequisites and do not imply that any objective is scientifically validated.

### 4.2 Preflight validation (no run directories created)

After the reference, sheet, sample files, and objective-specific inputs are present, run validation without `--out` to avoid writing a plan file:

```bash
pixi run -e default pfvarx validate \
  --config config/my-run.yaml \
  --samplesheet data/prepared/samples.tsv \
  --objectives diversity,mixed_infection
```

This checks the sample sheet, routes the data, checks paths and objective prerequisites, and previews the Stage 2 directory plan. It does not create the results, log, or work directories. Supplying `--out path/to/plan.json` intentionally writes a validation-plan file at that path.

### 4.3 Run paired FASTQ data with Snakemake

FASTQ alignment is performed by Snakemake; the in-process `pfvarx run` command does not align reads. First preview the Snakemake plan, then execute it:

```bash
pixi run -e workflow snakemake \
  --snakefile workflow/Snakefile \
  --configfile config/my-run.yaml \
  --profile workflow/profiles/local \
  --dry-run

pixi run -e workflow snakemake \
  --snakefile workflow/Snakefile \
  --configfile config/my-run.yaml \
  --profile workflow/profiles/local
```

The local profile sets its own core/job defaults. Snakemake creates the configured Stage 2 results, logs, and work directories as rules execute. For non-local execution, use only a deployment profile that has been configured and tested in your environment; the Slurm/Kubernetes profiles are not covered by the current audit evidence.

### 4.4 Run BAM/CRAM or pre-called VCF/BCF/gVCF

The default Pixi environment supports the in-process caller for BAM/CRAM and import for pre-called variant files. For example:

```bash
pixi run -e default pfvarx run \
  --config config/my-run.yaml \
  --samplesheet data/prepared/samples.tsv \
  --objectives diversity,relatedness
```

Use at least two samples for `relatedness` or `phylogenetics`. You may use the full Snakemake workflow for these routes too. `--out` overrides `outputs.root`; if an explicit `outputs.logs_dir` or `outputs.work_dir` is set, those configured destinations remain in effect. By default, run outputs are written under the configured `outputs.root`, with logs and intermediate files in its `logs/` and `intermediate/` subdirectories.

## 5. Expected Stage 2 output layout

The final-results root is `outputs.root` (or a CLI `--out` override). The actual per-sample evidence location differs by route: Snakemake writes evidence under `<results>/<sample>/`; the direct BAM/VCF CLI writes it under `<results>/variants/<sample>/`.

```text
<outputs.root>/
├── run_manifest.json             # objective status, inputs, stages and artifacts
├── run_summary.json              # direct CLI run summary (when produced)
├── routing_plan.json             # input routing decision
├── validation/                   # Snakemake sample-sheet validation record
├── config/effective_config.json  # resolved configuration snapshot (CLI run)
├── <sample>/                     # Snakemake evidence and VCF outputs
│   ├── <sample>.evidence.tsv
│   └── <sample>.pfvarx.vcf
├── variants/<sample>/             # direct CLI BAM/CRAM or VCF-import outputs
├── objectives/<objective>/        # selected objective outputs only
├── reports/                       # generated PDF/XLSX exports when requested
└── run_qc.json                    # unified QC summary when produced

<outputs.logs_dir>/                # defaults to <outputs.root>/logs
├── run.log
├── validate.log, report.log, ...  # workflow stage logs
└── <sample>/                      # per-sample trimming, alignment and call logs

<outputs.work_dir>/                # defaults to <outputs.root>/intermediate
└── <sample>/
    ├── trimmed/                   # FASTQ route: trimmed reads + fastp JSON/HTML reports
    └── aligned/                   # FASTQ route: BAM/SAM and index intermediates
```

Some outputs appear only when their corresponding stage ran. If you configure separate log/work paths, they are created at those locations instead of under the results root. A separate work path can be on scratch storage; preserve it until the workflow finishes and any needed intermediate files have been archived.

### Where to look first

- `run_manifest.json`: authoritative objective statuses, artifacts, warnings, and errors.
- `run_summary.json`: sample routes and run-level summary for direct CLI runs.
- `logs_dir`: workflow stdout/stderr logs, including errors from external tools.
- `<sample>.evidence.tsv` and `<sample>.pfvarx.vcf`: call evidence and emitted variants; review filters and uncertainty fields, not only the VCF.
- `objectives/<objective>/`: results for the objectives you selected.
- `reports/`: PDF/XLSX exports generated by the workflow/report command.

`completed` means the software stage completed; it is not scientific validation. `failed` means inspect the artifact's error and stage log. `unavailable` means a prerequisite was absent and must not be interpreted as a negative biological result.

## 6. Troubleshooting

| Problem | Likely cause | Action |
|---|---|---|
| `pfvarx init` says input/output paths overlap | An input folder was configured at or below `outputs.root`, `outputs.logs_dir`, or `outputs.work_dir` | Put raw/reference/prepared inputs in separate paths, then rerun `init`. |
| Config validation reports an invalid checksum | `reference.sha256` is missing, null, or not 64 hexadecimal characters | Recompute with `sha256sum` and paste the exact digest into the YAML. |
| Reference or `.fai` is missing | FASTA path is incorrect or index was not generated | Check `reference.fasta`; run `samtools faidx` in the workflow environment. |
| BWA reports missing index files | FASTQ reference was not indexed | Run `bwa-mem2 index` on the configured FASTA before the workflow. |
| FASTQ run does not align | The in-process `pfvarx run` command is not the FASTQ aligner | Use Snakemake in the `workflow` Pixi environment. |
| Sample files are reported missing | Relative paths are relative to the sample-sheet directory | Correct the path in the TSV or use an absolute path / `--base-dir`. |
| FASTQ pair is incomplete | Only one mate path is set or one file is absent | Supply readable paired `fastq_1` and `fastq_2` files. |
| Drug-resistance objective is unavailable | Matching GFF3 is absent, or the catalogue cannot be found | Set `annotation.gff` to a GFF3 matching the reference; check any catalogue override. |
| Metadata prerequisite fails | Required columns/sample rows are missing | Provide the `sample_id`-keyed metadata table with objective-required fields. |
| Multi-sample objective refuses to run | Fewer than two usable samples | Add samples or remove `relatedness`/`phylogenetics`. |
| Existing result files block a run | `outputs.overwrite` defaults to false | Choose a new unique `outputs.root`; only set overwrite true when replacing files is intentional. It does not clean old outputs. |
| Snakemake reports an external executable missing | Workflow environment was not installed or not activated through Pixi | Run `pixi install --locked -e workflow` and invoke with `pixi run -e workflow`. |
| Output says `ALT_UNCERTAIN`, `NO_CALL`, or `UNAVAILABLE` | The result is uncertain or the objective had no valid prerequisite/evidence | Read the evidence filters, manifest message, and logs; do not interpret these statuses as reference/no signal. |

For detailed exception meanings, see [`docs/TROUBLESHOOTING.md`](docs/TROUBLESHOOTING.md). For deployment constraints and tested environments, see [`docs/DEPLOYMENT.md`](docs/DEPLOYMENT.md).

## 7. Scientific scope and limitations

PfVarX's implemented scope is exploratory SNP and small-indel evidence for haploid *P. falciparum*. Structural variants, copy-number variation, diploid genotyping, and clinical interpretation are outside scope. The caller uses a conservative multiple-testing adjustment over the defined candidate universe; the adjustment is per sample, not cohort-wide, and can reduce sensitivity. Defaults have not been validated against an assay truth set. A `VALIDATED` config label alone does not establish validation or a limit of detection.

The `coi` objective does not claim a validated numeric complexity-of-infection estimate. Resistance matches are review candidates and depend on a matching reference annotation. Synthetic tests and end-to-end smoke runs verify software dataflow only; they are not evidence of sensitivity, specificity, calibration, or clinical validity.

For methods, see [`docs/METHODS.md`](docs/METHODS.md), [`docs/METHODOLOGY_CHAPTER.md`](docs/METHODOLOGY_CHAPTER.md), [`docs/OBJECTIVES.md`](docs/OBJECTIVES.md), and [`docs/VALIDATION.md`](docs/VALIDATION.md). The current software audit and explicit unresolved scientific risks are in [`docs/AUDIT_REPORT_2026-10-23.md`](docs/AUDIT_REPORT_2026-10-23.md).

## 8. Common commands

Run commands through Pixi from the project directory:

```bash
# Stage 1: create input directories only
pixi run -e default pfvarx init --config config/my-run.yaml

# List available objectives and parameter registry
pixi run -e default pfvarx objectives
pixi run -e default pfvarx parameters

# Stage 2: validate and run BAM/CRAM or VCF/BCF/gVCF
pixi run -e default pfvarx validate --config config/my-run.yaml --samplesheet data/prepared/samples.tsv
pixi run -e default pfvarx run --config config/my-run.yaml --samplesheet data/prepared/samples.tsv

# FASTQ Stage 2: run the Snakemake workflow
pixi run -e workflow snakemake --snakefile workflow/Snakefile --configfile config/my-run.yaml --profile workflow/profiles/local
```

## License

MIT. See [`LICENSE`](LICENSE).
