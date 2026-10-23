# PfVarX Architecture

PfVarX is a modular, configuration-driven platform for calling and interpreting
single-nucleotide polymorphisms (SNPs) and small insertions/deletions (INDELs,
1–50 bp) in *Plasmodium falciparum* from Illumina short-read data. This document
describes the layered architecture, the module map, the data-flow contracts
between layers, and the design invariants that the implementation enforces.

The architecture is deliberately organised as **five analytical layers (0–4)**
sitting on top of a **shared infrastructure core**. Each layer has a single
responsibility, a typed boundary, and an explicit failure mode. No layer reaches
across another; all cross-layer communication happens through immutable
dataclasses or on-disk artefacts with a documented schema.

---

## 1. Design principles

The implementation is governed by five invariants. They are not aspirations —
each is enforced in code and covered by tests.

1. **Evidence before inference.** Every call is derived from an explicit,
   persisted evidence record (reference/alternate/uninformative counts, depth,
   allele fraction, strand and position metrics, model identity and version).
   The evidence table is a first-class output, not a debugging by-product.
2. **Statistical honesty.** A raw likelihood-ratio statistic is never relabelled
   as a posterior probability, and an uncalibrated metric is never encoded as a
   VCF `QUAL` value. The writer requires `VALIDATED` status, a posterior, and an
   `ALT_PASS_VALIDATED` decision for QUAL; the current caller cannot produce that
   validated PASS without assay-validation/LoD provenance, so normal caller
   output leaves QUAL missing. This is tested by
   `tests/unit/test_writer.py::test_vcf_never_encodes_lrt_as_qual`.
3. **Fail closed.** External-tool adapters raise typed errors when a required
   executable, database, or assumption is absent. SnpEff selection is rejected
   by the configuration schema and per-call annotation uses only the internal
   engine; guarded wrappers for other optional tools do not imply that those
   analyses were run or scientifically validated.
4. **Provenance is mandatory.** Every run carries a configuration hash, a policy
   identifier, a parameter-registry hash, and per-artefact checksums. Two runs
   with identical inputs and configuration are byte-identical (verified by
   `tests/integration/test_end_to_end.py::test_rerun_is_idempotent`).
5. **Research-use honesty.** The platform advertises its own validation status.
   `pfvarx status` prints a banner of the form
   `RESEARCH USE ONLY — validation status: <STATUS>`; the status is computed as
   the weakest scientific status across all components, so it cannot be
   overstated.

---

## 2. Layer model

| Layer | Name | Responsibility | Principal modules |
|---|---|---|---|
| 0 | Input validation & routing | Detect file types, verify checksums, select the processing route, never realign a valid BAM/CRAM | `io/validators.py`, `io/routing.py`, `io/checksums.py` |
| 1 | Preprocessing & read QC | Trim, align (BWA-MEM2), sort, mark duplicates, UMI consensus | `preprocess/tools.py`, `preprocess/pipeline.py` |
| 2 | Statistical calling & variant QC | Extract evidence, fit count models, apply decision policy, emit VCF + evidence table | `calling/evidence.py`, `models/count_models.py`, `calling/policy.py`, `calling/writer.py` |
| 2.5 | Annotation | Consequence prediction and resistance-catalogue matching | `annotation/*` |
| 3 | Downstream analyses | Diversity, mixed infection, relatedness, phylogenetics, spatial-temporal | `downstream/*` |
| 4 | Unified QC, dashboard & export | Integrate QC metrics, render reports and dashboard | `qc/*`, `reporting/*` |

The shared infrastructure core comprises `config/` (schema, loader, parameter
registry), `release.py` (status matrix), `errors.py` (typed exception
hierarchy), and `cli.py` (command surface).

---

## 3. Module map

```
src/pfvarx/
├── __init__.py            # version, status_banner()
├── cli.py                 # 15 subcommands (argparse)
├── errors.py              # PfVarxError hierarchy
├── release.py             # RELEASE_STATUS matrix, overall_status()
├── synthetic.py           # deterministic FASTA/FASTQ/BAM/VCF generator
├── config/
│   ├── schema.py          # Pydantic v2 models, validators, config hashing
│   ├── loader.py          # YAML load, profile merge, precedence
│   └── parameters.py      # 57-entry parameter registry
├── io/
│   ├── validators.py      # file-type detection, samplesheet parsing
│   ├── routing.py         # route selection (alignment_reuse / variant_import / …)
│   ├── checksums.py       # SHA-256 helpers, provenance records
│   └── variant_import.py  # VCF/BCF/gVCF import path
├── preprocess/
│   ├── tools.py           # versioned wrappers (fastp, BWA-MEM2, samtools, Picard, fgbio)
│   └── pipeline.py        # stage composition, resume logic
├── calling/
│   ├── evidence.py        # CIGAR-aware pileup, fragment dedup, context metrics
│   ├── policy.py          # decision categories, filters, QC statuses
│   ├── caller.py          # orchestration of evidence → model → policy → output
│   └── writer.py          # VCF/BCF + 46-column evidence TSV
├── models/
│   └── count_models.py    # binomial, beta-binomial, LRT/bootstrap, Bayes, calibration
├── annotation/
│   ├── annotator.py       # internal per-call annotation; SnpEff not integrated
│   ├── consequence.py     # codon translation, variant consequence
│   ├── catalogue.py       # resistance-marker matching
│   └── gff.py             # GFF3 parsing
├── downstream/
│   ├── diversity.py       # π, Watterson's θ, Tajima's D proxies
│   ├── mixed_infection.py # intermediate-fraction binning
│   ├── relatedness.py     # IBS sharing, pairwise matrix
│   ├── phylogenetics.py   # NJ tree, Newick serialisation
│   └── spatial_temporal.py# grouping summaries
├── qc/
│   ├── metrics.py         # metric definitions, PASS/WARN/FAIL evaluation
│   └── integrator.py      # run-level QC aggregation
├── reporting/
│   ├── exporters.py       # CSV, XLSX, PNG, PDF
│   └── dashboard.py       # Streamlit file-backed dashboard
└── benchmark/
    └── runner.py          # simulation, ablation ladder, Wilson CIs
```

---

## 4. Layer contracts

### 4.1 Layer 0 — validation and routing

`io.validators.detect_type(path)` returns one of `fasta`, `fastq`,
`fastq.gz`, `bam`, `cram`, `vcf`, `gvcf`, `bcf`, or raises
`InputValidationError`. Detection combines suffix inspection with content
sniffing (BGZF magic bytes, `##fileformat` header lines) so that a mislabelled
file is rejected rather than silently misprocessed.

`io.routing.plan_run(samples, cfg)` returns a list of `RouteDecision` objects.
The route is one of:

- `alignment` — raw FASTQ requires alignment;
- `alignment_reuse` — a valid BAM/CRAM is reused as-is (never realigned);
- `variant_import` — a VCF/BCF/gVCF is imported through the variant-import
  route and never re-derived from reads;
- `rejected` — the input failed validation.

The routing contract is the mechanism by which the specification's
"reuse valid BAM/CRAM, never realign unnecessarily" and "VCF/BCF/gVCF follow
the variant-import route" requirements are honoured. A reference FASTA is
recognised as *not sample evidence* and cannot satisfy a sample slot.

### 4.2 Layer 1 — preprocessing

`preprocess.tools` exposes one adapter per external tool. Each adapter declares
the tool name, the minimum supported version, and the exact command line it
builds. `preprocess.pipeline` composes adapters into stages and records a
`StageResult` for each, including wall-clock time and exit status. When a
required binary is absent the adapter raises `ToolUnavailableError`, which the
Snakemake rules surface as a failed job rather than a silent skip.

### 4.3 Layer 2 — statistical calling

This is the scientific core and is implemented entirely in pure Python
(numpy/scipy/pysam), so it runs without any external binary.

1. **Evidence extraction** (`calling/evidence.py`) performs a CIGAR-aware
   pileup. Overlapping mate pairs are counted once; duplicate fragments,
   low-MAPQ reads, low-base-quality bases, and secondary alignments are
   excluded and counted separately. For each candidate it records reference,
   alternate, and other counts; informative depth; allele fraction with a
   Jeffreys credible interval; mean base quality; forward/reverse strand
   counts; mean read position; a Fisher exact strand-bias p-value; homopolymer
   run length; and left-alignment status.
2. **Model fitting** (`models/count_models.py`) evaluates binomial and
   beta-binomial likelihoods in log space. `evaluate_hypotheses` compares the
   null with a composite H1 (true AF at or above the configured boundary),
   maximizing over the allowed AF interval. The caller's inferential p-value is
   instead the exact one-sided binomial/beta-binomial ALT-count tail under the
   configured null, stored in legacy field `lrt_pvalue`; it applies Bonferroni
   using a complete pre-filter scan upper bound. That bound includes every
   configured possible small-indel sequence and can be very conservative. LRT
   bootstrap/mixture routines and `bh_fdr` remain model helpers, not the caller's
   inference path. Optional Bayesian output is exploratory and point-boundary
   only. Error rate and rho are treated as fixed configured null parameters;
   estimation uncertainty and context-specific errors are not represented.
3. **Decision policy** (`calling/policy.py`) maps evidence and model output to
   one of five decision categories. Hard gates and review filters are recorded
   with named reasons. Although configuration has a `VALIDATED` status, this
   caller consumes no verifiable assay-validation/LoD provenance and downgrades
   would-be PASS calls; it cannot produce a validated `ALT_PASS_VALIDATED`.
   Tests and synthetic demos establish software behaviour, not biological
   validation.
4. **Output** (`calling/writer.py`) writes a standards-compliant VCF/BCF plus a
   separate 46-column evidence TSV. Decision category and policy identifier are
   carried as `INFO` fields (`PVX_DECISION`, `PVX_POLICY_ID`); the raw LRT is
   never written to `QUAL`.

### 4.4 Layer 2.5 — annotation

Per-call `annotation.annotator` uses the internal engine only. The schema allows
`internal` or `auto` and rejects SnpEff selection; a SnpEff round-trip is not
wired into the per-call path. Consequence annotation requires the relevant GFF3
and FASTA inputs. Resistance markers are matched against
`config/catalogue/pf_resistance.yaml`; tests establish software behaviour, not
biological annotation accuracy.

### 4.5 Layer 3 — downstream

`downstream.diversity` computes nucleotide diversity, Watterson's theta, and
segregating-site counts. `downstream.mixed_infection` classifies samples by the
number of variants falling in the intermediate allele-fraction bin.
`downstream.relatedness` computes pairwise identity-by-state sharing and a
matrix. `downstream.phylogenetics` builds a neighbour-joining tree and
serialises Newick. DEploid, hmmIBD, and IQ-TREE are available through guarded
wrappers that fail closed.

### 4.6 Layer 4 — QC and reporting

`qc.metrics` evaluates each metric in `config/schemas/qc_ontology.json` against
a threshold with a declared direction (`min`, `max`, `range`, `none`) and
returns `PASS`, `WARN`, `FAIL`, or `INFO`. `qc.integrator` aggregates
per-sample metrics into a run-level summary. `reporting.exporters` produces CSV,
XLSX, PNG, and PDF; `reporting.dashboard` serves a file-backed Streamlit app.

---

## 5. Configuration and policy identity

Configuration is a Pydantic v2 model tree with `extra="forbid"`, so an unknown
key is a hard error rather than a silently ignored typo. Precedence is:

```
built-in defaults  <  mode profile (config/profiles/*.yaml)  <  run config  <  CLI overrides
```

Every loaded configuration computes a SHA-256 `config_hash()`. The parameter
registry likewise exposes `registry_hash()`. Together with the run mode these
form the policy identifier:

```
PVX-{mode}-{policy[:4]}-{config_hash[:12]}-{registry_hash[:8]}
```

The policy identifier is written into every VCF record and every QC record, so
any downstream artefact can be traced to the exact configuration and parameter
set that produced it.

---

## 6. Workflow orchestration

The Snakemake `Snakefile` imports the routing and configuration modules
directly. At parse time it validates the samplesheet, plans routes, resolves the
selected objectives, and constructs the target list from the actual samples and
objectives. Rules are included from `workflow/rules/*.smk`, one file per layer.
Because routing is computed before rule expansion, a sample that arrives as a
BAM never triggers the alignment rules, and a sample that arrives as a VCF never
triggers either alignment or calling — it goes straight to import.

Evidence production is the only route-specific part of the workflow. Once every
sample has an evidence table, a single `execute_objectives` rule hands off to
`pfvarx run --evidence-dir`, which runs exactly the objectives the user selected
and writes `run_manifest.json`. Annotation and downstream analyses are therefore
not separate workflow rules: they are objective runners inside the executor, and
they run only when their objective is selected.

Four executor profiles are provided: `local` (8 GB, 2–4 threads), `server`
(moderate parallelism), `hpc` (Slurm), and `cloud` (Kubernetes). See
`docs/DEPLOYMENT.md`.

---

## 7. Error model

All errors derive from `PfvarxError`. The hierarchy distinguishes:

- `ConfigError` — schema violation, missing required value, precedence conflict;
- `InputValidationError` — unreadable, mislabelled, or out-of-scope input;
- `RoutingError` — no valid route for a sample;
- `ModelError` — invalid count input or numerical failure;
- `ToolUnavailableError` — an external binary or database is absent.

The CLI maps these to distinct exit codes so that workflow engines can branch on
failure class.

---

## 8. Testing architecture

- **Unit tests** (`tests/unit/`) cover configuration, routing, evidence
  extraction, models, policy, writer, QC, downstream, reporting, annotation,
  variant import, and the CLI surface.
- **Integration tests** (`tests/integration/`) exercise the full
  calling → QC → report chain, atomic-write behaviour, idempotent re-runs, and
  the variant-import route.
- **Benchmark harness** (`benchmark/runner.py`) runs a simulation-based ablation
  ladder and reports precision, recall, F1, false positives in negative
  controls, and Wilson confidence intervals.

The suite is run with `pytest -q` and linted with `ruff check src tests`. See
`docs/VALIDATION.md` for the recorded results.
