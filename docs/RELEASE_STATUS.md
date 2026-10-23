# Release Status and Limitations

PfVarX tracks software behavior separately from scientific validity. Passing tests show that code paths behave as tested; they do not demonstrate call accuracy on biological samples. The machine-readable source is `src/pfvarx/release.py`, and `pixi run -e default status` prints the current summary.

## Status vocabulary and aggregate

| Status | Meaning |
|---|---|
| `NOT_IMPLEMENTED` | The component's scientific function is not validated or is not delivered for the stated use. |
| `IMPLEMENTED_NOT_VALIDATED` | Code exists, but no adequate scientific validation evidence is available. |
| `TESTED_NOT_SCIENTIFICALLY_VALIDATED` | Automated software tests pass; real-data scientific validation is absent. |
| `VALIDATED_WITH_RESTRICTIONS` | Validated for a defined scope with restrictions. |
| `RESEARCH_USE_READY` | Validated for research use within a stated scope. |

The aggregate status is deliberately the weakest scientific status among components. Current banner: `RESEARCH USE ONLY — validation status: NOT_IMPLEMENTED (no component is clinically validated)`. The banner is not configurable.

## Component matrix

| Component | Software status | Scientific status | Evidence / limitation |
|---|---|---|---|
| Configuration | TESTED_NOT_SCIENTIFICALLY_VALIDATED | IMPLEMENTED_NOT_VALIDATED | Pydantic validation, registry bounds, profile precedence, config hashing, input-only initialization, exact FASTA checksum verification, and configurable results/log/work roots are covered by tests. Some accepted settings remain unwired; see the audit report. |
| Input validation and routing | TESTED_NOT_SCIENTIFICALLY_VALIDATED | IMPLEMENTED_NOT_VALIDATED | FASTQ/BAM/CRAM and VCF routes are tested; BCF/gVCF support has not been independently end-to-end validated here. |
| Preprocessing and alignment | TESTED_NOT_SCIENTIFICALLY_VALIDATED | NOT_IMPLEMENTED | Pixi/Snakemake paired-FASTQ synthetic smoke run completed with fastp, BWA-MEM2, SAMtools, and Picard; no real-data accuracy or performance validation. |
| Evidence extraction | TESTED_NOT_SCIENTIFICALLY_VALIDATED | IMPLEMENTED_NOT_VALIDATED | CIGAR-aware evidence, fragment handling, and context metrics have unit tests; assay-specific calibration is absent. |
| Statistical models | TESTED_NOT_SCIENTIFICALLY_VALIDATED | IMPLEMENTED_NOT_VALIDATED | Composite alternative and exact one-sided count tails are tested; null calibration parameters are treated as fixed and are not empirically validated. |
| Decision policy | TESTED_NOT_SCIENTIFICALLY_VALIDATED | IMPLEMENTED_NOT_VALIDATED | Complete-family Bonferroni adjustment; no validated PASS without verifiable assay-validation/LoD provenance. Power and operating characteristics are unvalidated. |
| Annotation | TESTED_NOT_SCIENTIFICALLY_VALIDATED | NOT_IMPLEMENTED | Internal GFF3/CDS annotation and conservative catalogue matching have fixtures; no independent biological truth-set comparison. SnpEff is not integrated into the per-call path. |
| Downstream analyses | TESTED_NOT_SCIENTIFICALLY_VALIDATED | NOT_IMPLEMENTED | Selected summaries have unit tests; formal external tools remain guarded wrappers. |
| Unified QC and reporting | TESTED_NOT_SCIENTIFICALLY_VALIDATED | IMPLEMENTED_NOT_VALIDATED | Synthetic FASTQ workflow produced manifest, QC, PDF and XLSX outputs; no scientific validation of metrics. |
| Benchmarking | TESTED_NOT_SCIENTIFICALLY_VALIDATED | NOT_IMPLEMENTED | Simulation framework is software, not external-truth evidence. |

## Verified software scope

The current audit run used Pixi on native Linux ARM64. The locked `default` and `workflow` environments installed successfully; the default test task passed 196 tests (14 upstream Matplotlib/pyparsing deprecation warnings), Ruff passed, and Snakemake 9.27.1 was available. A local Snakemake dry run and synthetic paired-FASTQ runs completed, including one with separate configured results/log/work roots. Direct BAM and VCF-import runs also completed against synthetic fixtures. Stage 1 `pfvarx init` and exact FASTA checksum verification, including direct-import mismatch rejection, are covered by CLI/unit tests. See `docs/AUDIT_REPORT_2026-10-23.md` for command-level evidence and remaining gaps.

The workflow includes Slurm and Kubernetes executor plugin distributions, but no cluster deployment or cluster profile has been validated. Linux x86-64 and Intel macOS are declared lock targets but were not natively installed in this audit; native macOS ARM64 and Windows are not declared targets. Docker and Apptainer images were not built.

## Current scientific limitations

1. No empirical error-rate/overdispersion calibration, orthogonal truth set, validated limit of detection, or independent replication is supplied. The caller's p-values assume the configured null parameters are correct and fixed.
2. The caller tests a composite alternative (true AF at least the configured boundary) and reports an exact one-sided binomial/beta-binomial ALT-count tail. A complete pre-filter hypothesis bound is Bonferroni-adjusted per sample; the correction can be very conservative, particularly with the default 1–50 bp indel universe, and there is no cohort-wide adjustment. The legacy field name `lrt_pvalue` is retained for compatibility; in the caller it contains the exact count-tail p-value, not an LRT-calibrated p-value. Bayesian output is exploratory and uses a point AF boundary.
3. `ALT_PASS_VALIDATED` is unreachable without a validation/LoD provenance object, which is not currently represented or verified. A `VALIDATED` config label or manually supplied error rate is not evidence.
4. Broad PfKelch13/PfATP4 nonsynonymous catalogue descriptions can produce only `low_confidence_candidate` review entries. They are not variant-specific resistance determinations. The catalogue requires expert curation against a current, authoritative list.
5. GFF3 handling resolves transcript parents and CDS phase, but selects the lexicographically smallest transcript ID per gene. That deterministic engineering rule is not a biological canonical-isoform assertion. No independent annotator comparison was run.
6. UMI consensus calling is not established; the policy fails closed when consensus validity is unavailable. SnpEff is not an accepted per-call engine; configuration allows only `auto` and `internal`.
7. Structural variants, copy number, long-read/single-end processing, and clinical interpretation are outside the validated scope. The reference genome is not bundled.

## Intended use

Use only for exploratory research, method development, teaching, and software evaluation, with an independently reviewed reference and explicit review of evidence tables and uncertainty. Do not use for diagnosis, treatment selection, patient-facing decisions, or validated assay reporting.

Before scientific status can be promoted, the project needs assay-specific calibration, truth/dilution series, empirical error and overdispersion assessment, sensitivity/specificity and false-positive studies (including selection and indel-context effects), a validated LoD, external-annotator comparisons, independent replication, and a documented scope of validity. See `docs/VALIDATION.md` and `docs/RELEASE_CRITERIA.md`.
