# PfVarX end-to-end audit — 2026-10-23 (final verification)

## Executive verdict

The uploaded PfVarX 0.9.0 project has been corrected and re-tested as a **Pixi-managed research software pipeline**. On the current native Linux ARM64 Sandbox, the locked `default` and `workflow` environments installed successfully, all 196 automated tests passed, Ruff passed, and the local Snakemake workflow completed synthetic paired-FASTQ runs. The follow-up run also verified independent configured results, logs, and work directories, including fastp reports under the work root. Direct synthetic BAM and VCF-import runs completed in the earlier audit; Stage 1 initialization and exact FASTA checksum verification are covered by unit/CLI regressions.

This is **not scientific or clinical validation**. No real *P. falciparum* reads, external truth set, assay calibration data, validated limit of detection (LoD), or independent annotation truth set was supplied. The project is usable for exploratory analysis and software evaluation with user-supplied references and careful review; it is not supported for diagnosis, treatment decisions, validated assay reporting, or an unqualified production-science claim.

## Scope and evidence

The audit began from the supplied `pfvarx-0.9.0.zip`, extracted at `/home/ubuntu/work/pfvarx`. It reviewed Python source, tests, Pixi/Python manifests, lockfile, config/catalogue resources, Snakemake rules/profiles, scripts, and user documentation. All end-to-end fixtures created for verification were synthetic and stored under `/tmp`; none represents *P. falciparum* biological truth.

### Verification performed

| Check | Result | Exact scope / limitation |
|---|---|---|
| `pixi install --locked -e default` | **Passed** | Native Linux ARM64 Sandbox; Pixi Python 3.11 environment. |
| `pixi run -e default test` | **196 passed, 14 warnings** | 14 Matplotlib/pyparsing deprecation warnings from report-figure tests; zero failures. |
| `pixi run -e default ruff check src tests` | **Passed** | All Ruff checks passed after cleanup. |
| Stage 1 `pfvarx init` | **Passed in CLI/unit tests** | Configured input directories are created; results, logs, and work locations are left absent; overlap is rejected before creation. |
| Reference FASTA checksum preflight | **Passed in CLI/unit tests** | Exact file SHA-256 accepted; mismatches and setup placeholders rejected before analysis, including direct import before output creation. |
| `pixi install --locked -e workflow` | **Passed** | Native Linux ARM64 Sandbox; includes Snakemake, Bioconda tools, and executor plugins. |
| `pixi run -e workflow snakemake --version` | **Passed: 9.27.1** | Workflow environment. |
| Local Snakemake dry run | **Passed** | Nine-job DAG for one synthetic paired-FASTQ sample and the local profile. |
| Synthetic paired-FASTQ workflow | **Passed: 9/9 steps** | fastp → BWA-MEM2 → SAMtools/Picard → PfVarX call/objective → QC → PDF/XLSX report; follow-up run used a synthetic 50 kb reference and separate configured results/log/work roots, with fastp HTML/JSON outputs inside the work root. |
| Direct BAM route | **Passed** | Synthetic BAM input; `diversity` objective completed. |
| VCF-import route | **Passed** | Synthetic VCF input; `diversity` objective completed. |
| Docker/Apptainer, Slurm/Kubernetes cluster | **Not run** | Container build was not run; executor plugins are installed but no cluster/profile execution was possible. |
| Real-data/statistical validation | **Not run** | No real reads, truth set, controls, dilution series, assay calibration, validated LoD, or independent replicates. |

The lock declares `linux-64`, `linux-aarch64`, and `osx-64`. The actual native install/test above is for `linux-aarch64`; Linux x86-64 and Intel macOS were not natively tested in this audit. Native macOS ARM64 and Windows are not declared targets.

## Verified fixes and implementation changes

1. **Pixi reproducibility and deployment.** The current manifest uses `[workspace]`, keeps Pixi as the install/task mechanism, commits the lock, adds `linux-aarch64`, and resolves both default and workflow environments. Snakemake is constrained to 9.x; compatible Slurm and Kubernetes executor plugin distributions are locked. Plugins being installed only makes the executor interface available; it is not cluster validation.
2. **Workflow orchestration.** The Snakefile now takes the effective config path from Snakemake's `workflow.configfiles`, rather than an undefined Python variable. The local dry run and synthetic FASTQ run verified the local DAG and data flow. The paired-read group in BWA-MEM2 uses the actual sample wildcard.
3. **Statistical alternative.** The binomial alternative is fitted over true allele fraction (AF) at or above the configured boundary; beta-binomial optimization also maximizes over that range. High-AF variants no longer use only the fixed minimum-AF likelihood.
4. **Selection/multiplicity safety.** The caller reports exact one-sided binomial or beta-binomial ALT-count tail p-values under the configured null. It applies Bonferroni correction using an upper bound over the complete possible SNP and anchored-INDEL hypothesis universe in the requested scan span, before candidate filtering/truncation. Candidate-only BH inference is not used in the call decision. The output retains the legacy column name `lrt_pvalue`; despite that name, the current caller stores the exact count-tail p-value there.
5. **Validated-call gate.** The policy downgrades any legacy validated-PASS object and cannot emit `ALT_PASS_VALIDATED` without verifiable assay-validation and LoD provenance. Such provenance is not currently represented or checked; a YAML status label or manually entered error rate is not sufficient.
6. **Annotation correctness and conservative matching.** GFF3 CDS `Parent` relationships are resolved through transcripts to genes; phase consistency and initial partial-codon phase are handled, and minus-strand SNP alleles are complemented. Multiple transcripts are reduced deterministically to the lexicographically smallest transcript ID (an engineering tie-break, not a biological canonical-transcript assertion). Broad PfKelch13/PfATP4 nonsynonymous placeholders can only create `low_confidence_candidate` review entries, never inherit `established` evidence.
7. **Fail-fast unsupported annotation engine.** Config accepts only `annotation.engine: auto` or `internal`. SnpEff is present in the workflow environment but is not integrated into PfVarX per-call annotation.
8. **Earlier audit fixes retained.** These include multiallelic VCF allele-specific count handling, complete REF-span checks, VCF contig-dictionary-aware sort validation, indel-specific evidence thresholds/review gates, early objective/catalogue/GFF preflight, output-file overwrite protection, correct objective/run exit behavior, and Pixi-based user scripts/docs.
9. **Two-stage input/output layout.** Added `pfvarx init` to create input folders and declared resource parents only, with an overlap guard against all configured output locations. Added `outputs.logs_dir` and `outputs.work_dir` handling to direct CLI planning/run logs and to the Snakemake FASTQ/validation/call/report rules. The follow-up synthetic workflow verified that the manifest/results, logs, and alignment intermediates landed in three distinct configured roots.
10. **Reference identity enforcement.** Validation, full runs, direct calling, direct import, and direct annotation now recompute the configured FASTA's SHA-256 and fail before reference-dependent analysis on a setup placeholder, missing FASTA, or digest mismatch. The `init` command remains usable with the example placeholder so Stage 1 can create folders before the reference is placed.

## Remaining statistical and scientific risks

These are the unresolved requirements for scientific validation, not uncorrected claims of tested functionality.

- **Null calibration is conditional on fixed parameters.** Exact count-tail validity assumes the configured sequencing error rate and beta-binomial `rho` are correct. Uncertainty in their estimation and context-specific error (sequence context, mapping, homopolymers, read position, library/run variation) is not propagated. Empirical control data are required.
- **The family bound is deliberately conservative.** It includes all possible insertion/deletion sequences up to `candidate.max_small_indel_bp` (default 50 bp) across the requested span, so it can severely reduce power. In the low-coverage synthetic FASTQ smoke run, candidates remained `ALT_UNCERTAIN`/`NO_CALL`; adjusted p-values reached 1. This demonstrates conservative software behavior for that fixture, not sensitivity or an empirical false-discovery guarantee.
- **Multiplicity control is per sample.** The Bonferroni family is the requested scan for one sample; the pipeline does not combine all sample-level hypothesis families into a cohort-wide family-wise error guarantee.
- **No empirical FDR/sensitivity/specificity or LoD is established.** The correction is a Bonferroni upper bound under the configured count model; it is not a measured assay-level FDR. The Bayesian posterior remains exploratory and uses a point AF boundary, not a prior over the composite alternative.
- **Resistance interpretation remains research-only.** Broad K13/PfATP4 rules are not residue/domain-specific and produce low-confidence review candidates for nonsynonymous substitutions. Catalogue citations and marker scope require expert review against a current authoritative catalogue; no external truth comparison was performed.
- **Transcript choice is deterministic, not biological canonicalization.** Only one transcript consequence per gene is emitted. The parser has synthetic plus/minus/phase fixtures but no cross-check against an independent annotator or a curated Pf3D7 annotation release.
- **UMI mode is not established.** fgbio is installed for its preprocessing branch, but the caller fails closed unless valid consensus status is supplied; UMI consensus calling and assay performance have not been validated.
- **Deployment boundaries remain.** Slurm/Kubernetes plugins are locked but cluster profiles are untested; Docker/Apptainer builds are untested and their legacy definitions still install from `environment.yml` with micromamba rather than `pixi.lock`, so they are explicitly not part of the supported Pixi path. No performance/resource study on full parasite genomes was performed. CRAM, BCF, and gVCF were not independently end-to-end smoke-tested in this audit.
- **No reference is bundled.** Users must supply a pinned assembly, matching FASTA/checksum, indexes/resources, and annotation files when selected objectives require them.
- **Some schema-accepted settings are not runtime controls.** The revised README identifies the observed dormant/unwired fields, including `inputs.require_checksums`, `preprocessing.fastp_max_n_fraction` and `primer_bed`, `statistics.empirical_error_profile` and `bootstrap_replicates`, several extended callability/reporting settings, benchmark target values, `resources.concurrent_samples`, and output retention/atomic-write flags. Users must not assume those values alter execution until their code paths are implemented and tested. This remains a software/configuration completeness issue.

## Recommended use and release status

The project is **installable and software-smoke-tested for exploratory research on the verified Linux ARM64 environment**. Stage 1/Stage 2 directory behavior and independent result/log/work paths now have regression and synthetic workflow coverage. The Pixi lock also declares Linux x86-64 and Intel macOS, but those were not natively executed here. Users can analyze their own data only after reviewing the input/reference match, evidence table, filters, model assumptions, and uncertainty; all profiles remain provisional.

The overall scientific banner remains `NOT_IMPLEMENTED`, which is the intentionally conservative weakest-component status. Do not use output for clinical, treatment, or validated assay decisions. To promote scientific status, obtain independent truth data and controls, calibrate error/rho, validate LoD and performance stratified by variant class/context, evaluate the complete discovery/selection process empirically, cross-check annotations/catalogue provenance, and replicate independently.

See `docs/INSTALL.md` for Pixi setup and input preparation, `docs/VALIDATION.md` for command-level evidence, and `docs/RELEASE_CRITERIA.md` for promotion gates.
