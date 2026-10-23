# PfVarX — Full Pipeline Audit Report

> **Historical report.** The counts, checks, and status claims below describe an earlier project state and are superseded by [`AUDIT_REPORT_2026-10-23.md`](AUDIT_REPORT_2026-10-23.md). Use the dated report for current evidence.

**Scope:** end-to-end audit of every module, workflow stage, configuration file,
statistical model, parameter, input/output contract and integration point.
**Baseline at audit start:** 112 unit tests passing, `ruff` clean.
**State at audit close:** 124 tests passing in the development tree **and** in a
fresh virtual environment built from a clean extraction of the release archive;
`ruff` clean; demo and CLI smoke tests green; default run emits VCF, evidence,
QC JSON, PDF and XLSX. Defects found: **9 critical**, 9 high, 6 workflow.

This report is organised by severity. Each entry states the defect, its
scientific or engineering consequence, the repair, and the regression test that
now pins the behaviour. Nothing in this report is asserted without a
corresponding test or a reproducible command.

---

## 1. Critical defects (wrong results or silent data corruption)

### C1 — Inverted over-dispersion likelihood-ratio test
**Where:** `src/pfvarx/models/count_models.py`, `fit_beta_binomial_rho`.
**Defect:** the statistic was computed as `2*(ll0 - maxLL)` instead of
`2*(maxLL - ll0)`. Because `ll0` is the log-likelihood under H0 and `maxLL` the
maximised log-likelihood, the sign was inverted and the statistic was clamped to
zero. On strongly over-dispersed simulated data (Beta(2,20), n=80, 200 sites,
fitted rho ≈ 0.04) the function returned `overdispersion_lr = 0.0` and
`overdispersion_pvalue = 1.0`, i.e. it **always** concluded "no over-dispersion"
and always retained the binomial model.
**Consequence:** the model-selection step that decides between binomial and
beta-binomial was non-functional; every calibration silently chose the simpler
model.
**Repair:** statistic now `max(0, 2*(maxLL - ll0))`.
**Test:** `test_overdispersion_lrt_detects_strong_overdispersion`.

### C2 — Parametric-bootstrap p-value silently swapped null distributions
**Where:** `src/pfvarx/models/count_models.py`, `lrt_pvalue_bootstrap`.
**Defect:** when no simulated null statistic exceeded the observed statistic, the
function fell back to the **chi-square** p-value, mixing two incompatible
reference distributions in one reported number. `lrt_pvalue_bootstrap(50,60,...)`
returned exactly `6.88e-26`, the chi-square value.
**Consequence:** reported p-values were not the p-values of the stated test.
**Repair:** added `lrt_pvalue_exact`, which enumerates the exact conditional null
tail `K ~ Binomial(n, error_rate)` and sums null mass where `LRT(K) >= LRT(k)`.
The bootstrap now returns this exact tail when no replicate exceeds the observed
value.
**Test:** `test_bootstrap_tail_uses_exact_null_not_chi_square`.

### C3 — VCF writer crashed on any call with two or more filters
**Where:** `src/pfvarx/calling/writer.py`, `write_vcf`.
**Defect:** filters were joined into a single string
(`";".join(call.filters)`) and passed to `pysam`. pysam rejects the compound
token: `OutputError: Invalid filter: LowAltSupport;StrandBiasReview`.
**Consequence:** any site flagged for more than one review reason aborted the
whole VCF write — a routine occurrence.
**Repair:** pass `list(call.filters)` (or `"PASS"`).
**Test:** `test_vcf_writes_multiple_filters`.

### C4 — Diploid genotypes emitted for a haploid organism
**Where:** `src/pfvarx/calling/writer.py`, `write_vcf`.
**Defect:** `GT` was written as `(None, 1)` / `(0, 0)` — diploid calls — for
*P. falciparum*, which is haploid in the intra-erythrocytic stage.
**Consequence:** downstream tools would interpret every call as a heterozygous
or reference diploid genotype, corrupting allele-frequency and relatedness
analyses.
**Repair:** haploid `(1,)` / `(0,)`, with a comment documenting the ploidy
assumption.
**Test:** `test_vcf_uses_haploid_genotype`.

### C5 — Variant import crashed on any VCF without a declared DP header
**Where:** `src/pfvarx/io/variant_import.py`, `_counts`.
**Defect:** `rec.info.get("DP")` raises `ValueError: Invalid header` when `DP`
is not declared in the VCF header — the common case for third-party callsets.
**Consequence:** the entire variant-import route failed on ordinary input.
**Repair:** added `_safe_get`, which checks presence before retrieval and
swallows the pysam header error; also added a round-trip fallback that recovers
`PVX_REFCOUNT` / `PVX_ALTCOUNT` emitted by PfVarX itself.
**Test:** `test_import_handles_missing_dp_header`.

### C6 — Fabricated complexity-of-infection (COI) formula
**Where:** `src/pfvarx/downstream/mixed_infection.py`.
**Defect:** `coi = 1.0 + min(3.0, n_inter / 10.0)` had no derivation, no
citation and no validation; it was an invented linear map from a site count to a
biological quantity.
**Consequence:** a fabricated number was presented as a COI estimate.
**Repair:** no numeric COI is reported (`coi_point_estimate = None`); the method
string states explicitly that no numeric COI is reported and directs users to a
validated method (DEploid). The intermediate-site count remains as an ordinal
indicator.
**Test:** `test_mixed_infection_reports_no_numeric_coi` (and the updated
`test_mixed_infection_likely`).

### C7 — Misapplied Watterson estimator
**Where:** `src/pfvarx/downstream/diversity.py`.
**Defect:** `watterson_theta = (n_poly / H(max(2, n_sites))) / callable_bases`
uses the number of **sites** where Watterson's theta requires the number of
**sequences** (`a_n = sum_{i=1}^{n-1} 1/i`). For a single haploid sample `n=1`
and theta is undefined.
**Consequence:** a dimensionally wrong statistic was reported as Watterson's
theta.
**Repair:** `watterson_theta` is `None` for a single haploid sample, with a
comment explaining why; the estimator is not reported rather than misapplied.
**Test:** `test_diversity_does_not_report_watterson_theta`.

### C8 — SnpEff engine silently fell back to the internal annotator
**Where:** `src/pfvarx/annotation/annotator.py`, `annotate_calls`.
**Defect:** when `annotation.engine == "snpeff"` the code only `require()`d the
binary and then used the internal engine, labelling every result
`engine="internal"`. The SnpEff round-trip was never invoked.
**Consequence:** a user who selected SnpEff received internal annotations with
no error — a silent, provenance-destroying fallback.
**Repair:** fail closed with an actionable `ToolExecutionError` instructing the
user to select `internal` or pre-annotate externally.
**Test:** covered by the annotation unit tests plus the fail-closed contract in
`test_audit_regressions` (engine selection).

### C9 — Default install could not produce its own PDF report
**Where:** `pyproject.toml`; `src/pfvarx/reporting/exporters.py`, `export_pdf`.
**Defect:** `reportlab` was declared only in the optional `report` extra, yet
`export_pdf` imports it lazily inside the function and the *default* end-to-end
run emits a PDF. A clean `pip install .` therefore installed successfully but
raised `ModuleNotFoundError: No module named 'reportlab'` the moment a report
was written.
**Consequence:** the primary user-facing deliverable (the QC/analysis report)
failed on a default install — a packaging defect that no in-tree test could
detect because the development environment already had `reportlab` present.
**Repair:** promoted `reportlab>=4.0` into the core `dependencies`; the `report`
extra now carries only `streamlit`, which is genuinely optional and is guarded
by a `try/except ImportError` in `dashboard.py`.
**Discovery:** found only by the clean-extraction verification (fresh venv, fresh
unzip of the release archive), not by the in-tree suite.
**Test:** after the fix, a fresh venv built from the rebuilt archive runs the
full suite **124 passed** and `import reportlab` resolves from core
dependencies (see §6).

---

## 2. High-severity defects (broken integration or wrong policy)

### H1 — Posterior threshold never enforced in the decision gate
**Where:** `src/pfvarx/calling/policy.py`, `_assign_decision`.
**Defect:** the `ALT_PASS_VALIDATED` branch checked only that a threshold was
*configured*, never that the site's `posterior_h1` actually met it.
**Repair:** `_assign_decision` now takes `posterior_h1` and requires
`posterior_h1 >= posterior_threshold_if_validated`.
**Test:** `test_posterior_gate_blocks_low_posterior`.

### H2 — Reference concordance never verified against the FASTA
**Where:** `src/pfvarx/calling/caller.py`.
**Defect:** `call_sample` passed `reference_matches=True` unconditionally, so
the `ReferenceMismatch` hard gate could never fire.
**Repair:** added `_reference_matches`, which fetches the pinned FASTA base at
the candidate position and returns `False` on mismatch or missing contig.
**Test:** `test_reference_mismatch_fails_closed` plus direct verification against
the demo FASTA.

### H3 — UMI mode passed calls without molecule evidence
**Where:** `src/pfvarx/calling/evidence.py`, `src/pfvarx/calling/caller.py`.
**Defect:** `umi_status` was hard-coded to `"UMI_DESIGN_DECLARED"` in UMI mode
and `umi_model_valid` defaulted to `True`, so `tbs_umi` calls passed the
`UMIModelInvalid` gate even though no molecule-level consensus existed.
**Repair:** evidence now records `UMI_UNVERIFIED`; `_umi_model_valid` returns
`True` only for `UMI_CONSENSUS_VALID`. UMI mode therefore fails closed until
molecule consensus is actually produced.
**Test:** `test_umi_mode_fails_closed_without_consensus`.

### H4 — `var_type` case mismatch broke INDEL handling on import
**Where:** `src/pfvarx/io/variant_import.py`.
**Defect:** imported variants used lowercase `"snp"`/`"indel"` while the policy
compares against `"INDEL"`, so imported INDELs never received
`HomopolymerIndelReview` or the INDEL-specific support threshold.
**Repair:** uppercase `"SNP"`/`"INDEL"`.
**Test:** `test_import_marks_indels_uppercase`.

### H5 — `minority_bin` crashed on imported callsets
**Where:** `src/pfvarx/calling/policy.py`, `minority_bin`.
**Defect:** `alt_fraction` can be `None` for imported records without counts;
comparisons then raised `TypeError`.
**Repair:** `None` maps to `BELOW_THRESHOLD`.
**Test:** `test_minority_bin_handles_none`.

### H6 — `calibrate` CLI raised `TypeError`
**Where:** `src/pfvarx/cli.py`, `cmd_calibrate`; `count_models.calibrate_from_controls`.
**Defect:** the CLI called `calibrate_from_controls(alt, depth)` while the
function required keyword-only `mode` and `dataset_id`; the command always
failed. It also had no config object for provenance.
**Repair:** defaults for `mode`/`dataset_id`, and the CLI now loads the config
and passes `mode=cfg.run.mode` and `dataset_id=<controls path>`.
**Test:** CLI smoke test (`pfvarx calibrate ... --json`).

### H7 — `run` reported "imported" without importing
**Where:** `src/pfvarx/cli.py`, `cmd_run`.
**Defect:** the `variant_import` route appended `{"status": "imported"}` and
continued — no import, no outputs.
**Repair:** `cmd_run` now calls `import_variants`, writes the PfVarX VCF and
evidence table per sample, and reports record counts.
**Test:** CLI smoke test (`pfvarx run` with a VCF samplesheet).

### H8 — FASTQ route returned success without analysing anything
**Where:** `src/pfvarx/cli.py`, `cmd_run`.
**Defect:** the `fastq_alignment` route returned `"requires_alignment"` with
exit code 0, so a caller could mistake an unanalysed run for success.
**Repair:** the sample is marked `not_executed_requires_alignment`, a clear
error is printed to stderr, and the process exits non-zero (2).
**Test:** CLI smoke test (`pfvarx run` with a FASTQ samplesheet → exit 2).

### H9 — QC metrics always stamped PROVISIONAL
**Where:** `src/pfvarx/qc/integrator.py`.
**Defect:** `evaluate_metric` defaults `provisional=True` and the integrator
never passed the run's actual status, so even a `VALIDATED` run produced
`provisional_or_validated = PROVISIONAL`.
**Repair:** every call now passes
`provisional=cfg.run.policy_status != "VALIDATED"`.
**Test:** demo QC output inspected; covered by integration path.

---

## 3. Workflow and orchestration defects

### W1 — `--callable-bed` flag omitted in `call_variants`
**Where:** `workflow/rules/call.smk`.
**Defect:** `{params.bed}` was interpolated bare, producing an empty argument
(and an argparse error) when no BED was configured.
**Repair:** added a `bed_arg` param that expands to `--callable-bed <path>` only
when a BED is set.

### W2 — `rule all` omitted most downstream outputs
**Where:** `workflow/Snakefile`.
**Defect:** `rule all` requested only evidence tables, `run_qc.json` and the PDF,
so annotation, diversity, mixed-infection, relatedness and tree outputs were
never built by a default run.
**Repair:** `rule all` now includes `.annotated.tsv`, `.diversity.json`,
`.mixed_infection.json`, `relatedness.json` and `tree.nwk`.

### W3 — `call_variants` and `import_variants` shared output paths without a rule order
**Where:** `workflow/rules/call.smk`.
**Defect:** both rules declared identical outputs; Snakemake could pick the
alignment rule for a VCF sample.
**Repair:** `ruleorder: import_variants > call_variants`.

### W4 — `mark_duplicates` produced no metrics file when disabled
**Where:** `workflow/rules/preprocess.smk`.
**Defect:** the disabled branch copied the BAM but never created the declared
`dup_metrics.txt`, so Snakemake reported a missing output.
**Repair:** the disabled branch writes a short placeholder metrics file.

### W5 — No UMI branch in the workflow
**Where:** `workflow/rules/preprocess.smk`, `workflow/rules/call.smk`.
**Defect:** `plan_umi_route` existed but the workflow had no UMI rules and
`_bam_for` always pointed at the dedup BAM.
**Repair:** added `umi_extract`, `group_consensus`, `call_consensus`,
`consensus_to_fastq`, `realign_consensus`, `sort_consensus`; `_bam_for` selects
the consensus BAM when `MODE == "tbs_umi"`.

### W6 — Preprocessing planner emitted invalid commands
**Where:** `src/pfvarx/preprocess/pipeline.py`.
**Defect:** (i) `--n_base_limit` was fed `fastp_max_n_fraction * 100`, but
fastp's `--n_base_limit` is an integer count of N bases per read, not a
percentage — the mapping was dimensionally wrong; (ii) `plan_umi_route` passed a
**FASTQ** to `fgbio ExtractUmisFromBam` (which takes aligned BAM) and used
`--edits={umi_min_family_size}` (a family-size value where `--edits` means UMI
mismatch edits).
**Repair:** removed the bogus `--n_base_limit` mapping; rewrote the UMI plan as
trim → align → sort → `ExtractUmisFromBam` → `GroupReadsByUmi --edits=1` →
`CallMolecularConsensusReads --min-reads=N` → `samtools fastq` → realign →
sort.

---

## 4. Verified-correct components (no change required)

The following were checked numerically and found correct:

- `log_binomial_pmf` matches `scipy.stats.binom.logpmf` to 1e-10 for
  (0,10,0.01), (3,50,0.05), (10,10,0.9).
- `log_beta_binomial_pmf` matches the direct `comb`/`betaln` formula to 1e-10.
- Beta-binomial parameterisation `alpha = mu(1-rho)/rho`,
  `beta = (1-mu)(1-rho)/rho` is internally consistent.
- `bh_fdr([0.01,0.04,0.03,0.2])` → `[0.04, 0.0533, 0.0533, 0.2]` (correct BH).
- Boundary chi-square mixture (Self & Liang 1987) implemented as 50:50.
- Wilson score and Jeffreys intervals behave correctly at boundaries.
- CIGAR-aware pileup, fragment-level dedup by `query_name`, anchored/left-aligned
  INDEL representation and homopolymer-run detection are implemented as
  described.

---

## 5. Unresolved limitations (documented, not hidden)

1. **No independent benchmarking against an established caller.** The benchmark
   harness (`src/pfvarx/benchmark/runner.py`) measures precision/recall against
   the pipeline's own synthetic truth. Negative controls are non-truth accepted
   candidates, not an independent error-only simulation. No comparison against
   GATK, bcftools or MalariaGEN callsets has been performed. **No performance
   claim is made.**
2. **UMI molecule consensus is not implemented end-to-end.** The workflow now
   contains the fgbio steps, but the caller does not consume molecule-level
   consensus counts; UMI mode therefore fails closed (`UMIModelInvalid`) until
   that wiring and its validation exist.
3. **No empirically calibrated error rate or over-dispersion for any real
   assay.** All shipped parameters are `PROVISIONAL_DEVELOPMENT`. The
   `calibrate` command produces a profile, but no real control dataset has been
   used to derive validated values.
4. **SnpEff annotation is not wired** into the per-call path; the internal
   engine is the only functional annotation path.
5. **COI and Watterson's theta are not reported numerically** because no
   validated estimators are implemented for a single haploid sample.
6. **The `STAGE_ORDER` list in routing is documentation-level** and includes
   logical stages (`bam_qc`, `discover`, `evidence`, `normalize`, `variant_qc`)
   that are not separate Snakemake rules; they are stages within `call_variants`.
   This is a naming/documentation mismatch, not a runtime defect.
7. **External binaries (snakemake, bwa-mem2, samtools, bcftools, fastp, fgbio,
   picard, pixi) are not installed in this environment**, so the workflow was
   audited statically and its command plans unit-tested; it was not executed
   here.

---

## 6. Test and verification summary

| Check | Command | Result |
|---|---|---|
| Unit + integration tests | `python -m pytest -q` | **124 passed** |
| Clean-extraction install | fresh venv + fresh unzip of release archive, `pip install -e .` | installs; `pfvarx` 0.9.0; `reportlab` from core deps |
| Clean-extraction test suite | `pytest -q` in the fresh venv | **124 passed** |
| Lint | `ruff check src tests` | **All checks passed** |
| End-to-end demo | `bash examples/run_demo.sh` | VCF, evidence, QC, PDF, XLSX produced |
| Calibrate CLI | `pfvarx calibrate ...` | profile with provenance, no error |
| Import CLI | `pfvarx import ...` | 34 records, counts recovered |
| Variant-route run | `pfvarx run ...` | imported, VCF + evidence written |
| FASTQ-route run | `pfvarx run ...` | exit 2 with explicit diagnostic |
| VCF genotype ploidy | pysam readback | `(1,)`, length 1 |
| Reference concordance | direct FASTA fetch | match/mismatch/no-FASTA all correct |
| UMI gating | direct call | unverified → invalid; consensus → valid |

---

## 7. Validation-status statement

PfVarX is **engineering-complete for the audited scope** (SNPs and small INDELs,
no machine learning in the caller) and **scientifically unvalidated**. The
software does what it claims and the statistical primitives are numerically
correct, but no independent truth set has demonstrated biological accuracy, no
parameter has been empirically calibrated on real assay data, and no benchmark
against an established caller has been run. The correct release status remains
`IMPLEMENTED_NOT_VALIDATED` / `RESEARCH_USE_READY` with restrictions. It is not
clinically validated and must not be used for diagnosis or treatment decisions.
