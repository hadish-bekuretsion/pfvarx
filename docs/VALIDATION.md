# Validation and Verification

This record distinguishes software verification from scientific validation. Passing tests and synthetic workflow execution demonstrate code behavior on the tested fixtures; they do not establish variant accuracy, assay sensitivity/specificity, calibrated false-discovery rates, or clinical validity.

## 1. Current status

| Dimension | Status | Evidence / limit |
|---|---|---|
| Software correctness | TESTED_NOT_SCIENTIFICALLY_VALIDATED | 196 tests passed in the Pixi default environment; Ruff passed. |
| Pixi reproducibility | Verified on Linux ARM64 for this audit | Locked default and workflow environments installed natively. Other declared targets were not natively tested here. |
| Workflow software integration | Synthetic only | Local Snakemake dry run and paired-FASTQ end-to-end run completed; direct BAM reuse and VCF import runs completed on synthetic fixtures. |
| Statistical validity | Not scientifically validated | Exact count-tail p-values and Bonferroni bounds are implemented; error-rate/overdispersion uncertainty and empirical operating characteristics are unvalidated. |
| Annotation validity | Not scientifically validated | Transcript/phase and catalogue tests pass; no independent biological truth set or annotator comparison. |
| Clinical validity | Not established | No clinical use is supported. |

The status command's overall scientific banner remains `RESEARCH USE ONLY — validation status: NOT_IMPLEMENTED (no component is clinically validated)` because the aggregate uses the weakest scientific status.

## 2. Reproduction commands

From the repository root, with Pixi installed:

```bash
pixi install --locked -e default
pixi run -e default ruff check src tests
pixi run -e default test
pixi install --locked -e workflow
pixi run -e workflow snakemake --version
pixi run -e default status
```

Do not replace these with manually created Conda environments. The committed `pixi.lock` is authoritative for resolved package builds.

## 3. Audit evidence recorded 2026-10-23

| Check | Result | Scope |
|---|---|---|
| Native Pixi install: `pixi install --locked -e default` | Passed | Linux ARM64 Sandbox; Python 3.11. |
| Full test task: `pixi run -e default test` | **196 passed, 14 warnings** | Warnings are Matplotlib/pyparsing deprecations in figure export tests; no test failures. |
| Ruff: `pixi run -e default ruff check src tests` | Passed | All checks passed after import/loop cleanup. |
| Workflow install: `pixi install --locked -e workflow` | Passed | Linux ARM64 Sandbox. |
| Snakemake version | `9.27.1` | Workflow environment. |
| Snakemake local dry run | Passed | One synthetic paired-FASTQ sample; nine-job DAG built. |
| Synthetic paired-FASTQ pipeline | Passed | fastp, BWA-MEM2, SAMtools, Picard, PfVarX calling/objective, QC, PDF and XLSX outputs completed. Fixture used a synthetic 50 kb reference; it is not *P. falciparum* truth data. |
| Direct BAM reuse | Passed | Synthetic BAM fixture; `diversity` objective completed. |
| VCF import | Passed | Synthetic VCF fixture; `diversity` objective completed. |
| Stage 1 `pfvarx init` | Passed in CLI/unit regressions | Created configured input folders only; did not create results, logs, or work folders; rejected input/output overlap. |
| Reference FASTA checksum | Passed in CLI/unit regressions | Preflight accepted exact SHA-256 and rejected a mismatch or setup placeholder; direct-import mismatch fails before output creation. |
| Configurable results/log/work roots | Passed | Direct CLI test and a 9/9-step synthetic FASTQ workflow used separate results, logs, and work destinations. |
| Containers, HPC, real data | Not run | No Docker/Apptainer build, Slurm/Kubernetes cluster, real assay data, independent truth set, or LoD validation was available. |

In the paired-FASTQ fixture, the sample was intentionally low coverage and the conservative family adjustment caused candidates to remain `ALT_UNCERTAIN` or `NO_CALL`. That confirms fail-closed behavior for this test, not sensitivity. The configured 1–50 bp indel universe can make the Bonferroni adjustment very conservative.

## 4. What the tests establish—and do not establish

Unit/regression coverage includes configuration and Pixi-facing CLI behavior, input-only directory creation/path collision checks, configurable output roots, read evidence extraction, multiallelic VCF counts, exact REF-span checks, declared VCF contig ordering, composite-H1 likelihood fitting, exact count-tail p-values, family-bound adjustment, validated-PASS fail-closed behavior, GFF transcript-parent/phase handling, conservative generic resistance matching, and reporting outputs. The current complete suite is the source of the exact count above.

The synthetic route runs establish that selected local software components can pass data between FASTQ, BAM, VCF, calling, objectives, and report outputs. They do not prove any call is biologically true. In particular, no test here compared PfVarX to orthogonally validated *P. falciparum* genotypes or measured precision/recall on independent real data.

## 5. Statistical assumptions and scientific gaps

- The exact one-sided binomial or beta-binomial ALT-count tail is conditional on the configured null error rate and, for beta-binomial, configured `rho`. Uncertainty in estimating those parameters and context-specific errors are not propagated.
- The caller applies Bonferroni correction using an upper bound covering possible SNP and anchored-INDEL alleles in the full requested scan span, before candidate thresholds and truncation. This avoids candidate-only BH inference, but can lose substantial power; the bound does not equal an empirically measured false-discovery rate.
- The correction is per sample. No cohort-wide adjustment is made across multiple samples, so cohort-level family-wise error control is not claimed.
- `lrt_pvalue` is a legacy output field; in the current per-site caller it stores the exact count-tail p-value. Bootstrap/LRT helper functions remain available separately, but the final caller decision uses the exact count-tail path.
- Bayesian posterior output, when enabled, is exploratory and uses the configured point AF boundary; it cannot certify a validated call.
- No assay-validation or LoD provenance is represented/verified, so `ALT_PASS_VALIDATED` is not emitted.

## 6. Required work for scientific validation

Before any scientific status promotion, obtain and document: (1) orthogonally characterised positive and negative samples, (2) dilution series spanning the intended allele-fraction/depth range, (3) assay-specific error-rate and overdispersion calibration, (4) empirical calibration and null simulations through candidate discovery, (5) sensitivity, specificity, precision, recall, and false-positive estimates with confidence intervals stratified by variant type/context, (6) a validated LoD, (7) independent annotation comparison and catalogue provenance, and (8) replication on independent data. These cannot be inferred from a passing software test suite.

Until then, keep policy parameters provisional and use the tool only for research exploration with explicit human review. Clinical decisions and validated assay reporting are unsupported.
