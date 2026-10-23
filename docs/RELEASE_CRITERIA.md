# Release Criteria Checklist

This checklist separates software shipment evidence from scientific validity. The 2026-10-23 audit verified a locked Pixi install and synthetic software routes on Linux ARM64; it did not establish real-data accuracy or clinical validity.

## 1. Software release criteria

| # | Criterion | Current evidence | Status |
|---|---|---|---|
| S1 | Modular architecture | Layered modules with documented interfaces | Met for current scope |
| S2 | Configuration-driven | Pydantic schema, bounds and parameter registry | Met for current scope |
| S3 | Typed configuration | Pydantic v2; unknown fields rejected | Met |
| S4 | Provenance | Config/registry hash, policy ID and run manifest | Met in tested paths |
| S5 | Determinism | Seeded simulation and idempotence tests | Met in tested paths |
| S6 | Atomic writes/output safety | Atomic-write and file-collision regressions | Met in tested paths |
| S7 | Fail-closed wrappers | Missing/unverified tools and unsupported routes fail closed | Met in tested paths |
| S8 | Unit/regression tests | `pixi run -e default test`: 196 passed, 14 deprecation warnings | Met on Linux ARM64 |
| S9 | Integration | Synthetic FASTQ/Snakemake completed, including separate configured results/log/work roots; direct synthetic BAM and VCF routes completed; input-only `init` and direct-import checksum guard tested | Met for these fixtures only |
| S10 | Lint | `pixi run -e default ruff check src tests` passed | Met on Linux ARM64 |
| S11 | Documentation | Pixi install, usage, methods, parameters, validation and deployment docs | Updated; review changes before tagged release |
| S12 | Deployment profiles | Snakemake 9.27.1 and Slurm/Kubernetes executor distributions installed; local dry run passed | Cluster profiles not validated |
| S13 | Containers | Existing Docker and Apptainer definitions still use legacy `environment.yml`/micromamba, not `pixi.lock` | Not built/tested; not a supported Pixi deployment |
| S14 | Example data | Deterministic synthetic generator/demo | Met as software fixtures only |
| S15 | Honest status | Machine-readable status and research-use banner | Met |

The audit installed and ran the locked environments on Linux ARM64. Linux x86-64 and Intel macOS remain declared lock targets but were not natively verified here; native macOS ARM64 and Windows are not declared. See `docs/INSTALL.md`.

## 2. Scientific release criteria

No scientific criterion below is met by this software audit.

| # | Criterion | Evidence required | Status |
|---|---|---|---|
| T1 | Characterised truth set | Orthogonal genotypes/mixtures or validated multi-pipeline consensus | Not met |
| T2 | Truth provenance | Construction, source, limitations and version documented | Not met |
| T3 | Independent splits | Separate calibration and evaluation data | Not met |
| C1 | Error calibration | Assay-specific error rate fitted to controls | Not met |
| C2 | Overdispersion calibration | `beta_binomial_rho` fitted and uncertainty assessed | Not met |
| C3 | Calibration assessment | Reliability/null calibration on independent data | Not met |
| C4 | Validated LoD | Dilution series and prespecified operating point | Not met |
| P1 | SNP performance | Precision/recall/specificity with confidence intervals | Not met |
| P2 | Indel performance | Stratified by length, sequence context and repeat/homopolymer | Not met |
| P3 | Minority sensitivity | Across allele fraction, depth and assay conditions | Not met |
| P4 | False-positive rate | Negative controls and full selection/calling pipeline | Not met |
| P5 | External comparison | Same-read comparison with established callers and truth set | Not met |
| R1 | Independent replication | Independent laboratory or dataset | Not met |
| R2 | Robustness | Low depth, GC, repeats, contamination and mapping stress tests | Not met |
| R3 | Annotation validity | Independent consequence comparison and curated marker catalogue provenance | Not met |

## 3. Promotion rules

A component may be promoted only when evidence appropriate to its scope exists and `src/pfvarx/release.py` is updated. Software tests, synthetic truth, configuration labels, and manually supplied error rates do not qualify as assay validation. `ALT_PASS_VALIDATED` remains unavailable unless verifiable assay-validation/LoD provenance is represented and checked.

Clinical use is outside the platform's intended scope. No status in this checklist authorises clinical diagnosis or treatment decisions.

## 4. Reproduction commands

```bash
pixi install --locked -e default
pixi run -e default ruff check src tests
pixi run -e default test
pixi run -e default status
pixi install --locked -e workflow
pixi run -e workflow snakemake --version
```

Run real-data benchmarking only after defining and documenting its independent truth set, calibration split, and analysis plan. Never treat the bundled simulation as performance evidence.

## 5. Current verdict

| Axis | Verdict |
|---|---|
| Software | Pixi-based software and synthetic routes verified on Linux ARM64; other declared targets, cluster operation and container builds remain unverified. |
| Science | Not validated; provisional exploratory use only. |
| Overall | `NOT_IMPLEMENTED` scientific aggregate. |

PfVarX does not emit `ALT_PASS_VALIDATED` without verified validation/LoD provenance. Its p-value multiplicity adjustment is conservative and assumes configured null parameters. These gates reduce overclaiming but do not establish accuracy.
