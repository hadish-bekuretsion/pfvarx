# Benchmark & Ablation Results (synthetic)

**RESEARCH USE ONLY — synthetic benchmark, not clinical validation.**

Coverage: 60x · Replicates: 3 · Seed: 11

| Rung | Precision | Recall | F1 | TP | FP | FN |
|---|---|---|---|---|---|---|
| `binomial` | 0.126 | 0.889 | 0.221 | 16 | 111 | 2 |
| `beta_binomial` | 0.159 | 0.778 | 0.264 | 14 | 74 | 4 |
| `beta_binomial_lrt` | 1.000 | 0.611 | 0.759 | 11 | 0 | 7 |
| `beta_binomial_lrt_bayes` | 1.000 | 0.500 | 0.667 | 9 | 0 | 9 |
| `full_policy` | 1.000 | 0.500 | 0.667 | 9 | 0 | 9 |

## Sensitivity by truth allele fraction

| Rung | 0.02 | 0.05 | 0.10 | 0.20 | 0.30 | 0.50 |
|---|---|---|---|---|---|---|
| `binomial` | 0.33 (n=3) | 1.00 (n=3) | 1.00 (n=3) | 1.00 (n=3) | 1.00 (n=3) | 1.00 (n=3) |
| `beta_binomial` | 0.00 (n=3) | 0.67 (n=3) | 1.00 (n=3) | 1.00 (n=3) | 1.00 (n=3) | 1.00 (n=3) |
| `beta_binomial_lrt` | 0.00 (n=3) | 0.00 (n=3) | 0.67 (n=3) | 1.00 (n=3) | 1.00 (n=3) | 1.00 (n=3) |
| `beta_binomial_lrt_bayes` | 0.00 (n=3) | 0.00 (n=3) | 0.00 (n=3) | 1.00 (n=3) | 1.00 (n=3) | 1.00 (n=3) |
| `full_policy` | 0.00 (n=3) | 0.00 (n=3) | 0.00 (n=3) | 1.00 (n=3) | 1.00 (n=3) | 1.00 (n=3) |

## Interpretation

The ladder isolates each component's marginal contribution. The uncorrected
binomial rung maximises recall but admits a large number of false positives in
negative-control regions (precision 0.13). Adding beta-binomial over-dispersion
reduces these, and Benjamini–Hochberg FDR control removes the remainder entirely
(precision 1.00) at the cost of sensitivity at the lowest allele fractions. The
Bayesian and full-policy rungs add a calibrated posterior gate and review filters,
trading further sensitivity for specificity. No rung is scientifically validated;
these numbers characterise behaviour on simulated data only.

Reproduce with:

```bash
pfvarx benchmark --config config/config.yaml --replicates 3 --coverage 60 --seed 11
```
