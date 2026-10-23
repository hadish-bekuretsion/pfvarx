"""Synthetic benchmark and ablation framework.

The runner builds synthetic *P. falciparum*-like haplotypes with known truth,
simulates Illumina paired-end reads at a requested depth, calls variants with
each successive model configuration, and reports precision/recall/F1 with
Wilson confidence intervals.

Ablation ladder (each rung adds one component):

    binomial  ->  beta_binomial  ->  + LRT  ->  + Bayesian  ->  full policy

Results are written to ``out_dir`` as JSON plus a CSV summary so that the
numbers quoted in the documentation are always traceable to a reproducible run.
"""

from __future__ import annotations

import json
import random
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

from pfvarx.calling.caller import call_sample
from pfvarx.config.schema import PfVarXConfig
from pfvarx.models.count_models import wilson_interval
from pfvarx.synthetic import Variant, simulate_dataset

ABLATION_LADDER = [
    "binomial",
    "beta_binomial",
    "beta_binomial_lrt",
    "beta_binomial_lrt_bayes",
    "full_policy",
]


@dataclass
class SiteOutcome:
    chrom: str
    pos: int
    ref: str
    alt: str
    truth_fraction: float
    called: bool
    decision: str
    alt_fraction: float | None
    lrt_pvalue_adjusted: float | None
    is_truth: bool = True


@dataclass
class BenchmarkResult:
    config_label: str
    replicates: int
    coverage: int
    true_positives: int = 0
    false_positives: int = 0
    false_negatives: int = 0
    truth_sites: int = 0
    called_sites: int = 0
    per_fraction: dict[str, dict[str, float]] = field(default_factory=dict)
    outcomes: list[dict[str, Any]] = field(default_factory=list)

    def metrics(self) -> dict[str, Any]:
        tp, fp, fn = self.true_positives, self.false_positives, self.false_negatives
        precision = tp / (tp + fp) if (tp + fp) else 0.0
        recall = tp / (tp + fn) if (tp + fn) else 0.0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0
        p_lo, p_hi = wilson_interval(tp, tp + fp) if (tp + fp) else (0.0, 0.0)
        r_lo, r_hi = wilson_interval(tp, tp + fn) if (tp + fn) else (0.0, 0.0)
        return {
            "precision": precision,
            "precision_ci": [p_lo, p_hi],
            "recall": recall,
            "recall_ci": [r_lo, r_hi],
            "f1": f1,
            "true_positives": tp,
            "false_positives": fp,
            "false_negatives": fn,
            "truth_sites": self.truth_sites,
            "called_sites": self.called_sites,
        }

    def as_dict(self) -> dict[str, Any]:
        return {
            "config_label": self.config_label,
            "replicates": self.replicates,
            "coverage": self.coverage,
            "metrics": self.metrics(),
            "per_fraction": self.per_fraction,
            "outcomes": self.outcomes,
        }


def _design_variants(seq: str, rng: random.Random, n_sites: int) -> list[Variant]:
    """Place truth variants at a spread of allele fractions and positions.

    REF alleles are taken from the actual reference sequence so that the
    simulator never raises a reference-mismatch error.
    """
    length = len(seq)
    fractions = [0.50, 0.30, 0.20, 0.10, 0.05, 0.02]
    variants: list[Variant] = []
    step = max(20, length // (n_sites + 1))
    pos = step
    for i in range(n_sites):
        ref = seq[pos - 1]
        alt = rng.choice([b for b in "ACGT" if b != ref])
        frac = fractions[i % len(fractions)]
        variants.append(Variant(pos, ref, alt, frac))
        pos += step
        if pos >= length - 10:
            break
    return variants


# Each rung is a genuinely different *decision rule*, not just a different
# likelihood.  The ladder isolates the marginal contribution of each component:
#
#   binomial              raw binomial LRT, uncorrected alpha = 0.05
#   beta_binomial         over-dispersed likelihood (rho > 0), uncorrected
#   beta_binomial_lrt     + Benjamini-Hochberg FDR control at q <= 0.05
#   beta_binomial_lrt_bayes  + calibrated posterior gate (>= 0.90)
#   full_policy           + strand-bias / homopolymer review filters and a
#                           minimum alternate-observation gate
#
# This is what makes the ablation informative: adding components can trade
# sensitivity at low allele fractions against false positives in negative
# controls, and the table records that trade-off explicitly.

RUNG_RULES: dict[str, dict[str, Any]] = {
    "binomial": {
        "model": "binomial", "rho": None, "fdr": False, "bayes": False,
        "review": False, "min_alt": 2,
    },
    "beta_binomial": {
        "model": "beta_binomial", "rho": 0.001, "fdr": False, "bayes": False,
        "review": False, "min_alt": 2,
    },
    "beta_binomial_lrt": {
        "model": "beta_binomial", "rho": 0.001, "fdr": True, "bayes": False,
        "review": False, "min_alt": 2,
    },
    "beta_binomial_lrt_bayes": {
        "model": "beta_binomial", "rho": 0.001, "fdr": True, "bayes": True,
        "review": False, "min_alt": 2,
    },
    "full_policy": {
        "model": "beta_binomial", "rho": 0.001, "fdr": True, "bayes": True,
        "review": True, "min_alt": 3,
    },
}


def _configure(cfg: PfVarXConfig, rung: str) -> PfVarXConfig:
    """Return a copy of ``cfg`` adjusted for one ablation rung."""
    rule = RUNG_RULES[rung]
    data = cfg.model_dump()
    stats = data["statistics"]
    stats["count_model"] = rule["model"]
    stats["beta_binomial_rho"] = rule["rho"]
    stats["lrt_null_calibration"] = "chi2_mixture"
    stats["error_rate"] = stats.get("error_rate") or 0.005
    stats["bayesian_enabled"] = bool(rule["bayes"])
    if rule["bayes"]:
        stats["posterior_threshold_if_validated"] = 0.90
        stats["bayesian_prior_h1"] = 0.01
    return PfVarXConfig.model_validate(data)


def _accept(call, rung: str) -> bool:
    """Apply the rung-specific acceptance rule to one candidate call."""
    rule = RUNG_RULES[rung]
    ev = call.evidence
    if ev.alt_count < rule["min_alt"]:
        return False
    p = call.lrt_pvalue_adjusted if rule["fdr"] else call.lrt_pvalue
    if p is None or p > 0.05:
        return False
    if rule["bayes"] and (call.posterior_h1 is None or call.posterior_h1 < 0.90):
        return False
    return not (
        rule["review"]
        and ("StrandBiasReview" in call.filters or "HomopolymerIndelReview" in call.filters)
    )


def _evaluate_replicate(
    cfg: PfVarXConfig,
    *,
    rung: str,
    work_dir: Path,
    replicate: int,
    coverage: int,
    seed: int,
) -> list[SiteOutcome]:
    rng = random.Random(seed + replicate)
    length = 2000
    data_dir = work_dir / f"rep{replicate}"
    data_dir.mkdir(parents=True, exist_ok=True)
    from pfvarx.synthetic import write_reference

    ref_path = str(data_dir / "reference.fa")
    write_reference(ref_path, "Pf3D7_01_v3", length, seed=seed + replicate)
    import pysam

    with pysam.FastaFile(ref_path) as fa:
        seq = fa.fetch("Pf3D7_01_v3").upper()
    variants = _design_variants(seq, rng, n_sites=6)
    ds = simulate_dataset(
        str(data_dir), length=length, coverage=coverage, variants=variants, seed=seed + replicate
    )
    run_cfg = _configure(cfg, rung)
    res = call_sample(
        ds.bam_path,
        run_cfg,
        sample_id=f"rep{replicate}",
        fasta_path=ds.reference_path,
        output_dir=str(work_dir / "results" / rung),
    )
    truth = {(v.pos, v.ref, v.alt): v.fraction for v in variants}
    accepted = {
        (c.evidence.pos, c.evidence.ref, c.evidence.alt): c
        for c in res.calls
        if c.evidence.alt_count > 0 and _accept(c, rung)
    }
    outcomes: list[SiteOutcome] = []
    for (pos, ref, alt), frac in truth.items():
        hit = accepted.get((pos, ref, alt))
        raw = next(
            (c for c in res.calls
             if (c.evidence.pos, c.evidence.ref, c.evidence.alt) == (pos, ref, alt)),
            None,
        )
        outcomes.append(
            SiteOutcome(
                chrom=ds.contig,
                pos=pos,
                ref=ref,
                alt=alt,
                truth_fraction=frac,
                called=hit is not None,
                decision="ACCEPTED" if hit is not None else "REJECTED",
                alt_fraction=raw.evidence.alt_fraction if raw else None,
                lrt_pvalue_adjusted=raw.lrt_pvalue_adjusted if raw else None,
                is_truth=True,
            )
        )
    # Negative controls: accepted non-truth candidates are false positives.
    for key, call in accepted.items():
        if key in truth:
            continue
        pos, ref, alt = key
        outcomes.append(
            SiteOutcome(
                chrom=call.evidence.chrom,
                pos=pos,
                ref=ref,
                alt=alt,
                truth_fraction=0.0,
                called=True,
                decision="FALSE_POSITIVE",
                alt_fraction=call.evidence.alt_fraction,
                lrt_pvalue_adjusted=call.lrt_pvalue_adjusted,
                is_truth=False,
            )
        )
    return outcomes


def run_benchmark(
    cfg: PfVarXConfig,
    *,
    out_dir: str | Path = "benchmarks/results",
    replicates: int = 5,
    coverage: int = 60,
    ablation: bool = True,
    seed: int = 11,
) -> dict[str, Any]:
    """Run the benchmark across the ablation ladder and persist results."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    work = out / "_work"
    work.mkdir(parents=True, exist_ok=True)

    rungs = ABLATION_LADDER if ablation else ["full_policy"]
    results: list[BenchmarkResult] = []

    for rung in rungs:
        bench = BenchmarkResult(config_label=rung, replicates=replicates, coverage=coverage)
        frac_stats: dict[str, list[bool]] = {}
        for rep in range(replicates):
            outcomes = _evaluate_replicate(
                cfg, rung=rung, work_dir=work, replicate=rep, coverage=coverage, seed=seed
            )
            for o in outcomes:
                if o.is_truth:
                    bench.truth_sites += 1
                    if o.called:
                        bench.true_positives += 1
                        bench.called_sites += 1
                    else:
                        bench.false_negatives += 1
                    key = f"{o.truth_fraction:.2f}"
                    frac_stats.setdefault(key, []).append(o.called)
                elif o.called:
                    bench.false_positives += 1
                    bench.called_sites += 1
                bench.outcomes.append(asdict(o))
        for key, hits in frac_stats.items():
            n = len(hits)
            k = sum(hits)
            lo, hi = wilson_interval(k, n)
            bench.per_fraction[key] = {
                "n": n,
                "detected": k,
                "sensitivity": k / n if n else 0.0,
                "sensitivity_ci": [lo, hi],
            }
        results.append(bench)

    payload = {
        "banner": "RESEARCH USE ONLY — synthetic benchmark, not clinical validation",
        "coverage": coverage,
        "replicates": replicates,
        "seed": seed,
        "ablation_ladder": rungs,
        "results": [r.as_dict() for r in results],
    }
    (out / "benchmark_results.json").write_text(json.dumps(payload, indent=2, default=str))

    # Flat CSV summary
    lines = ["config_label,precision,recall,f1,true_positives,false_positives,false_negatives"]
    for r in results:
        m = r.metrics()
        lines.append(
            f"{r.config_label},{m['precision']:.4f},{m['recall']:.4f},{m['f1']:.4f},"
            f"{m['true_positives']},{m['false_positives']},{m['false_negatives']}"
        )
    (out / "benchmark_summary.csv").write_text("\n".join(lines) + "\n")
    return payload
