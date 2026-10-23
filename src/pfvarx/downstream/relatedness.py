"""Pairwise relatedness from identity-by-state sharing.

A transparent IBS-based estimate is always available; hmmIBD is wrapped for
production use when installed.  Neither is a substitute for a validated
transmission-inference workflow.
"""

from __future__ import annotations

from dataclasses import dataclass

from pfvarx.calling.policy import SiteCall


@dataclass
class RelatednessResult:
    sample_a: str
    sample_b: str
    n_shared_sites: int
    n_concordant: int
    n_discordant: int
    ibs0: float
    ibs1: float
    ibs2: float
    relatedness_estimate: float
    method: str

    def as_dict(self) -> dict:
        return self.__dict__.copy()


def _genotype(call: SiteCall) -> int:
    """Majority genotype: 0 = reference, 1 = alternate."""
    return 1 if call.evidence.alt_count > call.evidence.ref_count else 0


def pairwise_ibs(
    sample_a: str,
    calls_a: list[SiteCall],
    sample_b: str,
    calls_b: list[SiteCall],
    *,
    min_depth: int = 10,
) -> RelatednessResult:
    """Compute IBS sharing over sites present in both samples."""
    a = {(c.evidence.chrom, c.evidence.pos): c for c in calls_a if c.evidence.informative_depth >= min_depth}
    b = {(c.evidence.chrom, c.evidence.pos): c for c in calls_b if c.evidence.informative_depth >= min_depth}
    shared = sorted(set(a) & set(b))
    concordant = discordant = 0
    for key in shared:
        if _genotype(a[key]) == _genotype(b[key]):
            concordant += 1
        else:
            discordant += 1
    n = len(shared)
    # IBS2 = identical genotype, IBS0 = opposite homozygotes (undefined for
    # haploid data, so reported as 0), IBS1 = one shared allele.
    ibs2 = concordant / n if n else 0.0
    ibs0 = 0.0
    ibs1 = discordant / n if n else 0.0
    relatedness = ibs2 if n else 0.0
    return RelatednessResult(
        sample_a=sample_a,
        sample_b=sample_b,
        n_shared_sites=n,
        n_concordant=concordant,
        n_discordant=discordant,
        ibs0=ibs0,
        ibs1=ibs1,
        ibs2=ibs2,
        relatedness_estimate=relatedness,
        method="haploid IBS sharing (PROVISIONAL); hmmIBD available as optional wrapper",
    )


def relatedness_matrix(samples: dict[str, list[SiteCall]], *, min_depth: int = 10) -> list[RelatednessResult]:
    names = sorted(samples)
    out = []
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            out.append(pairwise_ibs(names[i], samples[names[i]], names[j], samples[names[j]], min_depth=min_depth))
    return out
