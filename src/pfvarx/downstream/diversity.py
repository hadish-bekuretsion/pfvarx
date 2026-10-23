"""Within- and between-sample diversity statistics.

All estimators are implemented directly from allele counts so that they are
transparent, testable and independent of optional third-party libraries.
Formulas follow Nei and Li (1979) and Nei (1987); see METHODS documentation.
"""

from __future__ import annotations

from dataclasses import dataclass

from pfvarx.calling.policy import SiteCall


@dataclass
class DiversitySummary:
    sample_id: str
    n_sites: int
    n_polymorphic: int
    proportion_polymorphic: float
    mean_alt_fraction: float
    expected_heterozygosity: float
    nucleotide_diversity: float
    watterson_theta: float | None
    segregating_sites: int

    def as_dict(self) -> dict:
        return self.__dict__.copy()


def _harmonic(n: int) -> float:
    return sum(1.0 / i for i in range(1, n + 1)) if n > 0 else 0.0


def diversity_summary(
    calls: list[SiteCall],
    *,
    sample_id: str = "sample",
    callable_bases: int | None = None,
) -> DiversitySummary:
    """Compute diversity statistics from called sites.

    ``expected_heterozygosity`` uses the unbiased estimator
    ``2p(1-p) * n/(n-1)`` per site, averaged over sites with depth >= 2.
    ``nucleotide_diversity`` is the mean per-site heterozygosity (pi).
    """
    usable = [c for c in calls if c.evidence.informative_depth >= 2]
    n_sites = len(usable)
    polymorphic = [c for c in usable if c.evidence.alt_count > 0 and c.evidence.ref_count > 0]
    n_poly = len(polymorphic)

    het_values = []
    for c in usable:
        n = c.evidence.informative_depth
        p = c.evidence.alt_fraction
        if n >= 2 and p is not None:
            het_values.append(2 * p * (1 - p) * n / (n - 1))
    mean_het = sum(het_values) / len(het_values) if het_values else 0.0

    # Watterson's theta is defined in terms of the number of sequences (n) and
    # the number of segregating sites, not the number of assayed sites.  A
    # single haploid sample provides n=1, for which theta is undefined, so no
    # value is reported rather than a misapplied estimate.
    theta = None

    alt_fractions = [c.evidence.alt_fraction for c in usable if c.evidence.alt_fraction is not None]
    return DiversitySummary(
        sample_id=sample_id,
        n_sites=n_sites,
        n_polymorphic=n_poly,
        proportion_polymorphic=(n_poly / n_sites) if n_sites else 0.0,
        mean_alt_fraction=(sum(alt_fractions) / len(alt_fractions)) if alt_fractions else 0.0,
        expected_heterozygosity=mean_het,
        nucleotide_diversity=mean_het,
        watterson_theta=theta,
        segregating_sites=n_poly,
    )
