"""Mixed-infection and complexity-of-infection summaries.

The platform is haploid, so true heterozygous genotypes are not expected.
Intermediate alternate fractions therefore indicate either mixed infection,
copy-number variation, mapping artefacts or sequencing error.  This module
summarises the evidence; it does not assert a biological cause.
"""

from __future__ import annotations

from dataclasses import dataclass

from pfvarx.calling.policy import (
    BIN_INTERMEDIATE,
    BIN_MINORITY,
    BIN_PREDOMINANT,
    SiteCall,
)


@dataclass
class MixedInfectionSummary:
    sample_id: str
    n_sites: int
    n_intermediate: int
    n_minority: int
    n_predominant: int
    max_minor_fraction: float
    mixed_infection_flag: str
    coi_point_estimate: float | None
    method: str

    def as_dict(self) -> dict:
        return self.__dict__.copy()


def mixed_infection_summary(calls: list[SiteCall], *, sample_id: str = "sample") -> MixedInfectionSummary:
    """Summarise intermediate-fraction evidence and give a provisional COI flag."""
    n_inter = sum(1 for c in calls if c.bin == BIN_INTERMEDIATE)
    n_minor = sum(1 for c in calls if c.bin == BIN_MINORITY)
    n_pred = sum(1 for c in calls if c.bin == BIN_PREDOMINANT)
    minor_fracs = [c.evidence.alt_fraction for c in calls if c.bin in (BIN_MINORITY, BIN_INTERMEDIATE)]
    max_minor = max(minor_fracs) if minor_fracs else 0.0

    if n_inter >= 3:
        flag = "MIXED_INFECTION_LIKELY"
    elif n_inter >= 1 or n_minor >= 3:
        flag = "MIXED_INFECTION_POSSIBLE"
    else:
        flag = "NO_EVIDENCE_OF_MIXED_INFECTION"

    # No validated COI model is implemented.  Reporting a numeric estimate from
    # the count of intermediate-fraction sites would be fabricated; the count
    # is exposed as an ordinal indicator and a numeric COI is left to a
    # validated method (e.g. DEploid).
    coi = None

    return MixedInfectionSummary(
        sample_id=sample_id,
        n_sites=len(calls),
        n_intermediate=n_inter,
        n_minority=n_minor,
        n_predominant=n_pred,
        max_minor_fraction=max_minor,
        mixed_infection_flag=flag,
        coi_point_estimate=coi,
        method="fraction-bin heuristic (PROVISIONAL); no numeric COI is reported; use DEploid for validated COI",
    )
