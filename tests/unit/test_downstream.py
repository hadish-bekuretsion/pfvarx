"""Downstream analysis tests: diversity, mixed infection, relatedness, trees."""

from __future__ import annotations

from pfvarx.calling.policy import (
    BIN_INTERMEDIATE,
    BIN_PREDOMINANT,
    SiteCall,
    SiteEvidence,
)
from pfvarx.downstream.diversity import diversity_summary
from pfvarx.downstream.mixed_infection import mixed_infection_summary
from pfvarx.downstream.phylogenetics import newick
from pfvarx.downstream.relatedness import pairwise_ibs, relatedness_matrix


def _call(pos: int, ref: int, alt: int, *, bin_: str = BIN_PREDOMINANT) -> SiteCall:
    depth = ref + alt
    ev = SiteEvidence(
        sample_id="s", chrom="c", pos=pos, ref="A", alt="T", var_type="SNP",
        ref_count=ref, alt_count=alt, depth=depth, informative_depth=depth,
        alt_fraction=(alt / depth) if depth else 0.0,
    )
    return SiteCall(evidence=ev, decision="ALT_DETECTED_EXPLORATORY", bin=bin_)


def test_diversity_summary_basic():
    calls = [_call(1, 50, 50), _call(2, 100, 0), _call(3, 80, 20)]
    s = diversity_summary(calls, sample_id="S1")
    assert s.sample_id == "S1"
    assert s.n_sites == 3
    assert s.n_polymorphic == 2
    assert 0 < s.nucleotide_diversity <= 0.5
    assert 0 <= s.proportion_polymorphic <= 1


def test_diversity_empty():
    s = diversity_summary([], sample_id="S1")
    assert s.n_sites == 0
    assert s.nucleotide_diversity == 0.0


def test_mixed_infection_likely():
    calls = [_call(i, 50, 50, bin_=BIN_INTERMEDIATE) for i in range(1, 5)]
    m = mixed_infection_summary(calls, sample_id="S1")
    assert m.mixed_infection_flag == "MIXED_INFECTION_LIKELY"
    # No validated COI model is implemented, so no numeric estimate is reported.
    assert m.coi_point_estimate is None
    assert "no numeric COI" in m.method


def test_mixed_infection_none():
    calls = [_call(1, 90, 10, bin_=BIN_PREDOMINANT)]
    m = mixed_infection_summary(calls, sample_id="S1")
    assert m.mixed_infection_flag == "NO_EVIDENCE_OF_MIXED_INFECTION"


def test_pairwise_ibs_identical_samples():
    a = [_call(1, 90, 10), _call(2, 90, 10)]
    b = [_call(1, 90, 10), _call(2, 90, 10)]
    r = pairwise_ibs("A", a, "B", b)
    assert r.n_shared_sites == 2
    assert r.ibs2 == 1.0
    assert r.relatedness_estimate == 1.0


def test_pairwise_ibs_disjoint_genotypes():
    a = [_call(1, 90, 10)]
    b = [_call(1, 10, 90)]
    r = pairwise_ibs("A", a, "B", b)
    assert r.n_discordant == 1
    assert r.ibs2 == 0.0


def test_relatedness_matrix_and_newick():
    samples = {
        "A": [_call(1, 90, 10), _call(2, 90, 10)],
        "B": [_call(1, 90, 10), _call(2, 90, 10)],
        "C": [_call(1, 10, 90), _call(2, 10, 90)],
    }
    results = relatedness_matrix(samples)
    assert len(results) == 3  # 3 choose 2
    tree = newick(sorted(samples), results)
    assert tree.endswith(";")
    assert "A" in tree and "B" in tree and "C" in tree
