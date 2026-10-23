"""Caller tests for complete hypothesis-family accounting."""

from __future__ import annotations

from pfvarx.calling.caller import _hypothesis_family_size


def test_family_size_bounds_all_alleles_in_requested_span(synthetic, cfg):
    family = _hypothesis_family_size(
        synthetic.bam_path,
        cfg,
        chrom=synthetic.contig,
        start=100,
        end=199,
    )
    max_indel = cfg.candidate.max_small_indel_bp
    per_anchor = 3 + 2 * sum(4**length for length in range(1, max_indel + 1))
    assert family == 100 * per_anchor
    assert family > 100 * 3  # substitutions alone; INDEL universe is also bounded
