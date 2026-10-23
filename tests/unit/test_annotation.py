"""Annotation tests: translation, consequence calling, resistance catalogue."""

from __future__ import annotations

import pysam

from pfvarx.annotation.catalogue import DEFAULT_CATALOGUE, load_catalogue
from pfvarx.annotation.consequence import annotate_variant
from pfvarx.annotation.gff import parse_gff3, revcomp, translate


def test_revcomp_is_involution():
    s = "ACGTACGT"
    assert revcomp(revcomp(s)) == s
    assert revcomp("AACCGG") == "CCGGTT"


def test_translate_standard_codons():
    assert translate("ATG") == "M"
    assert translate("TAA") == "*"
    assert translate("TTT") == "F"
    assert translate("GGG") == "G"


def test_translate_stops_at_stop_codon():
    assert translate("ATGTAAATG") == "M*M"


def test_default_catalogue_is_bundled_with_package():
    assert DEFAULT_CATALOGUE.name == "pf_resistance.yaml"
    assert DEFAULT_CATALOGUE.is_file()
    assert load_catalogue().markers


def test_resistance_catalogue_loads_and_has_markers():
    cat = load_catalogue()
    assert cat is not None
    assert len(cat.by_gene("PF3D7_0709000")) > 0
    # PfCRT K76T is a canonical, well-established marker.
    assert cat.by_gene("PF3D7_0709000")


def test_resistance_catalogue_match_returns_marker():
    cat = load_catalogue()
    hits = cat.match("PF3D7_0709000", "K76T")
    assert hits
    assert hits[0].drug
    assert cat.match("PF3D7_0709000", "p.K76T")


def test_resistance_catalogue_no_false_match():
    cat = load_catalogue()
    assert cat.match("PF3D7_0709000", "K999Z") == []


def _fasta(tmp_path, sequence: str) -> pysam.FastaFile:
    path = tmp_path / "ref.fa"
    path.write_text(f">Pf3D7\n{sequence}\n")
    pysam.faidx(str(path))
    return pysam.FastaFile(str(path))


def test_transcript_parent_resolves_to_gene_and_plus_cds_segments(tmp_path):
    fasta = _fasta(tmp_path, "ATGG" + "N" * 5 + "CTTAA")
    gff = tmp_path / "plus.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "Pf3D7\ttest\tgene\t1\t14\t.\t+\t.\tID=g1;Name=GENE1\n"
        "Pf3D7\ttest\tmRNA\t1\t14\t.\t+\t.\tID=tx1;Parent=g1\n"
        "Pf3D7\ttest\tCDS\t10\t14\t.\t+\t2\tParent=tx1\n"
        "Pf3D7\ttest\tCDS\t1\t4\t.\t+\t0\tParent=tx1\n"
    )
    model = parse_gff3(str(gff))
    cds = model.cds_by_gene["g1"]
    assert cds.transcript_id == "tx1"
    assert cds.gene_name == "GENE1"
    assert cds.segments == [(1, 4), (10, 14)]
    assert cds.cdna(fasta) == "ATGGCTTAA"
    result = annotate_variant(model, fasta, "Pf3D7", 10, "C", "T")[0]
    assert result.amino_acid_change == "p.A2V"
    fasta.close()


def test_minus_strand_multisegment_cds_and_phase_are_transcription_ordered(tmp_path):
    sequence = "N" * 9 + "TTAA" + "N" * 6 + "GCCAT"
    fasta = _fasta(tmp_path, sequence)
    gff = tmp_path / "minus.gff3"
    gff.write_text(
        "Pf3D7\ttest\tgene\t10\t24\t.\t-\t.\tID=gminus;Name=MINUS\n"
        "Pf3D7\ttest\tmRNA\t10\t24\t.\t-\t.\tID=txm;Parent=gminus\n"
        "Pf3D7\ttest\tCDS\t10\t13\t.\t-\t1\tParent=txm\n"
        "Pf3D7\ttest\tCDS\t20\t24\t.\t-\t0\tParent=txm\n"
    )
    model = parse_gff3(str(gff))
    cds = model.cds_by_gene["gminus"]
    assert cds.segments == [(20, 24), (10, 13)]
    assert cds.cdna(fasta) == "ATGGCTTAA"
    result = annotate_variant(model, fasta, "Pf3D7", 20, "G", "A")[0]
    assert result.amino_acid_change == "p.A2V"
    fasta.close()


def test_nonzero_initial_phase_skips_leading_partial_codon(tmp_path):
    fasta = _fasta(tmp_path, "CATGGCTTAA")
    gff = tmp_path / "partial.gff3"
    gff.write_text(
        "Pf3D7\ttest\tgene\t1\t10\t.\t+\t.\tID=gpartial\n"
        "Pf3D7\ttest\tmRNA\t1\t10\t.\t+\t.\tID=txp;Parent=gpartial\n"
        "Pf3D7\ttest\tCDS\t1\t10\t.\t+\t1\tParent=txp\n"
    )
    model = parse_gff3(str(gff))
    cds = model.cds_by_gene["gpartial"]
    assert cds.phase == 1
    assert cds.cdna(fasta) == "ATGGCTTAA"
    result = annotate_variant(model, fasta, "Pf3D7", 6, "C", "T")[0]
    assert result.amino_acid_change == "p.A2V"
    fasta.close()


def test_multiple_transcripts_select_lexicographically_smallest_id(tmp_path):
    fasta = _fasta(tmp_path, "ATGGCTTAA" + "N" * 10 + "ATGTTTTAA")
    gff = tmp_path / "multi.gff3"
    gff.write_text(
        "Pf3D7\ttest\tgene\t1\t28\t.\t+\t.\tID=gmulti\n"
        "Pf3D7\ttest\tmRNA\t20\t28\t.\t+\t.\tID=txZ;Parent=gmulti\n"
        "Pf3D7\ttest\tmRNA\t1\t9\t.\t+\t.\tID=txA;Parent=gmulti\n"
        "Pf3D7\ttest\tCDS\t20\t28\t.\t+\t0\tParent=txZ\n"
        "Pf3D7\ttest\tCDS\t1\t9\t.\t+\t0\tParent=txA\n"
    )
    model = parse_gff3(str(gff))
    assert model.cds_by_gene["gmulti"].transcript_id == "txA"
    assert model.cds_by_gene["gmulti"].cdna(fasta) == "ATGGCTTAA"
    fasta.close()


def test_generic_nonsynonymous_rules_are_only_low_confidence_candidates():
    cat = load_catalogue()
    # This substitution is not an explicit catalogue variant. It only matches
    # broad textual placeholders, which must be visibly downgraded.
    hits = cat.match("PF3D7_1343700", "p.K999R")
    assert hits
    assert all(hit.evidence_level == "low_confidence_candidate" for hit in hits)
    assert all(hit.as_dict()["evidence_level"] != "established" for hit in hits)
    assert cat.match("PF3D7_1343700", ".") == []
