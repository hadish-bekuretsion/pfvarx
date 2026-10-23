"""Input-detection and checksum tests."""

from __future__ import annotations

import gzip

import pytest

from pfvarx.errors import InputValidationError
from pfvarx.io.checksums import provenance, sha256_file, write_provenance
from pfvarx.io.validators import detect_type, validate_vcf


def test_detect_fasta(synthetic):
    assert detect_type(synthetic.reference_path) == "fasta"


def test_detect_fastq(synthetic):
    assert detect_type(synthetic.fastq_r1) == "fastq"


def test_detect_fastq_gz(tmp_path):
    p = tmp_path / "r1.fastq.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("@r1\nACGT\n+\nIIII\n")
    assert detect_type(p) == "fastq"


def test_detect_bam(synthetic):
    assert detect_type(synthetic.bam_path) == "bam"


def test_detect_vcf(synthetic):
    assert detect_type(synthetic.truth_vcf) == "vcf"


def test_detect_gvcf(tmp_path):
    p = tmp_path / "s.g.vcf"
    p.write_text(
        "##fileformat=VCFv4.2\n"
        "##ALT=<ID=NON_REF,Description=\"ref\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    assert detect_type(p) == "gvcf"


def test_detect_missing_raises(tmp_path):
    with pytest.raises(InputValidationError):
        detect_type(tmp_path / "nope.fa")


def test_detect_empty_raises(tmp_path):
    p = tmp_path / "empty.fa"
    p.write_text("")
    with pytest.raises(InputValidationError):
        detect_type(p)


def test_sha256_stable(synthetic):
    a = sha256_file(synthetic.reference_path)
    b = sha256_file(synthetic.reference_path)
    assert a == b
    assert len(a) == 64


def test_provenance_and_write(synthetic, tmp_path):
    rec = provenance(synthetic.reference_path, "reference")
    assert rec.role == "reference"
    assert rec.sha256
    out = write_provenance([rec], tmp_path / "prov.json")
    assert out.exists()


def test_validate_vcf_compares_full_ref_allele(tmp_path):
    import pysam

    reference = tmp_path / "reference.fa"
    reference.write_text(">chr1\nAACCGGTT\n")
    pysam.faidx(str(reference))
    vcf = tmp_path / "multibase_ref.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=8>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t2\t.\tAG\tT\t.\tPASS\t.\n"
    )

    result = validate_vcf(vcf, reference_fasta=reference)

    assert result.metadata["ref_mismatches"] == 1
    assert any(issue.code == "VCF_REF_MISMATCH" for issue in result.issues)


def test_validate_vcf_uses_header_contig_order(tmp_path):
    vcf = tmp_path / "header_order.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=100>\n"
        "##contig=<ID=chr2,length=100>\n"
        "##contig=<ID=chr10,length=100>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t2\t.\tA\tT\t.\tPASS\t.\n"
        "chr2\t2\t.\tA\tT\t.\tPASS\t.\n"
        "chr10\t2\t.\tA\tT\t.\tPASS\t.\n"
    )

    result = validate_vcf(vcf, check_ref=False)

    assert result.metadata["out_of_order"] == 0
    assert not any(issue.code == "VCF_UNSORTED" for issue in result.issues)
