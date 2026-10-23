"""Synthetic P. falciparum-like data generator for tests and benchmarks.

Produces a small reference, paired-end FASTQ, coordinate-sorted BAM and a truth
VCF with a known set of SNPs and small INDELs at controlled allele fractions.
The generator is deterministic given a seed so that benchmark splits are
reproducible.  It is deliberately simple (no quality-score recalibration, no
real error model) and is *only* suitable for software testing and for
demonstrating the analysis mechanics -- never for scientific validation.
"""

from __future__ import annotations

import gzip
import random
from dataclasses import dataclass, field
from pathlib import Path

import pysam

BASES = "ACGT"
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


@dataclass
class Variant:
    pos: int  # 1-based
    ref: str
    alt: str
    fraction: float

    @property
    def var_type(self) -> str:
        return "SNP" if len(self.ref) == len(self.alt) == 1 else "INDEL"


@dataclass
class SyntheticDataset:
    reference_path: str
    fastq_r1: str
    fastq_r2: str
    bam_path: str
    truth_vcf: str
    variants: list[Variant] = field(default_factory=list)
    contig: str = "Pf3D7_01_v3"
    length: int = 0


def _random_sequence(n: int, rng: random.Random) -> str:
    return "".join(rng.choice(BASES) for _ in range(n))


def _revcomp(seq: str) -> str:
    return "".join(COMPLEMENT[b] for b in reversed(seq))


def write_reference(path: str, contig: str, length: int, seed: int = 1) -> str:
    rng = random.Random(seed)
    seq = _random_sequence(length, rng)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write(f">{contig}\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i : i + 60] + "\n")
    return path


def _build_haplotype(ref: str, variants: list[Variant]) -> tuple[str, list[int]]:
    """Return (haplotype, ref_to_hap) where ref_to_hap[i] is the hap index of
    reference base i (0-based).  Deleted bases map to the anchor position;
    inserted bases are appended after the anchor."""
    by_pos = {v.pos - 1: v for v in variants}
    hap: list[str] = []
    ref_to_hap = [0] * len(ref)
    i = 0
    while i < len(ref):
        v = by_pos.get(i)
        if v is None:
            ref_to_hap[i] = len(hap)
            hap.append(ref[i])
            i += 1
            continue
        # anchor position
        ref_to_hap[i] = len(hap)
        if v.var_type == "SNP":
            hap.append(v.alt)
            i += 1
        elif len(v.alt) > len(v.ref):  # insertion
            hap.append(v.ref[0])
            hap.extend(v.alt[1:])
            i += 1
        else:  # deletion
            hap.append(v.ref[0])
            for j in range(i + 1, i + len(v.ref)):
                ref_to_hap[j] = len(hap) - 1
            i += len(v.ref)
    return "".join(hap), ref_to_hap


def simulate_dataset(
    out_dir: str,
    *,
    contig: str = "Pf3D7_01_v3",
    length: int = 2000,
    variants: list[Variant] | None = None,
    coverage: int = 60,
    read_len: int = 100,
    seed: int = 7,
    write_fastq: bool = True,
) -> SyntheticDataset:
    """Create reference, reads, BAM and truth VCF in ``out_dir``."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(seed)

    ref_path = str(out / "reference.fa")
    write_reference(ref_path, contig, length, seed=seed)
    with pysam.FastaFile(ref_path) as fa:
        ref_seq = fa.fetch(contig).upper()

    if variants is None:
        variants = []
    for v in variants:
        actual = ref_seq[v.pos - 1 : v.pos - 1 + len(v.ref)]
        if actual != v.ref:
            raise ValueError(f"variant at {v.pos} REF={v.ref} but reference has {actual}")


    n_fragments = max(1, int(length * coverage / read_len))
    reads: list[tuple[str, str, str, int]] = []  # name, seq, qual, ref_start
    counter = 0

    def _emit(start: int, seq: str) -> None:
        nonlocal counter
        if len(seq) < read_len:
            return
        chars = list(seq)
        for j, b in enumerate(chars):
            if rng.random() < 0.005:
                chars[j] = rng.choice([x for x in BASES if x != b])
        qual = "".join(chr(53 + rng.randint(0, 10)) for _ in chars)
        reads.append((f"read{counter}", "".join(chars), qual, start))
        counter += 1

    # Alternate-supporting fragments, allocated per variant at that variant's
    # own allele fraction, so that each truth site has the intended AF.
    n_alt_total = 0
    for v in variants:
        hap_v, map_v = _build_haplotype(ref_seq, [v])
        n_alt_v = max(1, int(n_fragments * v.fraction))
        n_alt_total += n_alt_v
        for _ in range(n_alt_v):
            start = rng.randint(0, max(0, length - read_len))
            hstart = map_v[start]
            _emit(start, hap_v[hstart : hstart + read_len])

    n_ref = max(0, n_fragments - n_alt_total)
    for _ in range(n_ref):
        start = rng.randint(0, max(0, length - read_len))
        _emit(start, ref_seq[start : start + read_len])

    r1_path = str(out / "sample_R1.fastq.gz")
    r2_path = str(out / "sample_R2.fastq.gz")
    if write_fastq:
        with gzip.open(r1_path, "wt") as f1, gzip.open(r2_path, "wt") as f2:
            for name, seq, qual, _ in reads:
                half = len(seq) // 2
                s1, s2 = seq[:half], _revcomp(seq[half:])
                q1, q2 = qual[:half], qual[half:]
                f1.write(f"@{name}/1\n{s1}\n+\n{q1}\n")
                f2.write(f"@{name}/2\n{s2}\n+\n{q2}\n")

    bam_path = str(out / "sample.bam")
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": contig, "LN": length}],
        "RG": [{"ID": "rg1", "SM": "sample"}],
    }
    tmp_bam = str(out / "sample.unsorted.bam")
    with pysam.AlignmentFile(tmp_bam, "wb", header=header) as bf:
        for name, seq, qual, start in reads:
            a = pysam.AlignedSegment()
            a.query_name = name
            a.query_sequence = seq
            a.flag = 0
            a.reference_id = 0
            a.reference_start = start
            a.mapping_quality = 60
            a.cigar = ((0, len(seq)),)
            a.query_qualities = pysam.qualitystring_to_array(qual)
            a.set_tag("RG", "rg1")
            bf.write(a)
    pysam.sort("-o", bam_path, tmp_bam)
    pysam.index(bam_path)
    Path(tmp_bam).unlink(missing_ok=True)

    truth_path = str(out / "truth.vcf")
    with open(truth_path, "w") as fh:
        fh.write("##fileformat=VCFv4.3\n")
        fh.write(f"##contig=<ID={contig},length={length}>\n")
        fh.write('##INFO=<ID=AF,Number=1,Type=Float,Description="Simulated allele fraction">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for v in sorted(variants, key=lambda x: x.pos):
            fh.write(f"{contig}\t{v.pos}\t.\t{v.ref}\t{v.alt}\t.\tPASS\tAF={v.fraction}\n")

    return SyntheticDataset(
        reference_path=ref_path,
        fastq_r1=r1_path,
        fastq_r2=r2_path,
        bam_path=bam_path,
        truth_vcf=truth_path,
        variants=variants,
        contig=contig,
        length=length,
    )
