"""Internal consequence prediction for SNPs and small INDELs.

This is a transparent, dependency-free fallback used when SnpEff is not
installed.  It is not a substitute for a curated annotation database and is
reported as ``INTERNAL`` provenance so downstream users can tell the
difference.
"""

from __future__ import annotations

from dataclasses import dataclass

from pfvarx.annotation.gff import CDS, GeneModel, translate

SEVERITY_ORDER = [
    "intergenic_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "intron_variant",
    "synonymous_variant",
    "missense_variant",
    "inframe_insertion",
    "inframe_deletion",
    "stop_lost",
    "start_lost",
    "stop_gained",
    "frameshift_variant",
    "coding_sequence_variant",
]


@dataclass
class Consequence:
    gene_id: str | None
    gene_name: str | None
    transcript_biotype: str
    consequence: str
    amino_acid_change: str | None
    codon_change: str | None
    protein_position: int | None
    severity_rank: int
    provenance: str = "INTERNAL"

    def as_dict(self) -> dict:
        return {
            "gene_id": self.gene_id or ".",
            "gene_name": self.gene_name or ".",
            "transcript_biotype": self.transcript_biotype,
            "consequence": self.consequence,
            "amino_acid_change": self.amino_acid_change or ".",
            "codon_change": self.codon_change or ".",
            "protein_position": self.protein_position if self.protein_position is not None else ".",
            "severity_rank": self.severity_rank,
            "provenance": self.provenance,
        }


def _rank(term: str) -> int:
    return SEVERITY_ORDER.index(term) if term in SEVERITY_ORDER else len(SEVERITY_ORDER)


def _annotate_cds_snp(cds: CDS, cdna: str, cds_idx: int, ref: str, alt: str) -> Consequence:
    # A variant in a leading partial codon has no complete codon position to
    # report. Do not invent an amino-acid change for it.
    cds_idx -= cds.phase
    if cds_idx < 0:
        term = "coding_sequence_variant"
        return Consequence(
            gene_id=cds.gene_id,
            gene_name=cds.gene_name,
            transcript_biotype="protein_coding",
            consequence=term,
            amino_acid_change=None,
            codon_change=None,
            protein_position=None,
            severity_rank=_rank(term),
        )
    if cds.strand == "-":
        complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
        ref = complement.get(ref.upper(), "N")
        alt = complement.get(alt.upper(), "N")
    codon_no = cds_idx // 3
    offset = cds_idx % 3
    ref_codon = cdna[codon_no * 3 : codon_no * 3 + 3]
    alt_codon = list(ref_codon)
    alt_codon[offset] = alt
    alt_codon = "".join(alt_codon)
    ref_aa = translate(ref_codon)
    alt_aa = translate(alt_codon)
    protein_pos = codon_no + 1

    if ref_aa == alt_aa:
        term = "synonymous_variant"
    elif ref_aa == "M" and protein_pos == 1:
        term = "start_lost"
    elif alt_aa == "*":
        term = "stop_gained"
    elif ref_aa == "*":
        term = "stop_lost"
    else:
        term = "missense_variant"

    aa_change = None
    if ref_aa != alt_aa:
        aa_change = f"p.{ref_aa}{protein_pos}{alt_aa}"
    return Consequence(
        gene_id=cds.gene_id,
        gene_name=cds.gene_name,
        transcript_biotype="protein_coding",
        consequence=term,
        amino_acid_change=aa_change,
        codon_change=f"{ref_codon}>{alt_codon}",
        protein_position=protein_pos,
        severity_rank=_rank(term),
    )


def annotate_variant(
    model: GeneModel,
    fasta,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
) -> list[Consequence]:
    """Return consequences for a variant, one per overlapping gene."""
    genes = model.genes_at(chrom, pos)
    if not genes:
        return [
            Consequence(
                gene_id=None,
                gene_name=None,
                transcript_biotype="intergenic",
                consequence="intergenic_variant",
                amino_acid_change=None,
                codon_change=None,
                protein_position=None,
                severity_rank=_rank("intergenic_variant"),
            )
        ]

    results: list[Consequence] = []
    for cds in genes:
        cdna = cds.cdna(fasta)
        genomic_idx = cds.genomic_to_cds_index(pos)
        if genomic_idx is None:
            continue
        idx = genomic_idx - cds.phase
        if len(ref) == 1 and len(alt) == 1:
            results.append(_annotate_cds_snp(cds, cdna, genomic_idx, ref, alt))
        else:
            delta = len(alt) - len(ref)
            if delta % 3 != 0:
                term = "frameshift_variant"
            elif delta > 0:
                term = "inframe_insertion"
            else:
                term = "inframe_deletion"
            protein_pos = idx // 3 + 1
            results.append(
                Consequence(
                    gene_id=cds.gene_id,
                    gene_name=cds.gene_name,
                    transcript_biotype="protein_coding",
                    consequence=term,
                    amino_acid_change=f"p.{protein_pos}_{protein_pos + abs(delta) // 3}del"
                    if delta < 0
                    else f"p.{protein_pos}_{protein_pos + delta // 3}ins",
                    codon_change=None,
                    protein_position=protein_pos,
                    severity_rank=_rank(term),
                )
            )
    return results
