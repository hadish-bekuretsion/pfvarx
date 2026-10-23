"""Minimal GFF3 parser and CDS model for internal consequence prediction.

CDS ``Parent`` links are resolved through transcript features to their gene.
When several coding transcripts exist for one gene, the lexicographically
smallest transcript ID is selected.  This is deterministic, not a biological
canonical-transcript assertion; the internal annotator currently emits one
consequence per gene.  CDS phase is checked across segments and the initial
partial-codon phase is applied when reconstructing the translated sequence.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from urllib.parse import unquote

from pfvarx.errors import InputValidationError

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def revcomp(seq: str) -> str:
    return "".join(COMPLEMENT.get(b, "N") for b in reversed(seq.upper()))


def translate(seq: str) -> str:
    seq = seq.upper()
    return "".join(CODON_TABLE.get(seq[i : i + 3], "X") for i in range(0, len(seq) - 2, 3))


@dataclass
class CDS:
    gene_id: str
    gene_name: str
    chrom: str
    strand: str
    segments: list[tuple[int, int]]  # 1-based inclusive, in transcription order
    product: str = ""
    transcript_id: str = ""
    phase: int = 0

    def cdna(self, fasta) -> str:
        parts = []
        for start, end in self.segments:
            piece = fasta.fetch(self.chrom, start - 1, end).upper()
            parts.append(piece if self.strand == "+" else revcomp(piece))
        return "".join(parts)[self.phase :]

    def genomic_to_cds_index(self, pos: int) -> int | None:
        """Map a 1-based genomic position to a 0-based cDNA index."""
        offset = 0
        for start, end in self.segments:
            if start <= pos <= end:
                if self.strand == "+":
                    return offset + (pos - start)
                return offset + (end - pos)
            offset += end - start + 1
        return None


@dataclass
class GeneModel:
    cds_by_gene: dict[str, CDS] = field(default_factory=dict)
    genes_by_chrom: dict[str, list[CDS]] = field(default_factory=dict)

    def add(self, cds: CDS) -> None:
        # parse_gff3 resolves transcript choice before add(); assigning here
        # keeps the mapping contract simple for callers constructing models.
        self.cds_by_gene[cds.gene_id] = cds
        self.genes_by_chrom.setdefault(cds.chrom, []).append(cds)

    def genes_at(self, chrom: str, pos: int) -> list[CDS]:
        return [g for g in self.genes_by_chrom.get(chrom, []) if g.genomic_to_cds_index(pos) is not None]


def parse_gff3(path: str) -> GeneModel:
    """Parse genes/CDS, resolving CDS -> transcript -> gene relationships.

    Unresolvable or phase-inconsistent transcript models are omitted rather
    than guessed.  Direct ``CDS Parent=gene`` records remain supported.
    """
    genes: dict[str, dict] = {}
    transcripts: dict[str, dict] = {}
    raw_cds: list[tuple[list[str], int, int, str, str, int | None, str | None]] = []
    try:
        with open(path) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 9:
                    continue
                chrom, _src, ftype, start, end, _score, strand, phase, attrs = cols
                if ftype not in ("gene", "CDS", "mRNA", "transcript"):
                    continue
                attr = {}
                for item in attrs.split(";"):
                    if "=" in item:
                        k, v = item.split("=", 1)
                        attr[k.strip()] = unquote(v.strip())
                if ftype == "gene":
                    gid = attr.get("ID") or attr.get("gene_id")
                    if not gid:
                        continue
                    genes[gid] = {
                        "name": attr.get("Name", gid),
                        "product": attr.get("product", ""),
                        "chrom": chrom,
                        "strand": strand,
                    }
                elif ftype in ("mRNA", "transcript"):
                    tid = attr.get("ID")
                    if not tid:
                        continue
                    parents = [p for p in attr.get("Parent", "").split(",") if p]
                    # GFF3 records need not be ordered. Retain parent IDs and
                    # resolve against the complete gene table below.
                    transcripts[tid] = {
                        "parents": parents,
                        "name": attr.get("Name", tid),
                        "product": attr.get("product", ""),
                        "chrom": chrom,
                        "strand": strand,
                    }
                elif ftype == "CDS":
                    parents = [p for p in attr.get("Parent", "").split(",") if p]
                    gene_hint = attr.get("gene_id") or attr.get("gene")
                    try:
                        phase_value = int(phase) if phase in ("0", "1", "2") else None
                        start_value, end_value = int(start), int(end)
                    except ValueError:
                        phase_value = None
                        start_value, end_value = int(start), int(end)
                    raw_cds.append((parents, start_value, end_value, chrom, strand, phase_value, gene_hint))
    except FileNotFoundError as exc:
        raise InputValidationError(
            f"GFF3 file not found: {path}",
            cause=str(exc),
            action="Provide a valid GFF3 annotation or disable internal annotation.",
        ) from exc

    # Expand each CDS parent into transcript models. The intermediate key
    # preserves distinct isoforms until deterministic gene-level selection.
    grouped: dict[tuple[str, str], list[tuple[int, int, str, str, int | None]]] = {}
    transcript_meta: dict[tuple[str, str], dict] = {}
    for parents, start, end, chrom, strand, phase_value, gene_hint in raw_cds:
        resolved: list[tuple[str, str, dict]] = []
        for parent in parents:
            if parent in transcripts:
                tx = transcripts[parent]
                gene_ids = [g for g in tx["parents"] if g in genes]
                if not gene_ids and gene_hint in genes:
                    gene_ids = [gene_hint]
                for gene_id in gene_ids:
                    resolved.append((gene_id, parent, tx))
            elif parent in genes:
                resolved.append((parent, parent, {}))
        if not resolved and gene_hint in genes:
            # Some GFF3 exporters retain transcript Parent IDs on CDS while
            # storing the gene relationship only as gene_id.
            tx_id = parents[0] if parents else gene_hint
            resolved.append((gene_hint, tx_id, transcripts.get(tx_id, {})))
        for gene_id, tx_id, tx in resolved:
            key = (gene_id, tx_id)
            grouped.setdefault(key, []).append((start, end, chrom, strand, phase_value))
            transcript_meta[key] = tx

    model = GeneModel()
    candidates: dict[str, list[CDS]] = {}
    for (gene_id, transcript_id), segs in grouped.items():
        gene_meta = genes[gene_id]
        tx_meta = transcript_meta[(gene_id, transcript_id)]
        strand = tx_meta.get("strand", segs[0][3])
        chrom = tx_meta.get("chrom", segs[0][2])
        ordered = sorted(segs, key=lambda s: s[0], reverse=(strand == "-"))
        if strand not in ("+", "-") or any(s[2] != chrom or s[3] != strand for s in ordered):
            continue
        phases = [s[4] for s in ordered]
        if any(p is None for p in phases):
            continue
        initial_phase = int(phases[0])
        coding_length = ordered[0][1] - ordered[0][0] + 1 - initial_phase
        phase_ok = coding_length >= 0
        for segment, actual_phase in zip(ordered[1:], phases[1:], strict=True):
            expected_phase = (3 - coding_length % 3) % 3
            if actual_phase != expected_phase:
                phase_ok = False
                break
            coding_length += segment[1] - segment[0] + 1
        if not phase_ok:
            continue
        candidates.setdefault(gene_id, []).append(
            CDS(
                gene_id=gene_id,
                gene_name=gene_meta.get("name", tx_meta.get("name", gene_id)),
                chrom=chrom,
                strand=strand,
                segments=[(s[0], s[1]) for s in ordered],
                product=gene_meta.get("product", tx_meta.get("product", "")),
                transcript_id=transcript_id,
                phase=initial_phase,
            )
        )
    for models in candidates.values():
        # Stable, explicit tie-break: transcript IDs sort lexicographically.
        chosen = min(models, key=lambda cds: cds.transcript_id)
        model.add(chosen)
    return model
