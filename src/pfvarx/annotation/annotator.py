"""Annotation orchestration: consequence prediction plus resistance matching."""

from __future__ import annotations

from dataclasses import dataclass

import pysam

from pfvarx.annotation.catalogue import ResistanceCatalogue, load_catalogue
from pfvarx.annotation.consequence import annotate_variant
from pfvarx.annotation.gff import GeneModel, parse_gff3
from pfvarx.calling.policy import SiteCall
from pfvarx.config.schema import PfVarXConfig
from pfvarx.errors import ToolExecutionError
from pfvarx.preprocess.tools import TOOLS, run_command


@dataclass
class AnnotatedSite:
    call: SiteCall
    consequences: list[dict]
    resistance: list[dict]
    engine: str

    def as_dict(self) -> dict:
        base = self.call.as_dict()
        base["annotation_engine"] = self.engine
        base["consequence"] = ";".join(c["consequence"] for c in self.consequences) or "."
        base["gene_id"] = ";".join(str(c["gene_id"]) for c in self.consequences) or "."
        base["gene_name"] = ";".join(str(c["gene_name"]) for c in self.consequences) or "."
        base["amino_acid_change"] = ";".join(str(c["amino_acid_change"]) for c in self.consequences) or "."
        base["resistance_drugs"] = ";".join(r["drug"] for r in self.resistance) or "."
        base["resistance_markers"] = ";".join(r["aa_change"] for r in self.resistance) or "."
        base["resistance_evidence"] = ";".join(r["evidence_level"] for r in self.resistance) or "."
        return base


def _snpeff_annotate(cfg: PfVarXConfig, vcf_path: str, out_vcf: str) -> None:
    spec = TOOLS["snpeff"]
    exe = spec.require()
    if not cfg.reference.snpeff_db:
        raise ToolExecutionError(
            "SnpEff annotation requested but no database is configured.",
            cause="reference.snpeff_db is empty.",
            action="Set reference.snpeff_db or use engine='internal'.",
        )
    run_command([exe, "-v", cfg.reference.snpeff_db, vcf_path], log_path=out_vcf + ".log")


def annotate_calls(
    calls: list[SiteCall],
    cfg: PfVarXConfig,
    *,
    gff_path: str | None = None,
    fasta_path: str | None = None,
    catalogue: ResistanceCatalogue | None = None,
) -> list[AnnotatedSite]:
    """Annotate calls using the internal engine (SnpEff is optional)."""
    engine = cfg.annotation.engine
    if engine == "snpeff":
        # The SnpEff round-trip is not wired into the per-call annotation path;
        # silently falling back to the internal engine would mislabel the
        # provenance.  Fail closed with an actionable message instead.
        TOOLS["snpeff"].require()
        raise ToolExecutionError(
            "annotation.engine='snpeff' is not implemented for per-call annotation",
            cause="SnpEff VCF round-trip is not wired into annotate_calls().",
            action="Set annotation.engine='internal' (implemented, not scientifically validated) or "
            "run SnpEff externally before import.",
        )

    cat = catalogue or load_catalogue()
    model = parse_gff3(gff_path) if gff_path else GeneModel()
    fasta = pysam.FastaFile(fasta_path) if fasta_path else None

    annotated: list[AnnotatedSite] = []
    try:
        for call in calls:
            ev = call.evidence
            if model.cds_by_gene and fasta is not None:
                cons = [
                    c.as_dict()
                    for c in annotate_variant(model, fasta, ev.chrom, ev.pos, ev.ref, ev.alt)
                ]
            else:
                cons = [
                    {
                        "gene_id": ".",
                        "gene_name": ".",
                        "transcript_biotype": "unannotated",
                        "consequence": "unannotated",
                        "amino_acid_change": ".",
                        "codon_change": ".",
                        "protein_position": ".",
                        "severity_rank": 99,
                        "provenance": "NONE",
                    }
                ]
            resistance: list[dict] = []
            for c in cons:
                if c["gene_id"] != ".":
                    for m in cat.match(str(c["gene_id"]), str(c["amino_acid_change"])):
                        resistance.append(m.as_dict())
            annotated.append(
                AnnotatedSite(call=call, consequences=cons, resistance=resistance, engine="internal")
            )
    finally:
        if fasta is not None:
            fasta.close()
    return annotated
