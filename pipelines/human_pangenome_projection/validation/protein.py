"""Protein-level quality checks for mapped coding transcripts."""

from dataclasses import dataclass, field
from difflib import SequenceMatcher
from typing import Dict, List, Optional

from src.fasta_handler import FastaHandler
from src.models import Gene, Transcript


@dataclass
class TranscriptProteinQC:
    """Protein-level QC result for one transcript."""
    transcript_id: str
    reference_transcript_id: Optional[str]
    status: str
    identity: float = 0.0
    coverage: float = 0.0
    internal_stop: bool = False
    reference_protein_length: int = 0
    target_protein_length: int = 0
    note: Optional[str] = None

    def to_dict(self) -> Dict:
        return {
            "transcript_id": self.transcript_id,
            "reference_transcript_id": self.reference_transcript_id,
            "status": self.status,
            "identity": self.identity,
            "coverage": self.coverage,
            "internal_stop": self.internal_stop,
            "reference_protein_length": self.reference_protein_length,
            "target_protein_length": self.target_protein_length,
            "note": self.note,
        }


@dataclass
class GeneProteinQC:
    """Protein-level QC summary for one gene."""
    gene_id: str
    mapped_from: Optional[str]
    status: str
    transcripts_checked: int = 0
    transcripts_ok: int = 0
    transcripts_with_internal_stop: int = 0
    mean_identity: float = 0.0
    mean_coverage: float = 0.0
    transcript_results: List[TranscriptProteinQC] = field(default_factory=list)

    def to_dict(self) -> Dict:
        return {
            "gene_id": self.gene_id,
            "mapped_from": self.mapped_from,
            "status": self.status,
            "transcripts_checked": self.transcripts_checked,
            "transcripts_ok": self.transcripts_ok,
            "transcripts_with_internal_stop": self.transcripts_with_internal_stop,
            "mean_identity": self.mean_identity,
            "mean_coverage": self.mean_coverage,
            "transcript_results": [t.to_dict() for t in self.transcript_results],
        }


def _translate_transcript(transcript: Transcript, fasta: FastaHandler) -> str:
    """Translate transcript CDS to amino-acid sequence (with terminal stop if present)."""
    if not transcript.cds_list:
        return ""
    sorted_cds = transcript.sorted_cds
    cds_sequences = []
    for cds in sorted_cds:
        # Defensive guard: malformed intervals should not crash whole-pipeline QC.
        if cds.start > cds.end:
            return ""
        try:
            cds_sequences.append(
                fasta.fetch(cds.seq_region, cds.start, cds.end, strand=transcript.strand)
            )
        except Exception:
            return ""
    phase = sorted_cds[0].phase or 0
    return fasta.translate_cds(cds_sequences, phase=phase)


def _protein_identity(reference_protein: str, target_protein: str) -> float:
    """Reference-normalized sequence identity using global matching blocks."""
    if not reference_protein:
        return 0.0
    matcher = SequenceMatcher(None, reference_protein, target_protein, autojunk=False)
    matches = sum(block.size for block in matcher.get_matching_blocks())
    return matches / max(1, len(reference_protein))


def evaluate_gene_protein_qc(
    original_gene: Gene,
    mapped_gene: Gene,
    ref_fasta: FastaHandler,
    target_fasta: FastaHandler,
    min_identity: float = 0.90,
    min_coverage: float = 0.90,
) -> GeneProteinQC:
    """Evaluate protein-level agreement for one mapped gene."""
    orig_by_tx_id = {
        tx.feature_id: tx for tx in original_gene.transcripts
    }
    tx_results: List[TranscriptProteinQC] = []

    for mapped_tx in mapped_gene.transcripts:
        if not mapped_tx.is_coding:
            continue

        orig_tx = orig_by_tx_id.get(mapped_tx.mapped_from or "")
        if orig_tx is None or not orig_tx.is_coding:
            tx_results.append(
                TranscriptProteinQC(
                    transcript_id=mapped_tx.feature_id,
                    reference_transcript_id=mapped_tx.mapped_from,
                    status="unmatched_reference_transcript",
                    note="No matching coding transcript in reference",
                )
            )
            continue

        ref_protein = _translate_transcript(orig_tx, ref_fasta)
        tgt_protein = _translate_transcript(mapped_tx, target_fasta)
        if not ref_protein or not tgt_protein:
            tx_results.append(
                TranscriptProteinQC(
                    transcript_id=mapped_tx.feature_id,
                    reference_transcript_id=orig_tx.feature_id,
                    status="untranslated",
                    reference_protein_length=len(ref_protein.rstrip("*")),
                    target_protein_length=len(tgt_protein.rstrip("*")),
                )
            )
            continue

        internal_stop = "*" in tgt_protein[:-1]
        ref_clean = ref_protein.rstrip("*")
        tgt_clean = tgt_protein.rstrip("*")
        identity = _protein_identity(ref_clean, tgt_clean)
        coverage = len(tgt_clean) / max(1, len(ref_clean))

        if internal_stop:
            status = "internal_stop"
        elif identity >= min_identity and coverage >= min_coverage:
            status = "ok"
        else:
            status = "diverged"

        tx_results.append(
            TranscriptProteinQC(
                transcript_id=mapped_tx.feature_id,
                reference_transcript_id=orig_tx.feature_id,
                status=status,
                identity=identity,
                coverage=coverage,
                internal_stop=internal_stop,
                reference_protein_length=len(ref_clean),
                target_protein_length=len(tgt_clean),
            )
        )

    checked = len(tx_results)
    if checked == 0:
        return GeneProteinQC(
            gene_id=mapped_gene.feature_id,
            mapped_from=mapped_gene.mapped_from,
            status="no_coding_transcripts",
            transcript_results=[],
        )

    identities = [t.identity for t in tx_results if t.reference_protein_length > 0]
    coverages = [t.coverage for t in tx_results if t.reference_protein_length > 0]
    transcripts_ok = sum(1 for t in tx_results if t.status == "ok")
    transcripts_internal_stop = sum(1 for t in tx_results if t.internal_stop)

    if transcripts_internal_stop > 0:
        gene_status = "internal_stop"
    elif transcripts_ok == checked:
        gene_status = "ok"
    elif any(t.status == "diverged" for t in tx_results):
        gene_status = "diverged"
    else:
        gene_status = "partial"

    return GeneProteinQC(
        gene_id=mapped_gene.feature_id,
        mapped_from=mapped_gene.mapped_from,
        status=gene_status,
        transcripts_checked=checked,
        transcripts_ok=transcripts_ok,
        transcripts_with_internal_stop=transcripts_internal_stop,
        mean_identity=(sum(identities) / len(identities)) if identities else 0.0,
        mean_coverage=(sum(coverages) / len(coverages)) if coverages else 0.0,
        transcript_results=tx_results,
    )


def run_protein_qc(
    original_genes: Dict[str, Gene],
    mapped_genes: Dict[str, Gene],
    ref_fasta: FastaHandler,
    target_fasta: FastaHandler,
    min_identity: float = 0.90,
    min_coverage: float = 0.90,
) -> Dict[str, GeneProteinQC]:
    """Run protein-level QC for all mapped genes with reference provenance."""
    results: Dict[str, GeneProteinQC] = {}
    for mapped_id, mapped_gene in mapped_genes.items():
        ref_id = mapped_gene.mapped_from
        if not ref_id:
            continue
        original_gene = original_genes.get(ref_id)
        if original_gene is None:
            continue
        try:
            results[mapped_id] = evaluate_gene_protein_qc(
                original_gene,
                mapped_gene,
                ref_fasta,
                target_fasta,
                min_identity=min_identity,
                min_coverage=min_coverage,
            )
        except Exception as exc:
            results[mapped_id] = GeneProteinQC(
                gene_id=mapped_gene.feature_id,
                mapped_from=mapped_gene.mapped_from,
                status="qc_error",
                transcript_results=[
                    TranscriptProteinQC(
                        transcript_id=t.feature_id,
                        reference_transcript_id=t.mapped_from,
                        status="qc_error",
                        note=str(exc),
                    )
                    for t in mapped_gene.transcripts
                ],
            )
    return results
