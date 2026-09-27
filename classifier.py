"""
variant_classifier.classifier
------------------------------
Core logic for predicting the consequence of a coding-sequence variant:
synonymous, missense, nonsense, or frameshift.

This is a simplified, from-scratch version of what production annotation
tools such as Ensembl VEP or ANNOVAR do. It's intended as a learning /
portfolio project, not a clinical-grade annotator.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from Bio import BiopythonWarning
from Bio.Seq import Seq

warnings.simplefilter("ignore", BiopythonWarning)


@dataclass
class Variant:
    """A single candidate variant, given as a 0-based substitution/indel
    against the reference coding sequence.

    alt="" represents a deletion of len(ref) bases starting at pos.
    len(alt) > len(ref) with ref="" (or ref shorter than alt) represents
    an insertion.
    """

    variant_id: str
    pos: int
    ref: str
    alt: str


def load_fasta(path: str | Path) -> str:
    """Load the first sequence record from a FASTA file as an uppercase string."""
    from Bio import SeqIO

    record = next(SeqIO.parse(str(path), "fasta"))
    return str(record.seq).upper()


def load_variants_csv(path: str | Path) -> list[Variant]:
    """Load variants from a CSV with columns: variant_id,pos,ref,alt."""
    df = pd.read_csv(path, keep_default_na=False, dtype={"ref": str, "alt": str})
    return [
        Variant(row.variant_id, int(row.pos), row.ref, row.alt)
        for row in df.itertuples(index=False)
    ]


def _trim_to_codon_boundary(seq: str) -> str:
    return seq[: len(seq) - (len(seq) % 3)]


def apply_variant(seq: str, variant: Variant) -> str:
    """Apply a single substitution, deletion, or insertion to `seq`."""
    pos, ref, alt = variant.pos, variant.ref, variant.alt
    observed = seq[pos : pos + len(ref)]
    if observed != ref:
        raise ValueError(
            f"{variant.variant_id}: reference mismatch at pos {pos} "
            f"(expected {ref!r}, found {observed!r})"
        )
    return seq[:pos] + alt + seq[pos + len(ref) :]


def classify_variant(reference_seq: str, variant: Variant) -> dict:
    """Classify a single variant's predicted effect on the translated protein.

    Returns a dict with predicted_consequence, clinically_flagged, and a
    preview of the mutant protein sequence.
    """
    ref_trimmed = _trim_to_codon_boundary(reference_seq)
    protein_ref = str(Seq(ref_trimmed).translate(to_stop=False))

    mutated_seq = apply_variant(reference_seq, variant)
    mutated_trimmed = _trim_to_codon_boundary(mutated_seq)
    protein_mut = str(Seq(mutated_trimmed).translate(to_stop=False))

    is_indel = len(variant.ref) != len(variant.alt)
    length_delta = len(variant.alt) - len(variant.ref)

    if is_indel and length_delta % 3 != 0:
        consequence = "frameshift"
    else:
        codon_idx = variant.pos // 3
        aa_ref = protein_ref[codon_idx] if codon_idx < len(protein_ref) else "?"
        aa_mut = protein_mut[codon_idx] if codon_idx < len(protein_mut) else "?"
        if aa_ref == aa_mut:
            consequence = "synonymous"
        elif aa_mut == "*":
            consequence = "nonsense"
        elif is_indel:
            consequence = "in-frame indel"
        else:
            consequence = "missense"

    clinically_flagged = consequence in ("nonsense", "frameshift")

    return {
        "variant_id": variant.variant_id,
        "pos": variant.pos,
        "ref": variant.ref,
        "alt": variant.alt,
        "predicted_consequence": consequence,
        "clinically_flagged": clinically_flagged,
        "mutant_protein_preview": protein_mut[:60],
    }


def classify_variants(reference_seq: str, variants: list[Variant]) -> pd.DataFrame:
    """Classify a batch of variants and return a results DataFrame."""
    results = [classify_variant(reference_seq, v) for v in variants]
    return pd.DataFrame(results)


def plot_consequences(df: pd.DataFrame, output_path: str | Path) -> None:
    """Bar chart of predicted consequence counts, saved to `output_path`."""
    import matplotlib.pyplot as plt

    colors = {
        "synonymous": "#8FBF8F",
        "missense": "#F2C14E",
        "nonsense": "#E15554",
        "frameshift": "#7A2E2E",
        "in-frame indel": "#6C8EBF",
    }
    counts = df["predicted_consequence"].value_counts()

    fig, ax = plt.subplots(figsize=(6, 4))
    bars = ax.bar(counts.index, counts.values, color=[colors.get(c, "#999") for c in counts.index])
    ax.set_title(f"Predicted Variant Consequences (n={len(df)} candidate variants)")
    ax.set_ylabel("Count")
    for bar in bars:
        h = bar.get_height()
        ax.annotate(str(int(h)), (bar.get_x() + bar.get_width() / 2, h), ha="center", va="bottom")
    plt.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
