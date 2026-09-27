"""Command-line interface: run the classifier on a FASTA + variants CSV."""

from __future__ import annotations

import argparse
from pathlib import Path

from .classifier import classify_variants, load_fasta, load_variants_csv, plot_consequences


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        prog="variant-classify",
        description="Classify candidate coding-sequence variants as synonymous, "
        "missense, nonsense, or frameshift.",
    )
    parser.add_argument("--fasta", required=True, help="Reference coding sequence (FASTA)")
    parser.add_argument("--variants", required=True, help="Variants CSV: variant_id,pos,ref,alt")
    parser.add_argument("--out", default="results", help="Output directory (default: results/)")
    args = parser.parse_args(argv)

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    reference_seq = load_fasta(args.fasta)
    variants = load_variants_csv(args.variants)

    df = classify_variants(reference_seq, variants)
    csv_path = out_dir / "variant_report.csv"
    chart_path = out_dir / "variant_consequence_chart.png"

    df.to_csv(csv_path, index=False)
    plot_consequences(df, chart_path)

    print(df.to_string(index=False))
    print(f"\nSaved: {csv_path}")
    print(f"Saved: {chart_path}")


if __name__ == "__main__":
    main()
