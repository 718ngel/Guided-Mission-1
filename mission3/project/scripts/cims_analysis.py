"""
CIMS analysis: Shannon entropy per position from samtools pileup.
Produces bedgraph files for UCSC Genome Browser visualization.
"""
import re
import math
import pandas as pd
from collections import Counter
from pathlib import Path

WORK = Path(__file__).parent.parent / "data" / "work"
RESULTS = Path(__file__).parent.parent / "results"
RESULTS.mkdir(exist_ok=True)

GENES = {
    "let7g":  ("chr9",  106056039, 106056126, "+"),
    "let7f1": ("chr13",  48691305,  48691393, "-"),
    "let7d":  ("chr13",  48689488,  48689590, "-"),
}

# pileup special-character pattern: insertion (+Nxxx), deletion (-Nxxx), structural tags
_STRIP = re.compile(r'[<>$*#^]|\+\d+[ACGTNacgtn]+|-\d+[ACGTNacgtn]+')


def clean_basereads(s: str) -> str:
    return _STRIP.sub("", s)


def count_bases(s: str) -> Counter:
    return Counter(c.upper() for c in s if c.upper() in "ACGTN.,")


def shannon_entropy(counts: Counter) -> float:
    total = sum(counts.values())
    if total == 0:
        return 0.0
    return -sum(
        (n / total) * math.log2(n / total)
        for n in counts.values()
        if n > 0
    )


def process_gene(name: str, chrom: str, start: int, end: int, strand: str) -> pd.DataFrame:
    pileup_path = WORK / f"CLIP-{name}-gene.pileup"
    cols = ["chrom", "pos", "_ref", "count", "basereads", "quals"]
    df = pd.read_csv(pileup_path, sep="\t", names=cols)

    df["matches"] = df["basereads"].apply(clean_basereads)
    df["base_counts"] = df["matches"].apply(count_bases)
    df["entropy"] = df["base_counts"].apply(shannon_entropy)
    df["total_reads"] = df["base_counts"].apply(lambda c: sum(c.values()))
    return df


def write_bedgraph(df: pd.DataFrame, name: str, track_label: str) -> Path:
    out = RESULTS / f"CLIP-{name}-entropy.bedgraph"
    header = f'track type=bedGraph name="{track_label}" visibility=full color=0,100,200\n'
    rows = []
    for _, row in df.iterrows():
        s = int(row["pos"]) - 1  # 1-based pileup → 0-based bedgraph
        e = int(row["pos"])
        rows.append(f"{row['chrom']}\t{s}\t{e}\t{row['entropy']:.6f}")
    out.write_text(header + "\n".join(rows) + "\n")
    return out


def main():
    summary = []
    for name, (chrom, start, end, strand) in GENES.items():
        print(f"\n--- Processing {name} ({chrom}:{start}-{end} {strand}) ---")
        df = process_gene(name, chrom, start, end, strand)

        max_ent = df["entropy"].max()
        max_pos = df.loc[df["entropy"].idxmax(), "pos"]
        print(f"  Positions: {len(df)}")
        print(f"  Max entropy: {max_ent:.4f} at pos {int(max_pos)}")
        print(f"  Mean coverage: {df['total_reads'].mean():.1f} reads")

        label = f"CLIP-{name} Shannon entropy"
        path = write_bedgraph(df, name, label)
        print(f"  Bedgraph written: {path.name}")
        summary.append({"gene": name, "positions": len(df),
                         "max_entropy": round(max_ent, 4),
                         "max_entropy_pos": int(max_pos),
                         "mean_coverage": round(df["total_reads"].mean(), 1)})

    print("\n=== Summary ===")
    summary_df = pd.DataFrame(summary)
    print(summary_df.to_string(index=False))
    summary_df.to_csv(RESULTS / "entropy_summary.csv", index=False)
    print(f"\nSummary saved to results/entropy_summary.csv")


if __name__ == "__main__":
    main()
