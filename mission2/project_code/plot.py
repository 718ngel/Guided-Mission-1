#!/usr/bin/env python3
"""
Mission 2 – Plot Figure S5A
RPF 5'-end distribution around start codons.

Run after run_pipeline.sh:
    python3 project_code/plot.py
"""

import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

ROOT      = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORK_DIR  = os.path.join(ROOT, "work")
RESULTS   = os.path.join(ROOT, "results")
os.makedirs(RESULTS, exist_ok=True)

INTER_TXT = os.path.join(WORK_DIR, "fivepcounts-filtered-RPF-siLuc.txt")

COLS = [
    "read_chr", "read_start", "read_end", "count",
    "exon_chr", "exon_start", "exon_end",
    "transcript_id", "startcodon_pos", "strand",
]
WINDOW = (-50, 150)
CHUNK  = 500_000

print("Loading intersection file …")
profile = pd.Series(dtype=float)
total_rows = 0

for chunk in pd.read_csv(
    INTER_TXT, sep="\t", header=None, names=COLS, chunksize=CHUNK,
):
    total_rows += len(chunk)
    chunk["rel_pos"] = chunk["read_start"] - chunk["startcodon_pos"]
    chunk = chunk[
        (chunk["rel_pos"] >= WINDOW[0]) & (chunk["rel_pos"] <= WINDOW[1])
    ]
    profile = profile.add(chunk.groupby("rel_pos")["count"].sum(), fill_value=0)

profile = profile.sort_index()

print(f"Total rows processed : {total_rows:,}")
print(f"Positions in window  : {len(profile)}")
print(f"Peak position        : {int(profile.idxmax())} nt  (count = {int(profile.max()):,})")

# ── Save CSV summary ──────────────────────────────────────────────────────────
csv_path = os.path.join(RESULTS, "rpf_profile_around_startcodon.csv")
pd.DataFrame({
    "position_nt":    profile.index.astype(int),
    "rpf_5end_count": profile.values.astype(int),
}).to_csv(csv_path, index=False)
print(f"CSV saved  → {csv_path}")

# ── Save stats summary ────────────────────────────────────────────────────────
window_data = profile[(profile.index >= -15) & (profile.index <= 100)]
stats = {
    "total_5end_counts_in_window":   int(profile.sum()),
    "peak_position_nt":              int(profile.idxmax()),
    "peak_count":                    int(profile.max()),
    "mean_count_per_position":       round(float(profile.mean()), 1),
    "positions_with_nonzero_count":  int((profile > 0).sum()),
}
stats_path = os.path.join(RESULTS, "summary_stats.txt")
with open(stats_path, "w") as f:
    f.write("RPF 5'-end profile around start codons – summary stats\n")
    f.write("=" * 55 + "\n")
    for k, v in stats.items():
        f.write(f"  {k:<40} {v}\n")
print(f"Stats saved → {stats_path}")

# ── Plot Figure S5A ───────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(14, 5))

bar_colors = ["#d62728" if int(i) == 0 else "#1f77b4" for i in profile.index]
ax.bar(profile.index, profile.values, width=1, color=bar_colors, alpha=0.85)

ax.axvline(0,   color="red",    linestyle="--", linewidth=1.5,
           label="Start codon (A of AUG, pos 0)")
ax.axvline(-15, color="orange", linestyle=":",  linewidth=1.2,
           label="−15 nt (expected ribosome 5'-end peak)")

ax.set_xlim(WINDOW)
ax.set_xlabel("Position relative to start codon A (nt)", fontsize=12)
ax.set_ylabel("RPF 5'-end count", fontsize=12)
ax.set_title(
    "Figure S5A – RPF 5'-end distribution around start codons\n"
    "(+ strand only · RPF-siLuc · GENCODE M27 · TSL-1 transcripts)",
    fontsize=13,
)
ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(3))
ax.legend(fontsize=10)
ax.grid(axis="y", alpha=0.3, linestyle="--")
plt.tight_layout()

png_path = os.path.join(RESULTS, "figure_S5A_rpf_startcodon.png")
pdf_path = os.path.join(RESULTS, "figure_S5A_rpf_startcodon.pdf")
plt.savefig(png_path, dpi=200)
plt.savefig(pdf_path)
print(f"PNG saved  → {png_path}")
print(f"PDF saved  → {pdf_path}")

print("\nDone. All results in:  results/")
