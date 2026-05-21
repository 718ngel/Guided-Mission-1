#!/usr/bin/env python3
"""
Analysis 1 – Translation Efficiency (TE) change upon Lin28a knockdown
        + Correlation with protein subcellular localization

Key correction vs. original step6/step7:
  Raw-count ratios are replaced by CPM-normalised values so that
  different library depths between RNA and RPF samples are accounted for.

Outputs (saved to ../output/):
  analysis1_TE_scatter.png       – ΔTE map: ΔRNA (x) vs ΔRPF (y), color = ΔTE
  analysis1_localization.png     – Violin plot of ΔTE per localisation group
  analysis1_results.csv          – Per-gene TE table (filtered genes only)
"""

import os, io, ssl, sys
import urllib.request
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy import stats

# ── Paths ─────────────────────────────────────────────────────────────────────
ROOT       = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR   = os.path.join(ROOT, "data")
OUTPUT_DIR = os.path.join(ROOT, "output")
os.makedirs(OUTPUT_DIR, exist_ok=True)

COUNTS_FILE   = os.path.join(DATA_DIR, "read-counts.txt")
LOCALIZE_URL  = "https://hyeshik.qbio.io/binfo/mouselocalization-20210507.txt"

# Sample column names (as they appear in featureCounts output)
RNA_CTRL  = "RNA-siLuc.bam"
RNA_KD    = "RNA-siLin28a.bam"
RPF_CTRL  = "RPF-siLuc.bam"
RPF_KD    = "RPF-siLin28a.bam"

PSEUDO    = 0.5   # CPM pseudocount – stabilises log2 for very low-expressed genes
MIN_RNA   = 1.0   # minimum CPM in both RNA conditions
MIN_RPF   = 5     # minimum raw reads in both RPF conditions

# ── 1. Load feature counts ─────────────────────────────────────────────────────
print("Loading featureCounts matrix …")
cnts = pd.read_csv(COUNTS_FILE, sep="\t", comment="#", index_col=0)
print(f"  {len(cnts):,} genes × {len(cnts.columns)} columns loaded")

# ── 2. CPM normalisation ───────────────────────────────────────────────────────
def to_cpm(series: pd.Series) -> pd.Series:
    return series / series.sum() * 1e6

cnts["rna_ctrl_cpm"] = to_cpm(cnts[RNA_CTRL])
cnts["rna_kd_cpm"]   = to_cpm(cnts[RNA_KD])
cnts["rpf_ctrl_cpm"] = to_cpm(cnts[RPF_CTRL])
cnts["rpf_kd_cpm"]   = to_cpm(cnts[RPF_KD])

# ── 3. Low-count filtering ─────────────────────────────────────────────────────
# Genes with very low expression produce unstable TE estimates.
mask = (
    (cnts["rna_ctrl_cpm"] >= MIN_RNA) &
    (cnts["rna_kd_cpm"]   >= MIN_RNA) &
    (cnts[RPF_CTRL]       >= MIN_RPF) &
    (cnts[RPF_KD]         >= MIN_RPF)
)
df = cnts[mask].copy()
print(f"  After filtering (RNA ≥ {MIN_RNA} CPM, RPF ≥ {MIN_RPF} raw reads): {len(df):,} genes")

# ── 4. Translation Efficiency ──────────────────────────────────────────────────
# TE = RPF_CPM / RNA_CPM  (pseudocount added to denominator only)
df["TE_ctrl"] = df["rpf_ctrl_cpm"] / (df["rna_ctrl_cpm"] + PSEUDO)
df["TE_kd"]   = df["rpf_kd_cpm"]   / (df["rna_kd_cpm"]   + PSEUDO)

df["delta_TE"]  = np.log2(df["TE_kd"])  - np.log2(df["TE_ctrl"])
df["delta_RNA"] = np.log2(df["rna_kd_cpm"] + PSEUDO) - np.log2(df["rna_ctrl_cpm"] + PSEUDO)
df["delta_RPF"] = np.log2(df["rpf_kd_cpm"] + PSEUDO) - np.log2(df["rpf_ctrl_cpm"] + PSEUDO)

print(f"\nΔTE  — mean: {df['delta_TE'].mean():.3f}  median: {df['delta_TE'].median():.3f}")
print(f"ΔRNA — mean: {df['delta_RNA'].mean():.3f}  median: {df['delta_RNA'].median():.3f}")
print(f"ΔRPF — mean: {df['delta_RPF'].mean():.3f}  median: {df['delta_RPF'].median():.3f}")

# ── 5. Scatter: ΔRNA (x) vs ΔRPF (y) coloured by ΔTE ────────────────────────
fig, ax = plt.subplots(figsize=(9, 8))
norm = mcolors.TwoSlopeNorm(vcenter=0, vmin=-2, vmax=2)
sc = ax.scatter(
    df["delta_RNA"], df["delta_RPF"],
    c=df["delta_TE"], cmap="RdBu_r", norm=norm,
    alpha=0.35, s=6, linewidths=0,
)
cbar = plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label("ΔTE (log2)", fontsize=11)

# Diagonal = ΔTE 0 line
lim = 6
ax.plot([-lim, lim], [-lim, lim], "k--", linewidth=0.8, alpha=0.5, label="ΔTE = 0")
ax.axhline(0, color="gray", linewidth=0.5, linestyle=":")
ax.axvline(0, color="gray", linewidth=0.5, linestyle=":")
ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)
ax.set_xlabel("Δlog₂(RNA CPM)  [siLin28a / siLuc]", fontsize=12)
ax.set_ylabel("Δlog₂(RPF CPM)  [siLin28a / siLuc]", fontsize=12)
ax.set_title("Translation Efficiency change upon Lin28a knockdown\n"
             f"(n = {len(df):,} genes, CPM-normalised, {MIN_RNA} CPM / {MIN_RPF} raw RPF filter)",
             fontsize=11)
ax.legend(fontsize=9)
plt.tight_layout()
out1 = os.path.join(OUTPUT_DIR, "analysis1_TE_scatter.png")
plt.savefig(out1, dpi=150)
plt.close()
print(f"\n[saved] {out1}")

# ── 6. Cellular localisation analysis ─────────────────────────────────────────
print("\nDownloading localisation data …")
ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE
with urllib.request.urlopen(LOCALIZE_URL, context=ctx) as resp:
    locdf = pd.read_csv(io.BytesIO(resp.read()), sep="\t")
print(f"  {len(locdf):,} entries loaded")

df["gene_id"] = df.index.str.split(".").str[0]
merged = df.merge(locdf[["gene_id", "type"]], on="gene_id", how="left")

LOC_ORDER = ["cytoplasm", "nucleus", "integral membrane"]
groups = {loc: merged[merged["type"] == loc]["delta_TE"].dropna() for loc in LOC_ORDER}
groups["unknown"] = merged[merged["type"].isna()]["delta_TE"].dropna()

for loc, g in groups.items():
    print(f"  {loc:20s}: n={len(g):4d}  mean ΔTE={g.mean():.3f}  median={g.median():.3f}")

# One-way ANOVA across the three annotated groups
f_stat, p_anova = stats.f_oneway(*[groups[k] for k in LOC_ORDER])
print(f"\nANOVA (cytoplasm vs nucleus vs membrane): F={f_stat:.3f}, p={p_anova:.4e}")

# Pairwise t-tests with Bonferroni correction
pairs = [("cytoplasm", "nucleus"), ("cytoplasm", "integral membrane"), ("nucleus", "integral membrane")]
print("\nPairwise t-tests (Bonferroni α = 0.0167):")
for a, b in pairs:
    t, p = stats.ttest_ind(groups[a], groups[b])
    print(f"  {a} vs {b}: t={t:.3f}, p={p:.4e}  {'*' if p < 0.0167 else 'ns'}")

# ── 7. Violin plot ─────────────────────────────────────────────────────────────
LOC_COLORS = {
    "cytoplasm":        "#1f77b4",
    "nucleus":          "#ff7f0e",
    "integral membrane":"#2ca02c",
    "unknown":          "#d3d3d3",
}
all_labels  = LOC_ORDER + ["unknown"]
data_arrays = [groups[k].values for k in all_labels]

fig, ax = plt.subplots(figsize=(10, 6))
parts = ax.violinplot(data_arrays, positions=range(len(all_labels)),
                      showmedians=True, showextrema=False)
for i, (pc, lbl) in enumerate(zip(parts["bodies"], all_labels)):
    pc.set_facecolor(LOC_COLORS[lbl])
    pc.set_alpha(0.7)
parts["cmedians"].set_color("black")
parts["cmedians"].set_linewidth(1.5)

for i, lbl in enumerate(all_labels):
    n = len(groups[lbl])
    ax.text(i, ax.get_ylim()[0] if ax.get_ylim()[0] else -4,
            f"n={n}", ha="center", va="bottom", fontsize=8, color="gray")

ax.axhline(0, color="black", linewidth=0.8, linestyle="--", alpha=0.5)
ax.set_xticks(range(len(all_labels)))
ax.set_xticklabels([l.replace(" ", "\n") for l in all_labels], fontsize=11)
ax.set_ylabel("ΔTE  (log₂ fold-change)", fontsize=12)
ax.set_title(
    f"ΔTE distribution by protein localisation\n"
    f"ANOVA: F={f_stat:.2f}, p={p_anova:.2e}",
    fontsize=11,
)
plt.tight_layout()
out2 = os.path.join(OUTPUT_DIR, "analysis1_localization.png")
plt.savefig(out2, dpi=150)
plt.close()
print(f"[saved] {out2}")

# ── 8. Save results table ──────────────────────────────────────────────────────
merged_out = merged[["gene_id", "type", "rna_ctrl_cpm", "rna_kd_cpm",
                      "rpf_ctrl_cpm", "rpf_kd_cpm",
                      "TE_ctrl", "TE_kd", "delta_TE", "delta_RNA", "delta_RPF"]]
csv_out = os.path.join(OUTPUT_DIR, "analysis1_results.csv")
merged_out.to_csv(csv_out, index=True)
print(f"[saved] {csv_out}")
print("\nDone – Analysis 1 complete.")
