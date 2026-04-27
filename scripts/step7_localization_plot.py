#!/usr/bin/env python3
"""
Step 7: Scatter plot colored by protein subcellular localization.
Reproduces Figure 5B / S6A from the paper.
"""

import os
import ssl
import io
import urllib.request
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

ROOT = os.path.join(os.path.dirname(__file__), "..")
DATA_DIR = os.path.join(ROOT, "data")
OUTPUT_DIR = os.path.join(ROOT, "output")
os.makedirs(OUTPUT_DIR, exist_ok=True)

LOCALIZATION_URL = "https://hyeshik.qbio.io/binfo/mouselocalization-20210507.txt"

LOC_COLORS = {
    "cytoplasm":        "#1f77b4",   # blue
    "nucleus":          "#ff7f0e",   # orange
    "integral membrane":"#2ca02c",   # green
}

# --- Load counts ---
cnts = pd.read_csv(
    os.path.join(DATA_DIR, "read-counts.txt"),
    sep="\t", comment="#", index_col=0
)

cnts["clip_enrichment"] = cnts["CLIP-35L33G.bam"] / cnts["RNA-control.bam"]
cnts["rden_change"] = (
    cnts["RPF-siLin28a.bam"] / cnts["RNA-siLin28a.bam"]
) / (
    cnts["RPF-siLuc.bam"] / cnts["RNA-siLuc.bam"]
)

# Strip version suffix (e.g. ENSMUSG00000102693.2 -> ENSMUSG00000102693)
cnts["gene_id"] = cnts.index.str.split(".").str[0]

# --- Download localization data ---
print("Downloading localization data...")
ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE
with urllib.request.urlopen(LOCALIZATION_URL, context=ctx) as r:
    mouselocal = pd.read_csv(io.BytesIO(r.read()), sep="\t")
print(f"Localization entries: {len(mouselocal)}")

# --- Merge ---
merged = cnts.merge(mouselocal[["gene_id", "type"]], on="gene_id", how="left")

# Keep only rows with finite log2 values
x = np.log2(merged["clip_enrichment"])
y = np.log2(merged["rden_change"])
finite = np.isfinite(x) & np.isfinite(y)
merged = merged[finite].copy()
x = x[finite]
y = y[finite]
print(f"Genes plotted: {len(merged)}  (with localization: {merged['type'].notna().sum()})")

# --- Plot ---
fig, ax = plt.subplots(figsize=(8, 8))

# Unknown / no localization (gray background)
mask_none = merged["type"].isna()
ax.scatter(x[mask_none], y[mask_none],
           c="lightgray", alpha=0.25, s=5, linewidths=0, label="unknown")

# Each localization category
for loc, color in LOC_COLORS.items():
    mask = merged["type"] == loc
    ax.scatter(x[mask], y[mask],
               c=color, alpha=0.6, s=10, linewidths=0, label=loc)

ax.axhline(0, color="black", linewidth=0.8, linestyle="--", alpha=0.5)
ax.axvline(0, color="black", linewidth=0.8, linestyle="--", alpha=0.5)
ax.set_xlabel("log2(CLIP Enrichment)", fontsize=12)
ax.set_ylabel("log2(Ribosome Density Change)", fontsize=12)
ax.set_title("CLIP Enrichment vs Ribosome Density Change\n(colored by protein localization)", fontsize=12)

legend_handles = [
    mpatches.Patch(color=c, label=loc) for loc, c in LOC_COLORS.items()
] + [mpatches.Patch(color="lightgray", label="unknown")]
ax.legend(handles=legend_handles, loc="upper left", fontsize=10, framealpha=0.8)

out = os.path.join(OUTPUT_DIR, "scatter_localization.png")
plt.savefig(out, dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved: {out}")
