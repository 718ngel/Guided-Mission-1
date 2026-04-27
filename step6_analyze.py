#!/usr/bin/env python3
"""
Step 6: Load featureCounts results and create basic scatter plot.
Computes CLIP enrichment and ribosome density change, then plots log2 values.
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "output")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Load counts
counts_path = os.path.join(DATA_DIR, "read-counts.txt")
cnts = pd.read_csv(counts_path, sep="\t", comment="#", index_col=0)
print(f"Loaded {len(cnts)} genes, columns: {list(cnts.columns)}")

# Compute metrics
cnts["clip_enrichment"] = cnts["CLIP-35L33G.bam"] / cnts["RNA-control.bam"]
cnts["rden_change"] = (
    cnts["RPF-siLin28a.bam"] / cnts["RNA-siLin28a.bam"]
) / (
    cnts["RPF-siLuc.bam"] / cnts["RNA-siLuc.bam"]
)

# Filter to finite values only
mask = np.isfinite(np.log2(cnts["clip_enrichment"])) & np.isfinite(np.log2(cnts["rden_change"]))
cnts_valid = cnts[mask]
print(f"Genes with finite log2 values: {len(cnts_valid)}")

# Basic scatter plot
fig, ax = plt.subplots(figsize=(8, 8))
ax.scatter(
    np.log2(cnts_valid["clip_enrichment"]),
    np.log2(cnts_valid["rden_change"]),
    alpha=0.4,
    s=5,
    color="steelblue",
)
ax.set_xlabel("log2(CLIP Enrichment)", fontsize=12)
ax.set_ylabel("log2(Ribosome Density Change)", fontsize=12)
ax.set_title("CLIP Enrichment vs Ribosome Density Change", fontsize=13)
ax.axhline(y=0, color="gray", linestyle="--", linewidth=0.8, alpha=0.7)
ax.axvline(x=0, color="gray", linestyle="--", linewidth=0.8, alpha=0.7)

out = os.path.join(OUTPUT_DIR, "scatter_basic.png")
plt.savefig(out, dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved: {out}")
