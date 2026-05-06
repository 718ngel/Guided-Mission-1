#!/usr/bin/env python3
"""
Mission 2 – RPF 5'-end distribution around start codons
Adapted from Colab notebook for local VS Code use.

Usage:
    python scripts/analysis.py

Prerequisites (install once):
    conda install -y -c bioconda bedtools bioawk samtools
    pip install pandas matplotlib
"""

import os
import subprocess
import sys
import pandas as pd
import matplotlib.pyplot as plt

# ── Paths ─────────────────────────────────────────────────────────────────────
ROOT       = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR   = os.path.join(ROOT, "data")
WORK_DIR   = os.path.join(ROOT, "work")
OUTPUT_DIR = os.path.join(ROOT, "output")

os.makedirs(WORK_DIR,   exist_ok=True)
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── Helper ────────────────────────────────────────────────────────────────────
def run(cmd: str) -> None:
    print(f"\n$ {cmd}\n")
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")


# ── Step 1: Extract start codons (+ strand, TSL 1) from GTF ──────────────────
print("=== Step 1: extract start codons ===")
run(
    f"grep '\tstart_codon\t.*\t+\t.*transcript_support_level \"1\"'"
    f" {DATA_DIR}/gencode.gtf"
    f" | sed -e 's/\t[^\t]*transcript_id \"\\([^\"]*\\)\".*$/\t\\1/g'"
    f" > {WORK_DIR}/gencode-start.gtf"
)

# ── Step 2: Extract + strand exons from GTF ───────────────────────────────────
print("=== Step 2: extract + strand exons ===")
run(
    f"grep '\texon\t.*\t+\t' {DATA_DIR}/gencode.gtf"
    f" | sed -e 's/\t[^\t]*transcript_id \"\\([^\"]*\\)\".*$/\t\\1/g'"
    f" > {WORK_DIR}/gencode-plusexon.gtf"
)

# ── Step 3: Find exons containing start codons ────────────────────────────────
print("=== Step 3: exons containing start codons ===")
run(
    f"bedtools intersect"
    f" -a {WORK_DIR}/gencode-start.gtf"
    f" -b {WORK_DIR}/gencode-plusexon.gtf"
    f" -wa -wb"
    f" | awk -F'\\t' -v OFS='\\t'"
    f"   '$9 == $18 {{ print $10, $13-1, $14, $18, $4-1, $16; }}'"
    f" | sort -k1,1 -k2,3n -k4,4"
    f" > {WORK_DIR}/gencode-exons-containing-startcodon.bed"
)

# ── Step 4: Filter RPF BAM (+ strand, read length >= 25 nt) ──────────────────
print("=== Step 4: filter RPF BAM ===")
run(
    f"(samtools view -H {DATA_DIR}/RPF-siLuc.bam;"
    f" samtools view -F20 {DATA_DIR}/RPF-siLuc.bam"
    f" | bioawk -c sam '{{ if (length($seq) >= 25) print $0; }}')"
    f" | samtools view -b -o {WORK_DIR}/filtered-RPF-siLuc.bam"
)

# ── Step 5: Compute 5'-end coverage (bedgraph) ────────────────────────────────
print("=== Step 5: 5'-end coverage ===")
run(
    f"bedtools genomecov"
    f" -ibam {WORK_DIR}/filtered-RPF-siLuc.bam"
    f" -bg -5"
    f" > {WORK_DIR}/fivepcounts-RPF-siLuc.bed"
)

# ── Step 6: Keep only 5'-ends inside start-codon exons ───────────────────────
print("=== Step 6: intersect with start-codon exons ===")
run(
    f"bedtools intersect"
    f" -a {WORK_DIR}/fivepcounts-RPF-siLuc.bed"
    f" -b {WORK_DIR}/gencode-exons-containing-startcodon.bed"
    f" -wa -wb -nonamecheck"
    f" > {WORK_DIR}/fivepcounts-filtered-RPF-siLuc.txt"
)

# ── Step 7: Plot Figure S5A ───────────────────────────────────────────────────
print("=== Step 7: plot Figure S5A ===")

# Column layout after -wa -wb intersect:
#  0  read_chr | 1 read_start (0-based 5' end) | 2 read_end | 3 count
#  4  exon_chr | 5 exon_start                  | 6 exon_end
#  7  transcript_id | 8 startcodon_pos (0-based) | 9 strand
COLS = [
    "read_chr", "read_start", "read_end", "count",
    "exon_chr", "exon_start", "exon_end",
    "transcript_id", "startcodon_pos", "strand",
]

WINDOW = (-50, 100)   # nt window around start codon
CHUNK  = 500_000

profile = pd.Series(dtype=float)

for chunk in pd.read_csv(
    f"{WORK_DIR}/fivepcounts-filtered-RPF-siLuc.txt",
    sep="\t", header=None, names=COLS, chunksize=CHUNK,
):
    chunk["rel_pos"] = chunk["read_start"] - chunk["startcodon_pos"]
    chunk = chunk[
        (chunk["rel_pos"] >= WINDOW[0]) & (chunk["rel_pos"] <= WINDOW[1])
    ]
    agg = chunk.groupby("rel_pos")["count"].sum()
    profile = profile.add(agg, fill_value=0)

profile = profile.sort_index()

fig, ax = plt.subplots(figsize=(12, 4))
ax.bar(profile.index, profile.values, width=1, color="steelblue", alpha=0.85)
ax.axvline(0, color="red", linestyle="--", linewidth=1.5, label="Start codon (A of AUG)")
ax.set_xlabel("Position relative to start codon (nt)")
ax.set_ylabel("RPF 5'-end count")
ax.set_title("RPF 5'-end distribution around start codons\n(+ strand, RPF-siLuc, TSL-1 transcripts)")
ax.legend()
plt.tight_layout()

out_path = os.path.join(OUTPUT_DIR, "figure_S5A.png")
plt.savefig(out_path, dpi=150)
print(f"\nFigure saved → {out_path}")
plt.show()
