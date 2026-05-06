#!/usr/bin/env python3
"""
Mission 2 – RPF 5'-end distribution around start codons
Full pipeline: shell steps + Figure S5A plot.

Usage:  python project_code/pipeline.py
"""

import os
import subprocess
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")          # headless rendering (no display needed)
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ── Tool paths (Homebrew samtools takes priority over broken conda version) ──
SAMTOOLS = "/usr/local/bin/samtools"
BEDTOOLS  = "/Users/angela/opt/anaconda3/bin/bedtools"
BIOAWK    = "/Users/angela/opt/anaconda3/bin/bioawk"

ROOT       = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR   = os.path.join(ROOT, "data")
WORK_DIR   = os.path.join(ROOT, "work")
RESULTS    = os.path.join(ROOT, "results")

for d in (WORK_DIR, RESULTS):
    os.makedirs(d, exist_ok=True)

GTF     = os.path.join(DATA_DIR, "gencode.gtf")
RPF_BAM = os.path.join(DATA_DIR, "RPF-siLuc.bam")

# ── Helper ────────────────────────────────────────────────────────────────────
def run(cmd: str, desc: str = "") -> None:
    label = f"  [{desc}]" if desc else ""
    print(f"\n>>> {label}\n$ {cmd}\n")
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 1 – Extract start codons (+ strand, TSL 1)
# ═══════════════════════════════════════════════════════════════════════════════
START_GTF = os.path.join(WORK_DIR, "gencode-start.gtf")
print("=" * 60)
print("STEP 1 – Extract start codons (+ strand, TSL 1)")
print("=" * 60)
run(
    f"grep $'\\tstart_codon\\t' {GTF}"
    f" | grep $'\\t+\\t'"
    f" | grep 'transcript_support_level \"1\"'"
    f" | sed -e 's/\\t[^\\t]*transcript_id \"\\([^\"]*\\)\".*$/\\t\\1/g'"
    f" > {START_GTF}",
    "gencode-start.gtf"
)
run(f"wc -l {START_GTF}", "line count")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 2 – Extract + strand exons
# ═══════════════════════════════════════════════════════════════════════════════
EXON_GTF = os.path.join(WORK_DIR, "gencode-plusexon.gtf")
print("\n" + "=" * 60)
print("STEP 2 – Extract + strand exons")
print("=" * 60)
run(
    f"grep $'\\texon\\t' {GTF}"
    f" | grep $'\\t+\\t'"
    f" | sed -e 's/\\t[^\\t]*transcript_id \"\\([^\"]*\\)\".*$/\\t\\1/g'"
    f" > {EXON_GTF}",
    "gencode-plusexon.gtf"
)
run(f"wc -l {EXON_GTF}", "line count")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 3 – Find exons containing a start codon (same transcript)
# ═══════════════════════════════════════════════════════════════════════════════
SC_EXON_BED = os.path.join(WORK_DIR, "gencode-exons-containing-startcodon.bed")
print("\n" + "=" * 60)
print("STEP 3 – Exons containing start codons")
print("=" * 60)
run(
    f"{BEDTOOLS} intersect"
    f" -a {START_GTF} -b {EXON_GTF} -wa -wb"
    f" | awk -F'\\t' -v OFS='\\t'"
    f"   '$9 == $18 {{ print $10, $13-1, $14, $18, $4-1, $16 }}'"
    f" | sort -k1,1 -k2,3n -k4,4"
    f" > {SC_EXON_BED}",
    "gencode-exons-containing-startcodon.bed"
)
run(f"wc -l {SC_EXON_BED}", "line count")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 4 – Filter RPF BAM (+ strand, read length ≥ 25 nt)
# ═══════════════════════════════════════════════════════════════════════════════
FILT_BAM = os.path.join(WORK_DIR, "filtered-RPF-siLuc.bam")
print("\n" + "=" * 60)
print("STEP 4 – Filter RPF BAM (+ strand, length ≥ 25 nt)")
print("=" * 60)
run(
    f"({SAMTOOLS} view -H {RPF_BAM};"
    f" {SAMTOOLS} view -F20 {RPF_BAM}"
    f" | {BIOAWK} -c sam '{{ if (length($seq) >= 25) print $0 }}')"
    f" | {SAMTOOLS} view -b -o {FILT_BAM}",
    "filtered-RPF-siLuc.bam"
)
run(f"ls -lh {FILT_BAM}", "file size")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 5 – 5'-end coverage (bedgraph)
# ═══════════════════════════════════════════════════════════════════════════════
FIVEP_BED = os.path.join(WORK_DIR, "fivepcounts-RPF-siLuc.bed")
print("\n" + "=" * 60)
print("STEP 5 – 5'-end coverage (bedgraph)")
print("=" * 60)
run(
    f"{BEDTOOLS} genomecov -ibam {FILT_BAM} -bg -5 > {FIVEP_BED}",
    "fivepcounts-RPF-siLuc.bed"
)
run(f"wc -l {FIVEP_BED}", "line count")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 6 – Intersect 5'-ends with start-codon exons
# ═══════════════════════════════════════════════════════════════════════════════
INTER_TXT = os.path.join(WORK_DIR, "fivepcounts-filtered-RPF-siLuc.txt")
print("\n" + "=" * 60)
print("STEP 6 – Intersect 5'-ends with start-codon exons")
print("=" * 60)
run(
    f"{BEDTOOLS} intersect"
    f" -a {FIVEP_BED} -b {SC_EXON_BED}"
    f" -wa -wb -nonamecheck"
    f" > {INTER_TXT}",
    "fivepcounts-filtered-RPF-siLuc.txt"
)
run(f"wc -l {INTER_TXT}", "line count")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 7 – Compute profile & plot Figure S5A
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("STEP 7 – Build profile and plot Figure S5A")
print("=" * 60)

COLS = [
    "read_chr", "read_start", "read_end", "count",
    "exon_chr", "exon_start", "exon_end",
    "transcript_id", "startcodon_pos", "strand",
]
WINDOW = (-50, 150)
CHUNK  = 500_000

profile = pd.Series(dtype=float)
total_rows = 0

print("  Loading intersection file in chunks …")
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
print(f"  Total rows processed : {total_rows:,}")
print(f"  Positions in window  : {len(profile)}")
print(f"  Peak position        : {profile.idxmax()} nt (count={int(profile.max()):,})")

# ── Summary table (save to results/) ─────────────────────────────────────────
summary = pd.DataFrame({
    "position_nt": profile.index.astype(int),
    "rpf_5end_count": profile.values.astype(int),
})
csv_path = os.path.join(RESULTS, "rpf_profile_around_startcodon.csv")
summary.to_csv(csv_path, index=False)
print(f"\n  Profile CSV saved → {csv_path}")

# ── Figure ────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(14, 5))

colors = ["#d62728" if i == 0 else "#1f77b4" for i in profile.index]
ax.bar(profile.index, profile.values, width=1, color=colors, alpha=0.85)

ax.axvline(0,   color="red",    linestyle="--", linewidth=1.5, label="Start codon A (pos 0)")
ax.axvline(-15, color="orange", linestyle=":",  linewidth=1.2, label="−15 nt (ribosome P-site)")

ax.set_xlim(WINDOW)
ax.set_xlabel("Position relative to start codon A (nt)", fontsize=12)
ax.set_ylabel("RPF 5'-end count", fontsize=12)
ax.set_title(
    "RPF 5'-end distribution around start codons\n"
    "(+ strand, RPF-siLuc, GENCODE M27, TSL-1 transcripts)",
    fontsize=13,
)
ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(3))
ax.legend(fontsize=10)
ax.grid(axis="y", alpha=0.3)

plt.tight_layout()

png_path = os.path.join(RESULTS, "figure_S5A_rpf_startcodon.png")
pdf_path = os.path.join(RESULTS, "figure_S5A_rpf_startcodon.pdf")
plt.savefig(png_path, dpi=200)
plt.savefig(pdf_path)
print(f"  Figure PNG saved → {png_path}")
print(f"  Figure PDF saved → {pdf_path}")

print("\n" + "=" * 60)
print("Pipeline complete. Results in:  results/")
print("=" * 60)
