#!/usr/bin/env python3
"""
Analysis 2 – Ribosome Triplet Phasing Verification
  (extends mission2 meta-gene plot with frame quantification + FFT)

Data required:
  mission2/work/fivepcounts-filtered-RPF-siLuc.txt  (already produced by mission2 pipeline)

Outputs (saved to ../output/):
  analysis2_metagene.png     – RPF 5'-end meta-gene profile around start codon
  analysis2_phasing.png      – Frame 0/1/2 bar chart with chi-square p-value
  analysis2_fft.png          – FFT power spectrum + autocorrelation
  analysis2_stats.txt        – Numerical summary
"""

import os, sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats, signal

# ── Paths ──────────────────────────────────────────────────────────────────────
ROOT        = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR  = os.path.join(ROOT, "results")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# The filtered 5'-end counts file is produced by mission2 pipeline
DATA_FILE = os.path.join(ROOT, "..", "..", "mission2", "work",
                         "fivepcounts-filtered-RPF-siLuc.txt")

if not os.path.exists(DATA_FILE):
    sys.exit(f"ERROR: data file not found:\n  {DATA_FILE}\n"
             "Run the mission2 pipeline first.")

# ── Parameters ─────────────────────────────────────────────────────────────────
WINDOW = (-50, 100)   # nt window around start codon
CHUNK  = 500_000

COLS = [
    "read_chr", "read_start", "read_end", "count",
    "exon_chr", "exon_start", "exon_end",
    "transcript_id", "startcodon_pos", "strand",
]

# ── 1. Build meta-gene profile ─────────────────────────────────────────────────
print("Building meta-gene profile …")
profile = pd.Series(dtype=float)

for chunk in pd.read_csv(DATA_FILE, sep="\t", header=None, names=COLS,
                         chunksize=CHUNK):
    chunk["rel_pos"] = chunk["read_start"] - chunk["startcodon_pos"]
    chunk = chunk[(chunk["rel_pos"] >= WINDOW[0]) &
                  (chunk["rel_pos"] <= WINDOW[1])]
    agg = chunk.groupby("rel_pos")["count"].sum()
    profile = profile.add(agg, fill_value=0)

profile = profile.sort_index()
full_index = pd.RangeIndex(WINDOW[0], WINDOW[1] + 1)
profile = profile.reindex(full_index, fill_value=0)

total_reads = int(profile.sum())
print(f"  Positions: {len(profile)}  |  Total accumulated reads: {total_reads:,}")

# ── 2. Plot 1 – meta-gene bar chart ───────────────────────────────────────────
fig, ax = plt.subplots(figsize=(13, 4))
colors = np.where(profile.index % 3 == 0, "#e34a33", "#9ecae1")
ax.bar(profile.index, profile.values, width=1.0,
       color=colors, edgecolor="none", alpha=0.9)
ax.axvline(0, color="black", linewidth=1.5, linestyle="--", label="Start codon (A of AUG)")
ax.set_xlabel("Position relative to start codon (nt)", fontsize=12)
ax.set_ylabel("RPF 5'-end read count", fontsize=12)
ax.set_title("Meta-gene RPF 5'-end profile around start codon\n"
             "(+ strand, RPF-siLuc, TSL-1; red bars = frame 0 positions)",
             fontsize=11)
ax.legend(fontsize=10)
plt.tight_layout()
out1 = os.path.join(OUTPUT_DIR, "analysis2_metagene.png")
plt.savefig(out1, dpi=150)
plt.close()
print(f"[saved] {out1}")

# ── 3. Frame phasing quantification (CDS region: pos 0 → 99) ─────────────────
# Ribosome occupies the start codon at pos 0.
# Frame 0 = in-frame codons: 0,3,6,9 …
# Frame 1 = +1 shift: 1,4,7,10 …
# Frame 2 = +2 shift: 2,5,8,11 …
cds = profile[(profile.index >= 0) & (profile.index <= 99)]
frame_counts = np.array([
    cds[cds.index % 3 == f].sum() for f in range(3)
])
frame_labels = ["Frame 0 (in-frame)", "Frame 1 (+1 shift)", "Frame 2 (+2 shift)"]
frame_pct    = 100 * frame_counts / frame_counts.sum()

print("\nFrame distribution (0–99 nt after start codon):")
for lbl, cnt, pct in zip(frame_labels, frame_counts, frame_pct):
    print(f"  {lbl}: {cnt:,.0f} reads  ({pct:.1f}%)")

chi2, p_chi2 = stats.chisquare(frame_counts)
print(f"\nChi-square test (expected uniform 33.3% each):")
print(f"  χ² = {chi2:.2f},  p = {p_chi2:.4e}")

# Plot 2 – frame phasing bar chart
fig, ax = plt.subplots(figsize=(6, 5))
bars = ax.bar(frame_labels, frame_pct,
              color=["#e34a33", "#9ecae1", "#9ecae1"],
              edgecolor="black", linewidth=0.8)
ax.axhline(100/3, color="gray", linewidth=1.2, linestyle="--", label="Expected (uniform)")
ax.set_ylabel("% of CDS reads", fontsize=12)
ax.set_title(f"Ribosome Triplet Phasing (0–99 nt after AUG)\n"
             f"χ² = {chi2:.1f},  p = {p_chi2:.2e}", fontsize=11)
ax.set_ylim(0, max(frame_pct) * 1.2)
ax.legend(fontsize=9)
for bar, pct in zip(bars, frame_pct):
    ax.text(bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.5, f"{pct:.1f}%",
            ha="center", va="bottom", fontsize=10)
plt.tight_layout()
out2 = os.path.join(OUTPUT_DIR, "analysis2_phasing.png")
plt.savefig(out2, dpi=150)
plt.close()
print(f"[saved] {out2}")

# ── 4. FFT + Autocorrelation ───────────────────────────────────────────────────
y = profile.values.astype(float)
y_centered = y - y.mean()

# FFT (one-sided power spectrum)
fft_power = np.abs(np.fft.rfft(y_centered)) ** 2
freqs     = np.fft.rfftfreq(len(y))
periods   = np.where(freqs > 0, 1.0 / freqs, np.inf)

# Find the peak near period = 3
mask_3 = (periods >= 2.5) & (periods <= 3.5)
peak_power  = fft_power[mask_3].max()
total_power = fft_power[1:].sum()  # exclude DC
period3_fraction = peak_power / total_power * 100
print(f"\nFFT: period-3 signal captures {period3_fraction:.1f}% of total spectral power")

# Autocorrelation (lag 0 – 30)
acf_full = np.correlate(y_centered, y_centered, mode="full")
acf      = acf_full[len(acf_full) // 2:]
acf      = acf / acf[0]         # normalise to 1 at lag 0
lags     = np.arange(len(acf))

# Prominent peak at lag 3 confirms periodicity
print(f"Autocorrelation at lag 3: {acf[3]:.4f}"
      f"  (lag 1: {acf[1]:.4f}, lag 2: {acf[2]:.4f})")

# Plot 3 – FFT + autocorrelation side by side
fig, axes = plt.subplots(1, 2, figsize=(13, 4))

# Left: FFT
ax = axes[0]
valid = freqs > 0
ax.plot(1.0 / freqs[valid], fft_power[valid], color="steelblue", linewidth=1)
ax.axvline(3, color="#e34a33", linewidth=1.5, linestyle="--",
           label="Period = 3 nt (triplet codon)")
ax.set_xlim(1, 20)
ax.set_xlabel("Period (nt)", fontsize=12)
ax.set_ylabel("FFT Power", fontsize=12)
ax.set_title(f"FFT Power Spectrum\n(period-3 power: {period3_fraction:.1f}% of total)", fontsize=11)
ax.legend(fontsize=9)

# Right: Autocorrelation
ax = axes[1]
ax.bar(lags[:31], acf[:31], width=0.8,
       color=["#e34a33" if lag % 3 == 0 else "#9ecae1" for lag in range(31)],
       edgecolor="none")
ax.axhline(0, color="black", linewidth=0.5)
ax.set_xlabel("Lag (nt)", fontsize=12)
ax.set_ylabel("Autocorrelation", fontsize=12)
ax.set_title(f"Autocorrelation (lag 0–30)\n"
             f"lag-3 ACF = {acf[3]:.3f}", fontsize=11)
ax.set_xlim(0, 30)

plt.tight_layout()
out3 = os.path.join(OUTPUT_DIR, "analysis2_fft.png")
plt.savefig(out3, dpi=150)
plt.close()
print(f"[saved] {out3}")

# ── 5. Write numerical summary ─────────────────────────────────────────────────
stats_out = os.path.join(OUTPUT_DIR, "analysis2_stats.txt")
with open(stats_out, "w") as fh:
    fh.write("=== Analysis 2: Ribosome Triplet Phasing ===\n\n")
    fh.write(f"Window: {WINDOW[0]} to {WINDOW[1]} nt around start codon\n")
    fh.write(f"Total reads in profile: {total_reads:,}\n\n")
    fh.write("Frame phasing (0–99 nt, CDS region):\n")
    for lbl, cnt, pct in zip(frame_labels, frame_counts, frame_pct):
        fh.write(f"  {lbl}: {cnt:,.0f}  ({pct:.1f}%)\n")
    fh.write(f"\nChi-square test: χ² = {chi2:.4f}, p = {p_chi2:.6e}\n")
    fh.write(f"\nFFT period-3 power fraction: {period3_fraction:.2f}%\n")
    fh.write(f"Autocorrelation at lag 3: {acf[3]:.4f}\n")
print(f"[saved] {stats_out}")
print("\nDone – Analysis 2 complete.")
