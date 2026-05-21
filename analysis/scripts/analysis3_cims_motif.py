#!/usr/bin/env python3
"""
Analysis 3 – CLIP-seq CIMS Shannon Entropy + Cross-linking Site Characterisation
             + k-mer Motif Enrichment at Lin28a Binding Sites

Key improvement vs. mission3 cims_analysis.py:
  - Adds deletion-rate analysis (true CIMS biomarker)
  - Identifies high-confidence cross-linking sites by combining entropy + deletion rate
  - Extracts read sequences at CIMS sites from BAM files via samtools
  - Computes 4-mer frequency enrichment (CIMS vs. background reads)

Data required:
  mission3/project/data/work/CLIP-{let7g,let7f1,let7d}-gene.pileup
  mission3/project/data/work/CLIP-{let7g,let7f1,let7d}.bam

Outputs (saved to ../output/):
  analysis3_entropy_deletion.png   – entropy + deletion-rate per position
  analysis3_cims_sites.png         – CIMS site VENN / scatter
  analysis3_kmer_enrichment.png    – top enriched 4-mers at CIMS vs background reads
  analysis3_cims_sites.csv         – table of detected CIMS sites
"""

import os, re, io, sys, math, subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from collections import Counter
from pathlib import Path
from scipy import stats

# ── Paths ──────────────────────────────────────────────────────────────────────
ROOT        = Path(__file__).resolve().parent.parent
WORK3       = ROOT.parent / "mission3" / "project" / "data" / "work"
OUTPUT_DIR  = ROOT / "output"
OUTPUT_DIR.mkdir(exist_ok=True)

GENES = {
    "let7g":  ("chr9",  106056039, 106056126, "+"),
    "let7f1": ("chr13",  48691305,  48691393, "-"),
    "let7d":  ("chr13",  48689488,  48689590, "-"),
}

# Thresholds for CIMS site calling
ENTROPY_QUANTILE  = 0.90   # top 10% of entropy → candidate CIMS
DEL_RATE_MIN      = 0.05   # ≥ 5% deletion rate required as confirmatory signal
MIN_DEPTH         = 5      # ignore positions with fewer reads

# ── Pileup parsing helpers ─────────────────────────────────────────────────────
_STRIP = re.compile(r'[<>$*#^]|\+\d+[ACGTNacgtn]+|-\d+[ACGTNacgtn]+')

def clean_basereads(s: str) -> str:
    return _STRIP.sub("", s)

def count_bases(s: str) -> Counter:
    return Counter(c.upper() for c in s if c.upper() in "ACGT.,")

def shannon_entropy(counts: Counter) -> float:
    total = sum(counts.values())
    if total == 0:
        return 0.0
    return -sum((n / total) * math.log2(n / total) for n in counts.values() if n > 0)

def deletion_rate(raw: str, depth: int) -> float:
    return raw.count("*") / depth if depth > 0 else 0.0

def mismatch_rate(counts: Counter) -> float:
    total = sum(counts.values())
    if total == 0:
        return 0.0
    mismatches = sum(v for k, v in counts.items() if k in "ACGT")
    return mismatches / total

# ── 1. Process pileup files ────────────────────────────────────────────────────
print("Processing pileup files …")
all_dfs = []

for name, (chrom, start, end, strand) in GENES.items():
    pileup_path = WORK3 / f"CLIP-{name}-gene.pileup"
    if not pileup_path.exists():
        print(f"  WARNING: {pileup_path.name} not found – skipping")
        continue

    cols = ["chrom", "pos", "_ref", "depth", "basereads", "quals"]
    df = pd.read_csv(pileup_path, sep="\t", names=cols)

    df["clean_reads"] = df["basereads"].apply(clean_basereads)
    df["base_counts"] = df["clean_reads"].apply(count_bases)
    df["entropy"]     = df["base_counts"].apply(shannon_entropy)
    df["del_rate"]    = df.apply(lambda r: deletion_rate(r["basereads"], r["depth"]), axis=1)
    df["mis_rate"]    = df["base_counts"].apply(mismatch_rate)
    df["gene"]        = name
    df["strand"]      = strand

    # Filter out very low-coverage positions
    df = df[df["depth"] >= MIN_DEPTH].copy()
    all_dfs.append(df)
    print(f"  {name}: {len(df)} positions  (depth ≥ {MIN_DEPTH})")

if not all_dfs:
    sys.exit("No pileup data loaded – check mission3 paths.")

combined = pd.concat(all_dfs, ignore_index=True)
print(f"\nTotal positions: {len(combined)}")

# ── 2. CIMS site calling ───────────────────────────────────────────────────────
ent_thresh = combined["entropy"].quantile(ENTROPY_QUANTILE)
print(f"\nEntropy threshold (top {100*(1-ENTROPY_QUANTILE):.0f}%): {ent_thresh:.4f}")

# CIMS site = high entropy AND confirmed by deletion rate
combined["is_high_entropy"] = combined["entropy"] >= ent_thresh
combined["is_cims"]         = combined["is_high_entropy"] & \
                              (combined["del_rate"] >= DEL_RATE_MIN)

n_ent  = combined["is_high_entropy"].sum()
n_cims = combined["is_cims"].sum()
print(f"High-entropy sites (entropy ≥ {ent_thresh:.3f}): {n_ent}")
print(f"Confirmed CIMS sites (+ del_rate ≥ {DEL_RATE_MIN}): {n_cims}")

cims_df = combined[combined["is_cims"]].copy()
bg_df   = combined[~combined["is_high_entropy"]].copy()

# Mann-Whitney test: deletion rate at CIMS vs. background
u_stat, p_mwu = stats.mannwhitneyu(cims_df["del_rate"], bg_df["del_rate"],
                                    alternative="greater")
print(f"\nMann-Whitney U (deletion rate CIMS > bg): U={u_stat:.0f}, p={p_mwu:.4e}")

# ── 3. Plot – entropy + deletion rate per gene ─────────────────────────────────
gene_list = list(GENES.keys())
fig = plt.figure(figsize=(13, 3.5 * len(gene_list)))
gs  = gridspec.GridSpec(len(gene_list), 1, hspace=0.55)

for i, name in enumerate(gene_list):
    gdf = combined[combined["gene"] == name]
    if gdf.empty:
        continue
    ax1 = fig.add_subplot(gs[i])
    ax2 = ax1.twinx()

    ax1.fill_between(gdf["pos"], gdf["entropy"], alpha=0.6,
                     color="#3182bd", label="Shannon entropy")
    ax2.plot(gdf["pos"], gdf["del_rate"] * 100, color="#e34a33",
             linewidth=0.8, alpha=0.9, label="Deletion rate (%)")

    # Mark CIMS sites
    cims_here = gdf[gdf["is_cims"]]
    if not cims_here.empty:
        ax1.scatter(cims_here["pos"], cims_here["entropy"],
                    c="red", s=25, zorder=5, label="CIMS site")

    ax1.set_xlabel("Genomic position", fontsize=10)
    ax1.set_ylabel("Shannon entropy (bits)", color="#3182bd", fontsize=10)
    ax2.set_ylabel("Deletion rate (%)", color="#e34a33", fontsize=10)
    ax1.set_title(f"Lin28a CLIP CIMS – {name}  ({GENES[name][0]}:{GENES[name][1]}–{GENES[name][2]})",
                  fontsize=11)

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper right", fontsize=8)

out1 = str(OUTPUT_DIR / "analysis3_entropy_deletion.png")
plt.savefig(out1, dpi=150, bbox_inches="tight")
plt.close()
print(f"\n[saved] {out1}")

# ── 4. Entropy vs. deletion-rate scatter (all genes) ──────────────────────────
fig, ax = plt.subplots(figsize=(7, 6))
GENE_COLORS = {"let7g": "#1f77b4", "let7f1": "#ff7f0e", "let7d": "#2ca02c"}
for name in gene_list:
    gdf = combined[combined["gene"] == name]
    ax.scatter(gdf["entropy"], gdf["del_rate"] * 100,
               alpha=0.4, s=8, label=name, color=GENE_COLORS.get(name))

ax.axhline(DEL_RATE_MIN * 100, color="gray", linewidth=1, linestyle="--",
           label=f"del_rate = {DEL_RATE_MIN*100:.0f}%")
ax.axvline(ent_thresh, color="gray", linewidth=1, linestyle=":",
           label=f"entropy threshold = {ent_thresh:.2f}")
ax.set_xlabel("Shannon entropy (bits)", fontsize=12)
ax.set_ylabel("Deletion rate (%)", fontsize=12)
ax.set_title("Shannon entropy vs deletion rate at all CLIP positions\n"
             f"(n={len(combined):,} positions, depth ≥ {MIN_DEPTH})", fontsize=11)
ax.legend(fontsize=9)
plt.tight_layout()
out2 = str(OUTPUT_DIR / "analysis3_cims_scatter.png")
plt.savefig(out2, dpi=150)
plt.close()
print(f"[saved] {out2}")

# ── 5. k-mer enrichment from BAM reads ────────────────────────────────────────
# Strategy: extract read sequences from BAM files for CIMS regions vs. all reads,
# then compute 4-mer frequencies.
SAMTOOLS = "samtools"
K = 4

def extract_reads_from_bam(bam_path: Path, region: str) -> list[str]:
    """Return list of read sequences overlapping `region` (samtools view format)."""
    cmd = [SAMTOOLS, "view", str(bam_path), region]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        seqs = []
        for line in result.stdout.splitlines():
            fields = line.split("\t")
            if len(fields) >= 10:
                seqs.append(fields[9].upper())
        return seqs
    except Exception as e:
        print(f"    samtools failed for {region}: {e}")
        return []

def kmer_counts(sequences: list[str], k: int = 4) -> Counter:
    cnt = Counter()
    for seq in sequences:
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            if all(b in "ACGT" for b in kmer):
                cnt[kmer] += 1
    return cnt

print("\nExtracting reads for k-mer analysis …")
cims_reads = []
bg_reads   = []

for name, (chrom, gstart, gend, strand) in GENES.items():
    bam_path = WORK3 / f"CLIP-{name}.bam"
    if not bam_path.exists():
        print(f"  BAM not found: {bam_path.name}")
        continue

    gene_cims = cims_df[cims_df["gene"] == name]

    # CIMS reads: one samtools view call per CIMS site (1-based region)
    for _, row in gene_cims.iterrows():
        pos = int(row["pos"])
        region = f"{chrom}:{pos}-{pos}"
        reads = extract_reads_from_bam(bam_path, region)
        cims_reads.extend(reads)

    # Background reads: all reads from the gene region
    all_gene_region = f"{chrom}:{gstart}-{gend}"
    all_reads = extract_reads_from_bam(bam_path, all_gene_region)
    bg_reads.extend(all_reads)

print(f"  CIMS reads: {len(cims_reads):,}  |  Background reads: {len(bg_reads):,}")

if cims_reads and bg_reads:
    cims_kmers = kmer_counts(cims_reads, K)
    bg_kmers   = kmer_counts(bg_reads, K)

    # Normalise to frequency
    cims_total = sum(cims_kmers.values()) or 1
    bg_total   = sum(bg_kmers.values()) or 1

    all_kmers = set(cims_kmers) | set(bg_kmers)
    enrich = {}
    for km in all_kmers:
        f_cims = cims_kmers.get(km, 0) / cims_total
        f_bg   = bg_kmers.get(km, 0)   / bg_total
        if f_bg > 0:
            enrich[km] = f_cims / f_bg
        else:
            enrich[km] = np.inf

    enrich_series = pd.Series(enrich).replace([np.inf, -np.inf], np.nan).dropna().sort_values(ascending=False)
    top = enrich_series.head(20)

    # Known Lin28a binding motifs to highlight
    LIN28A_MOTIFS = ["GGAG", "GAGA", "AGGA", "GAAG", "AAGN"]
    highlight = {km for km in top.index if any(m[:K] in km or km in m for m in LIN28A_MOTIFS)}

    fig, ax = plt.subplots(figsize=(10, 5))
    bar_colors = ["#e34a33" if km in highlight else "#9ecae1" for km in top.index]
    ax.barh(top.index[::-1], top.values[::-1], color=bar_colors[::-1],
            edgecolor="black", linewidth=0.5)
    ax.axvline(1, color="gray", linewidth=1, linestyle="--", label="Enrichment = 1 (no change)")
    ax.set_xlabel(f"{K}-mer frequency enrichment (CIMS / background)", fontsize=11)
    ax.set_title(f"Top 20 enriched {K}-mers at Lin28a CIMS sites\n"
                 f"(red = overlaps known Lin28a binding motif GGAG/GAGA)", fontsize=11)
    ax.legend(fontsize=9)
    plt.tight_layout()
    out3 = str(OUTPUT_DIR / "analysis3_kmer_enrichment.png")
    plt.savefig(out3, dpi=150)
    plt.close()
    print(f"[saved] {out3}")

    print(f"\nTop 5 enriched {K}-mers at CIMS sites:")
    for km, val in enrich_series.head(5).items():
        flag = "** Lin28a motif **" if km in highlight else ""
        print(f"  {km}: {val:.2f}x  {flag}")
else:
    print("  Skipping k-mer plot (no reads extracted – check BAM paths or samtools)")

# ── 6. Save CIMS site table ────────────────────────────────────────────────────
cims_out_cols = ["gene", "chrom", "pos", "depth", "entropy", "del_rate", "mis_rate"]
cims_out      = cims_df[cims_out_cols].sort_values(["gene", "entropy"], ascending=[True, False])
csv_out = str(OUTPUT_DIR / "analysis3_cims_sites.csv")
cims_out.to_csv(csv_out, index=False)
print(f"\n[saved] {csv_out}")

# Print summary
print("\n=== CIMS site summary ===")
for name in gene_list:
    sub = cims_out[cims_out["gene"] == name]
    if not sub.empty:
        print(f"  {name}: {len(sub)} sites  "
              f"(max entropy={sub['entropy'].max():.3f}  "
              f"max del_rate={sub['del_rate'].max():.3f})")

print("\nDone – Analysis 3 complete.")
