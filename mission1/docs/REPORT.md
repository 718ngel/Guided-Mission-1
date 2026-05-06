# Guided Mission 1 — Project Report

**Course:** 생물정보학 및 실습 1 (Bioinformatics and Lab 1)  
**Institution:** Seoul National University  
**Date:** 2026-04-27  

---

## 1. Project Overview

This project analyzes the role of **Lin28a** in translational regulation using three types of sequencing data:

| Data type | Description |
|-----------|-------------|
| **CLIP-seq** | Identifies RNA targets bound by Lin28a protein |
| **RNA-seq** | Measures mRNA abundance (control and siLin28a knockdown) |
| **RPF-seq** (Ribosome Footprint) | Measures actively translated regions (proxy for translation efficiency) |

The goal is to reproduce the scatter plots from the original paper (Figure 4D and Figure 5B), showing which Lin28a-bound mRNAs show altered ribosome density upon Lin28a knockdown.

---

## 2. Pipeline Summary

```
BAM files (6 × ~1 GB)
        ↓
   featureCounts  ←  GENCODE vM27 annotation (gencode.gtf, 842 MB)
        ↓
  read-counts.txt  (55,359 genes × 6 samples)
        ↓
  Compute metrics:
    clip_enrichment = CLIP-35L33G / RNA-control
    rden_change     = (RPF-siLin28a / RNA-siLin28a) / (RPF-siLuc / RNA-siLuc)
        ↓
  scatter_basic.png          ← Step 6
  scatter_localization.png   ← Step 7  (colored by UniProt subcellular location)
```

---

## 3. Data Integrity (Step 4)

All 12 input files passed MD5 checksum verification.

| File | Size | MD5 |
|------|------|-----|
| CLIP-35L33G.bam | 1.3 GB | ✅ PASS |
| CLIP-35L33G.bam.bai | 3.0 MB | ✅ PASS |
| RNA-control.bam | 957 MB | ✅ PASS |
| RNA-control.bam.bai | 2.2 MB | ✅ PASS |
| RNA-siLin28a.bam | 1.2 GB | ✅ PASS |
| RNA-siLin28a.bam.bai | 2.6 MB | ✅ PASS |
| RNA-siLuc.bam | 936 MB | ✅ PASS |
| RNA-siLuc.bam.bai | 2.5 MB | ✅ PASS |
| RPF-siLin28a.bam | 703 MB | ✅ PASS |
| RPF-siLin28a.bam.bai | 2.4 MB | ✅ PASS |
| RPF-siLuc.bam | 1.0 GB | ✅ PASS |
| RPF-siLuc.bam.bai | 2.5 MB | ✅ PASS |

---

## 4. featureCounts Results (Step 5)

**Tool:** featureCounts v2.0.1 (Subread package)  
**Annotation:** GENCODE mouse vM27 — 841,952 features → 55,359 meta-features (genes)  
**Mode:** Single-end reads, meta-feature level (gene)

### Read Assignment Summary

| Sample | Total alignments | Assigned | % Assigned | Multi-mapping | No feature |
|--------|-----------------|----------|------------|---------------|------------|
| CLIP-35L33G | 38,880,853 | 13,630,945 | **35.1%** | 17,781,701 | 4,669,517 |
| RNA-control | 24,971,058 | 5,200,836 | **20.8%** | 5,511,892 | 1,395,594 |
| RNA-siLin28a | 35,108,178 | 12,338,280 | **35.1%** | 18,072,245 | 2,035,101 |
| RNA-siLuc | 28,117,241 | 9,748,787 | **34.7%** | 13,999,092 | 1,661,248 |
| RPF-siLin28a | 30,095,551 | 10,694,162 | **35.5%** | 16,675,170 | 720,876 |
| RPF-siLuc | 42,002,435 | 12,932,285 | **30.8%** | 26,408,669 | 752,105 |

> **Note:** The high multi-mapping rate (~40–60%) is expected for CLIP-seq and RPF-seq data, as many reads map to repetitive elements or highly similar gene family members. RNA-control has the lowest assignment rate (20.8%) mainly due to a high unmapped fraction (12.6M reads), suggesting lower sequencing quality or a different library preparation protocol.

### Expressed Genes per Sample

| Sample | Genes with read count > 0 |
|--------|--------------------------|
| CLIP-35L33G | 27,594 |
| RNA-control | 27,582 |
| RNA-siLin28a | 27,316 |
| RNA-siLuc | 25,425 |
| RPF-siLin28a | 21,985 |
| RPF-siLuc | 23,253 |

---

## 5. Quantitative Analysis (Steps 6 & 7)

### Metric Definitions

- **CLIP Enrichment** = CLIP-35L33G counts / RNA-control counts  
  → High values indicate the gene's mRNA is preferentially bound by Lin28a.

- **Ribosome Density Change (rden_change)** = (RPF-siLin28a / RNA-siLin28a) / (RPF-siLuc / RNA-siLuc)  
  → Values > 1 (log2 > 0): translation increases when Lin28a is knocked down.  
  → Values < 1 (log2 < 0): translation decreases when Lin28a is knocked down.

Genes with zero counts in any sample produce undefined (Inf/NaN) log2 values and are excluded from plots.

**Genes with finite log2 values (used for plotting): 16,201 out of 55,359**

### Distribution of log2(CLIP Enrichment)

| Statistic | Value |
|-----------|-------|
| Mean | 0.74 |
| Median | 0.80 |
| Std dev | 1.64 |
| Min | −9.61 |
| Max | 11.38 |

> The positive mean and median indicate that Lin28a is generally enriched on most expressed mRNAs, consistent with its role as a broad RNA-binding protein.

### Distribution of log2(Ribosome Density Change)

| Statistic | Value |
|-----------|-------|
| Mean | −0.79 |
| Median | −0.83 |
| Std dev | 1.01 |
| Min | −6.29 |
| Max | 6.00 |

> The negative mean/median shows that Lin28a knockdown **globally reduces translation efficiency**. This is the expected result: Lin28a promotes translation of its target mRNAs, so removing it lowers ribosome occupancy on average.

---

## 6. Protein Localization Analysis (Step 7)

Subcellular localization data was obtained from UniProt (via `mouselocalization-20210507.txt`, 9,523 mouse genes).

### Localization Breakdown (among 16,279 plotted genes)

| Localization | Gene count | % of plotted |
|---|---|---|
| Unknown / not annotated | 9,162 | 56.3% |
| Nucleus | 3,501 | 21.5% |
| Integral membrane | 1,897 | 11.7% |
| Cytoplasm | 1,719 | 10.6% |

### Mean log2 Values by Localization

| Localization | Mean log2(CLIP enrich.) | Mean log2(rden change) |
|---|---|---|
| Cytoplasm | 0.67 | −0.94 |
| Nucleus | 0.77 | −1.03 |
| Integral membrane | **1.43** | **−0.27** |

> **Key observation:** Integral membrane proteins show **higher CLIP enrichment** (mean log2 = 1.43 vs ~0.7 for others) and a **much smaller decrease in ribosome density** (mean log2 = −0.27 vs ~−1.0 for others) upon Lin28a knockdown. This suggests that Lin28a may have a particularly important role in the translational regulation of mRNAs encoding membrane-targeted proteins — consistent with findings in the original paper linking Lin28a to ER-associated translation.

---

## 7. Output Files

| File | Description |
|------|-------------|
| `data/read-counts.txt` | Raw gene-level read counts (55,359 genes × 6 samples) — generated locally |
| `data/read-counts.txt.summary` | featureCounts assignment summary |
| `output/scatter_basic.png` | log2(CLIP enrichment) vs log2(rden change), all genes |
| `output/scatter_localization.png` | Same scatter, colored by subcellular localization |
| `step6_analyze.py` | Script for basic analysis and scatter plot |
| `step7_localization_plot.py` | Script for localization-colored scatter plot |

---

## 8. How to Reproduce

```bash
# Prerequisites: conda environment with featureCounts, Python with pandas/matplotlib/numpy

# Step 5 — Count reads
cd data/
conda run -n subread_env featureCounts -a gencode.gtf -o read-counts.txt \
    CLIP-35L33G.bam RNA-control.bam RNA-siLin28a.bam RNA-siLuc.bam \
    RPF-siLin28a.bam RPF-siLuc.bam

# Step 6 — Basic scatter plot
cd ..
python3 step6_analyze.py

# Step 7 — Localization-colored scatter plot
python3 step7_localization_plot.py
```

Output figures are saved in `output/`.
