# Analysis 2 — Ribosome Triplet Phasing Verification

**Course:** 생물정보학 및 실습 1  
**Institution:** Seoul National University  
**Reference:** `analysis/docs/REPORT.md` → Section "Analysis 2"  
**Script:** `analysis/analysis2/analysis2_triplet_phasing.py`

---

## Overview

This analysis verifies the **3-nucleotide (triplet) periodicity** of ribosome-protected fragments (RPF) around the start codon. Because the ribosome advances along the mRNA exactly one codon (3 nt) at a time, the 5′ ends of RPF reads should pile up at every third position relative to the AUG start codon. Confirming this triplet phasing serves two purposes:

1. **Quality control** — a strong period-3 signal indicates a high-quality ribosome profiling library with accurate ribosome footprinting.
2. **Biological validation** — the phasing directly mirrors the mechanics of translation elongation, confirming that the sequenced fragments are genuine ribosome footprints rather than random RNA degradation products.

The analysis extends the Mission 2 meta-gene pipeline by adding quantitative frame phasing (chi-square test), FFT power spectrum analysis, and autocorrelation to provide orthogonal evidence for the periodicity.

---

## Biological Background

The ribosome reads mRNA in non-overlapping triplets called **codons**. During translation elongation, the ribosome's P-site (peptidyl-tRNA site) occupies one codon while the A-site decodes the next. Because the ribosome moves exactly 3 nt per elongation cycle, the 5′ ends of RPF reads — which mark the upstream boundary of the ribosome footprint — should appear at positions 0, 3, 6, 9 … relative to the start codon AUG (these are **Frame 0** positions). Frames 1 and 2 (positions 1, 4, 7 … and 2, 5, 8 …) should have substantially fewer reads.

A chi-square test quantifies whether the frame 0:1:2 read distribution deviates significantly from the uniform null (33.3% each). FFT and autocorrelation independently confirm the period-3 signal in the full meta-gene profile.

---

## Input Data

| File | Description |
|------|-------------|
| `mission2/work/fivepcounts-filtered-RPF-siLuc.txt` | Tab-separated table of RPF 5′-end counts intersected with start-codon windows; produced by the Mission 2 `bedtools intersect` pipeline |

**Required columns (no header):**

| Column | Content |
|--------|---------|
| `read_chr` | Chromosome of the RPF read |
| `read_start` | Genomic start position of the RPF read (0-based) |
| `read_end` | Genomic end position |
| `count` | Number of reads at this position |
| `exon_chr` | Chromosome of the annotated exon |
| `exon_start` | Exon start |
| `exon_end` | Exon end |
| `transcript_id` | Ensembl transcript ID |
| `startcodon_pos` | Genomic position of the first base (A) of AUG |
| `strand` | Strand (`+` or `-`) |

> **Prerequisite:** The Mission 2 pipeline (bedtools filtering of RPF reads against TSL-1 start-codon windows) must have been run before this analysis.

---

## Step-by-Step Pipeline

### Step 1 — Load and build the meta-gene profile

**What:**  
Read `fivepcounts-filtered-RPF-siLuc.txt` in chunks (`chunksize = 500,000` rows) to stay within laptop RAM.  
For each chunk, compute the **relative position**:

```
rel_pos = read_start − startcodon_pos
```

Keep only positions within the window `[−50, +100]` nt, then sum read counts per relative position. Accumulate partial sums across chunks into a single `pandas.Series` indexed by `rel_pos`.

**Why chunked reading:**  
The file contains ~340,000 rows. Loading it all at once would require ~3–4 GB RAM. Chunked `groupby` + `Series.add(fill_value=0)` accumulation keeps peak memory under 500 MB.

**Window choice:**  
`−50` to `+100` nt captures 50 nt of 5′ UTR (where the ribosome is scanning/initiating) and the first 100 nt of coding sequence (where frame phasing is strongest and most interpretable).

---

### Step 2 — Plot the meta-gene bar chart (Figure 2A)

**What:**  
Draw a bar chart: x-axis = relative position (−50 to +100), y-axis = cumulative RPF 5′-end read count. Colour Frame 0 positions (index % 3 == 0) **red** and all other positions **blue**. Mark the start codon (position 0) with a vertical dashed line.

**What to look for:**  
Red bars should be consistently taller than blue bars starting from position 0, visually demonstrating the 3 nt periodicity.

**Output:** `output/analysis2_metagene.png`

---

### Step 3 — Frame phasing quantification and chi-square test (Figure 2B)

**What:**  
Restrict the profile to the CDS window `[0, 99]` nt (positions 0–99 after the AUG). Assign each position to a frame:

| Frame | Positions |
|-------|-----------|
| Frame 0 (in-frame) | 0, 3, 6, 9, … |
| Frame 1 (+1 shift) | 1, 4, 7, 10, … |
| Frame 2 (+2 shift) | 2, 5, 8, 11, … |

Sum read counts in each frame. Compute the percentage each frame contributes to the total CDS reads. Run a **chi-square goodness-of-fit test** against the null hypothesis of a uniform 33.3%:33.3%:33.3% distribution.

Draw a bar chart with three bars (Frame 0/1/2), adding the chi-square statistic and p-value in the title. Annotate each bar with its percentage.

**Why the CDS window only (not the full −50 to +100 range):**  
The −50 to 0 nt region contains the 5′ UTR and initiation complex; ribosome occupancy there is irregular and does not follow the regular reading-frame pattern. Restricting to 0–99 nt isolates genuine elongation-phase triplet phasing.

**Expected result:** Frame 0 ≥ 60–70% of CDS reads; chi-square p < 10⁻¹⁰.

**Output:** `output/analysis2_phasing.png`

---

### Step 4 — FFT power spectrum + autocorrelation (Figure 2C)

**What (FFT):**  
Apply `numpy.fft.rfft` to the mean-centred full meta-gene profile (−50 to +100 nt). Compute the one-sided power spectrum. Convert frequency to period (nt). Identify the peak power in the range `[2.5, 3.5]` nt and compute the fraction of total spectral power (excluding DC component) that falls in this period-3 band.

**What (autocorrelation):**  
Compute the full cross-correlation of the mean-centred profile with itself using `numpy.correlate`, normalise to lag-0 = 1, and plot lags 0–30. A prominent peak at lag 3 (and multiples thereof) confirms that the signal repeats every 3 nt.

Draw a two-panel figure: left panel = FFT power spectrum (x: period in nt, x-range 1–20 nt), right panel = autocorrelation bar chart with Frame 0 lags coloured red.

**Why FFT + autocorrelation in addition to chi-square:**  
- Chi-square tests whether the *aggregate* frame distribution is non-uniform but cannot detect *where* the periodicity starts or how far it extends.  
- FFT decomposes the full meta-gene signal into frequency components, quantifying what fraction of the total signal variance is explained by the period-3 component.  
- Autocorrelation at lag 3 confirms that the signal truly repeats at multiples of 3, providing independent evidence that is not sensitive to the choice of CDS window.

**Expected result:**  
- Period-3 band captures > 30% of total spectral power  
- Autocorrelation at lag 3 > 0.4

**Output:** `output/analysis2_fft.png`

---

### Step 5 — Write numerical summary

**What:**  
Write `output/analysis2_stats.txt` containing:
- Window size used
- Total reads in profile
- Frame 0/1/2 counts and percentages
- Chi-square statistic and p-value
- FFT period-3 power fraction
- Autocorrelation coefficients at lag 1, 2, 3

**Output:** `output/analysis2_stats.txt`

---

## Expected Outputs

| File | Description |
|------|-------------|
| `output/analysis2_metagene.png` | Meta-gene bar chart (−50 to +100 nt). Frame-0 positions (red) should show consistently higher peaks, with the pattern beginning sharply at position 0 (AUG). |
| `output/analysis2_phasing.png` | Frame 0/1/2 bar chart with percentage labels. Frame 0 should be dominant (~60–70%). Chi-square p-value annotated in title. |
| `output/analysis2_fft.png` | Two-panel: FFT power spectrum (clear spike at period = 3 nt) + autocorrelation function (peaks at lags 3, 6, 9). |
| `output/analysis2_stats.txt` | Numerical summary of all statistics for reporting. |

---

## How to Run

```bash
# From the project root
source analysis/.venv/bin/activate

# Ensure mission2 pipeline output exists
ls mission2/work/fivepcounts-filtered-RPF-siLuc.txt

# Run Analysis 2
python analysis/analysis2/analysis2_triplet_phasing.py
```

All output files are written to `analysis/output/` (one level up from `analysis2/`).

---

## Parameters Reference

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Window | −50 to +100 nt | Captures 5′ UTR scanning region + first 100 nt of CDS |
| CDS phasing window | 0 to +99 nt | Avoids UTR noise; covers ~33 codons |
| Chunk size | 500,000 rows | Keeps RAM < 500 MB |
| FFT period mask | 2.5–3.5 nt | Captures the period-3 peak ± 0.5 nt tolerance |
| Autocorrelation range | 0–30 lags | Covers 10 full codon cycles |

---

## Interpretation Guide

| Metric | Strong phasing | Weak phasing |
|--------|---------------|--------------|
| Frame 0 % | ≥ 60% | < 45% |
| Chi-square p | < 10⁻¹⁰ | > 0.01 |
| FFT period-3 power | > 30% of total | < 10% |
| ACF at lag 3 | > 0.4 | < 0.1 |

Strong phasing (all metrics in the left column) indicates a high-quality RPF library where ribosome positions are accurately captured and the data faithfully represents biological translation elongation.
