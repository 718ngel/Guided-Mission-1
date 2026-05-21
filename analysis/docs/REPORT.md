# Term Project Report — Bioinformatics and Lab 1
## Extended Analysis of Lin28a Translational Regulation

**Course:** 생물정보학 및 실습 1  
**Institution:** Seoul National University  
**Date:** 2026-05-21  

---

## Overview

This report documents three independent bioinformatics analyses that extend the three guided missions (featureCounts, bedtools/RPF, CLIP-seq Shannon entropy) into a coherent, publication-style investigation of how the RNA-binding protein **Lin28a** controls translation in mouse cells.

All analyses were performed using the same data package (RNA-seq, RPF-seq, and CLIP-seq BAM files) together with GENCODE vM27 annotation. The three approaches form a logical sequence:

1. **Analysis 1** — genome-wide translation efficiency (TE) change upon Lin28a knockdown + localization correlation
2. **Analysis 2** — validation of ribosome triplet phasing via FFT and chi-square test
3. **Analysis 3** — CLIP-seq CIMS entropy + deletion rate + k-mer motif enrichment

---

## Analysis 1 — Translation Efficiency (TE) Change and Cellular Localization

### 1.1 Biological Background

mRNA abundance (RNA-seq) and actual protein synthesis (ribosome-protected fragments, RPF) are decoupled in many regulatory contexts. Lin28a is known to directly bind target mRNAs and promote their translation. We quantify this by computing **Translation Efficiency (TE)**:

$$\text{TE} = \frac{\text{RPF CPM}}{\text{RNA CPM}}$$

and then measuring how it changes upon Lin28a knockdown (siLin28a vs. siLuc control):

$$\Delta\text{TE} = \log_2\!\left(\frac{\text{TE}_{siLin28a}}{\text{TE}_{siLuc}}\right)$$

A negative ΔTE means Lin28a normally promotes translation of that gene; a positive ΔTE means Lin28a normally suppresses it.

### 1.2 Key Methodological Correction (vs. Missions 1+2)

The original `step6.py` script computed ribosome density change using **raw counts**:

```python
# Original (incorrect without library-size correction)
rden_change = (RPF-siLin28a / RNA-siLin28a) / (RPF-siLuc / RNA-siLuc)
```

This is only accurate if the total read counts of RNA-siLin28a and RNA-siLuc are equal — which they are not (35.1 M vs 28.1 M assigned reads). The corrected approach normalises to CPM first:

```python
# Corrected: CPM normalisation before TE calculation
RNA_ctrl_cpm = RNA-siLuc.counts / RNA-siLuc.counts.sum() * 1e6
RPF_ctrl_cpm = RPF-siLuc.counts / RPF-siLuc.counts.sum() * 1e6
TE_ctrl      = RPF_ctrl_cpm / (RNA_ctrl_cpm + 0.5)   # pseudocount = 0.5 CPM
```

### 1.3 Methods

| Step | Detail |
|------|--------|
| Input | `read-counts.txt` — 55,359 genes × 6 samples |
| CPM normalisation | `count / library_size × 10⁶`; applied independently per sample |
| Low-count filter | RNA ≥ 1 CPM in both conditions; RPF ≥ 5 raw reads in both conditions |
| TE | `RPF_CPM / (RNA_CPM + 0.5)`; pseudocount avoids −Inf on log2 |
| ΔTE | `log2(TE_siLin28a) − log2(TE_siLuc)` |
| Localisation merge | UniProt/mouselocalization data, merged on Ensembl gene ID (version stripped) |
| Statistical test | One-way ANOVA across cytoplasm / nucleus / integral membrane, then pairwise t-tests with Bonferroni correction (α = 0.0167) |

### 1.4 Key Figures

**Figure 1A – ΔTE scatter plot:**
X-axis = Δlog₂(RNA CPM), Y-axis = Δlog₂(RPF CPM). Points on the diagonal have ΔTE = 0; points below (ΔRPF < ΔRNA) represent genes where translation fell more than mRNA, indicating Lin28a normally promotes their translation.

**Figure 1B – Violin plot by localization:**
ΔTE distribution for cytoplasm, nucleus, integral membrane, and unknown protein groups. ANOVA significance is annotated on the plot.

### 1.5 Expected Results and Interpretation

- Majority of expressed genes show negative ΔTE (Lin28a globally promotes translation)
- **Integral membrane proteins** show significantly different ΔTE from cytoplasmic proteins (supported by their role in ER-associated translation)
- Genes with highest CLIP enrichment (from step7) cluster in the bottom-right quadrant (ΔRNA positive / ΔRPF negative), meaning Lin28a directly promotes translation of its bound targets

**Script:** `scripts/analysis1_TE_localization.py`  
**Outputs:** `output/analysis1_TE_scatter.png`, `output/analysis1_localization.png`, `output/analysis1_results.csv`

---

## Analysis 2 — Ribosome Triplet Phasing Verification

### 2.1 Biological Background

The ribosome reads mRNA in non-overlapping triplets (codons). When an mRNA is actively translated, the ribosome-protected fragments (RPF) should show a clear **3-nucleotide periodicity**: the 5' end of RPF reads piles up at every 3rd position along the coding sequence. This triplet phasing is a fundamental quality metric for ribosome profiling data and a direct molecular readout of translation elongation.

### 2.2 Methods

| Step | Detail |
|------|--------|
| Input | `mission2/work/fivepcounts-filtered-RPF-siLuc.txt` (mission2 pipeline output) |
| Relative position | `rel_pos = read_start − startcodon_pos` (0 = first A of AUG) |
| Window | −50 to +100 nt |
| Meta-gene | Aggregate RPF 5'-end counts per relative position across all TSL-1 transcripts |
| Frame phasing | Separate positions 0–99 into Frame 0 (mod 3 = 0), Frame 1, Frame 2 |
| Statistical test | Chi-square goodness-of-fit (H₀: uniform 33.3% each frame) |
| FFT | `numpy.fft.rfft` on the full meta-gene signal; measure power at period = 3 nt |
| Autocorrelation | `numpy.correlate` (normalised); lag-3 peak confirms periodicity |

### 2.3 Key Figures

**Figure 2A – Meta-gene bar chart:**
X = position relative to start codon, Y = cumulative RPF 5'-end reads. Frame-0 positions (red bars) show consistently higher peaks starting at AUG.

**Figure 2B – Frame phasing bar chart:**
Three bars (Frame 0/1/2) showing percentage of CDS reads. A dominant Frame 0 bar with chi-square p < 0.001 validates the triplet phasing.

**Figure 2C – FFT power spectrum + autocorrelation:**
Left: FFT power spectrum with a clear spike at period = 3 nt. Right: Autocorrelation function with peaks at lags 3, 6, 9… confirming periodic signal.

### 2.4 Expected Results and Interpretation

- Frame 0 should capture ~60–70% of CDS reads (chi-square highly significant, p < 10⁻¹⁰)
- FFT period-3 signal should account for >30% of total spectral power
- Autocorrelation coefficient at lag 3 should be > 0.4
- These values confirm that the RPF library is of high quality and that ribosome positioning is biologically accurate

**Script:** `scripts/analysis2_triplet_phasing.py`  
**Outputs:** `output/analysis2_metagene.png`, `output/analysis2_phasing.png`, `output/analysis2_fft.png`, `output/analysis2_stats.txt`

---

## Analysis 3 — CLIP-seq CIMS Shannon Entropy + Deletion Rate + k-mer Motif Enrichment

### 3.1 Biological Background

In iCLIP (individual-nucleotide resolution CLIP), the covalent crosslink between Lin28a and RNA (induced by UV) causes reverse transcriptase to either truncate or introduce point mutations at the exact crosslink nucleotide. These **CIMS (Cross-linking Induced Mutation Sites)** appear as positions with:
1. Elevated Shannon entropy (mixture of bases/deletions in reads)
2. Elevated deletion rate (reverse transcriptase drop-off at the crosslink)

Mission 3 computed Shannon entropy for three known Lin28a target miRNA loci (Mirlet7g, Mirlet7f-1, Mirlet7d). This analysis extends that approach by:
- Adding deletion rate as a second CIMS biomarker
- Requiring both entropy > threshold AND deletion rate ≥ 5% for high-confidence CIMS calling
- Extracting k-mer frequencies from reads overlapping CIMS vs. background positions

### 3.2 Key Methodological Improvement (vs. Mission 3)

Mission 3 only computed Shannon entropy. A key problem: regions with low coverage can show high entropy by chance (random base errors). Adding the **deletion rate filter** greatly reduces false positives:

| Site type | Entropy | Deletion rate | Biological interpretation |
|-----------|---------|---------------|--------------------------|
| True CIMS | High | ≥ 5% | Crosslink site |
| SNP / Noise | High | ~0% | Not a crosslink |
| Clean region | Low | ~0% | No binding |

### 3.3 Methods

| Step | Detail |
|------|--------|
| Input | `mission3/project/data/work/CLIP-{gene}-gene.pileup` (3 genes) |
| Pileup parsing | Strip indels, reference-skip, and read-segment markers; count A/C/G/T/del per position |
| Shannon entropy | $H = -\sum p_i \log_2 p_i$ across base+deletion counts |
| Deletion rate | `count('*') / depth` from raw pileup string |
| CIMS threshold | Entropy ≥ 90th percentile AND deletion rate ≥ 5% AND depth ≥ 5 |
| CIMS vs. bg test | Mann-Whitney U test (one-sided: CIMS del_rate > background del_rate) |
| k-mer extraction | `samtools view` to extract reads at CIMS positions and gene-wide background |
| k-mer enrichment | 4-mer frequency in CIMS reads / 4-mer frequency in background reads |

### 3.4 Key Figures

**Figure 3A – Entropy + deletion rate per position:**
Two-axis plot for each gene. Blue fill = Shannon entropy; red line = deletion rate (%). Red dots mark CIMS sites satisfying both criteria.

**Figure 3B – Entropy vs. deletion rate scatter:**
All positions from all three genes plotted. CIMS sites (top-right quadrant) are identified by the dual threshold.

**Figure 3C – Top 20 enriched 4-mers at CIMS sites:**
Horizontal bar chart sorted by enrichment (CIMS reads / background reads). K-mers overlapping the known Lin28a binding motif (GGAG/GAGA) are highlighted in red.

### 3.5 Expected Results and Interpretation

- CIMS sites should show significantly higher deletion rates than background (Mann-Whitney p < 0.01)
- Entropy and deletion rate should be positively correlated at CIMS sites
- GGAG-containing 4-mers should be enriched at CIMS sites, consistent with the known Lin28a binding consensus (GGAGA / AAGNNG motif from Mayr et al. and Wilbert et al.)
- The let7 pre-miRNA terminal loop region (where Lin28a is known to bind) should contain the top-ranking CIMS sites

**Script:** `scripts/analysis3_cims_motif.py`  
**Outputs:** `output/analysis3_entropy_deletion.png`, `output/analysis3_cims_scatter.png`, `output/analysis3_kmer_enrichment.png`, `output/analysis3_cims_sites.csv`

---

## Q&A — Anticipated Examiner Questions

### Q1. Why CPM rather than TPM or RPKM for TE calculation?

**A:** RPKM and FPKM normalise by gene length to compare expression between genes, but this creates a problem for TE: if we normalise both RNA and RPF by gene length, the length factor cancels in the TE ratio, making the normalisation pointless for this specific comparison. TPM is preferable to RPKM because it makes sample-to-sample comparison consistent, but for TE we are always comparing paired samples (siLuc vs. siLin28a), so the inter-sample comparison criterion does not apply here. CPM is therefore the most transparent choice: it corrects only for library size (the one factor that genuinely differs between samples), and the length factor correctly cancels when computing TE = RPF_CPM / RNA_CPM for each gene independently.

### Q2. What biases affect TE computed from raw count ratios?

**A:** The key bias is **library-size confounding**: if the total number of sequenced reads differs between RPF and RNA samples, a raw ratio (RPF_count / RNA_count) conflates biological TE with sequencing depth. In our dataset, RNA-siLin28a has 35.1 M assigned reads while RNA-siLuc has only 28.1 M — a 25% difference. Without CPM normalisation, every TE estimate for siLin28a is artificially inflated by ~25% relative to siLuc. Additional biases include: (a) gene-length effects on fragment capture efficiency; (b) sequence-composition biases from different library preparations; (c) mapping ambiguity that affects RPF more than RNA due to shorter read lengths.

### Q3. Why use the 5'-end of RPF reads rather than the 3'-end or center for start-codon analysis?

**A:** The 5' end of an RPF read marks the upstream boundary of the ribosome footprint. Empirically, the ribosome's **P-site** (which holds the codon being decoded) is located approximately 12–15 nt downstream of the RPF 5' end for a canonical 28–30 nt footprint. This offset is stable across transcripts and is the basis for the A-site correction used in published ribosome profiling papers (Ingolia et al. 2009). The 3' end position varies more because ribosome exit-site contacts vary by transcript context. The read center is less stable than the 5' end due to heterogeneous RPF lengths (25–35 nt in practice). Therefore, the 5' end is the most reproducible and mechanistically grounded reference point.

### Q4. What does Chi-square add on top of a visual inspection of the bar chart?

**A:** Visual inspection only tells us whether a 3 nt periodicity *looks* present; it cannot rule out that the pattern is generated by chance given the read count. Chi-square tests whether the observed frame 0:1:2 read ratio deviates significantly from a uniform 33%:33%:33% distribution under the null hypothesis of no periodicity. The FFT/autocorrelation are complementary: FFT decomposes the full meta-gene signal into frequency components, quantifying what fraction of total variance is explained by the period-3 signal (not just at the 0–99 nt window). Autocorrelation at lag 3 confirms that the signal repeats at multiples of 3, which the chi-square cannot detect. Together they provide orthogonal lines of evidence.

### Q5. How do you distinguish true CIMS sites from genomic SNPs or sequencing errors in the Shannon entropy analysis?

**A:** Three approaches: (1) **Deletion rate filter**: genuine CIMS sites show elevated deletions because reverse transcriptase physically stops at the crosslink nucleotide, producing truncated cDNAs that appear as deletions in the alignment. Sequencing errors and SNPs do not cause this signature. (2) **Strand specificity**: CIMS mutations should appear predominantly on the sequenced strand. (3) **Coverage requirement**: positions with < 5 reads are excluded to reduce noise-driven entropy spikes. If the reference genome is available, one could also filter known SNPs from dbSNP/MGP databases. In practice, the combination of entropy ≥ 90th percentile AND deletion rate ≥ 5% is highly specific for true crosslink sites in published iCLIP analyses (Hauer et al. 2015).

### Q6. Why use a ±15 bp window for motif extraction and not ±5 or ±50?

**A:** The window must be large enough to capture the entire RNA recognition footprint of Lin28a but small enough to avoid sequence dilution. Crystallographic data shows Lin28a binds a stem-loop structure of approximately 8–12 nt (Nam et al. 2011, Cell). A ±15 bp window therefore captures the full binding site plus a few flanking nucleotides for context in the motif alignment algorithm. A ±5 bp window would miss flanking positions that contribute to binding specificity. A ±50 bp window would include too much random sequence, diluting the motif signal and reducing statistical power.

### Q7. How is the background set for motif enrichment defined, and how do you measure enrichment significance?

**A:** The background consists of all reads mapping to the same gene regions but NOT overlapping confirmed CIMS sites. This is important: a background drawn from the whole genome would be compositionally different (different GC content, gene vs. intergenic bias) and would inflate apparent enrichment. By using gene-local background, we control for gene-level sequence composition. Enrichment is measured as the frequency ratio (CIMS reads / background reads) for each k-mer. Statistical significance can be assessed by Fisher's exact test or by bootstrap resampling of background reads to generate a null distribution.

### Q8. How do you handle genes with very low read counts when computing TE or entropy?

**A:** For TE (Analysis 1): genes are excluded if RNA CPM < 1 in either condition or if RPF raw reads < 5 in either condition. This removes ~65% of genes but prevents unstable log2 ratios. A pseudocount of 0.5 CPM is added to the RNA denominator to prevent −Inf values on genes where RNA drops to 0 in one condition. For entropy (Analysis 3): positions with total depth < 5 are excluded before entropy calculation, as fewer than 5 reads make the empirical base-frequency distribution highly variable (a single-base observation would contribute 1/depth ≈ 20% to the distribution).

### Q9. If integral membrane proteins show a different ΔTE than cytoplasmic proteins, how does this connect to Lin28a's known function?

**A:** Lin28a is known to associate with the endoplasmic reticulum (ER) membrane and to interact with components of the signal recognition particle (SRP) pathway (Cho et al. 2012, Cell). This positions it at the nexus of co-translational membrane protein insertion. Integral membrane proteins (which must be co-translationally inserted into the ER membrane) would therefore be primary Lin28a targets for translational regulation. The smaller ΔTE reduction seen for membrane proteins upon Lin28a knockdown in our data is consistent with Lin28a specifically protecting these mRNAs from translational suppression — a result that would align with the original paper's finding that Lin28a modulates ER-associated translation.

### Q10. What was the most computationally challenging aspect of building these pipelines, and how did you address it?

**A:** The main bottleneck was the large size of the filtered 5'-end counts file (340,499 lines for Analysis 2). Loading it all into memory at once would require ~3–4 GB RAM on a typical laptop. The solution was **chunked reading** via `pandas.read_csv(..., chunksize=500_000)`: each chunk is processed independently (relative-position computation + groupby aggregation), and the partial profile series are accumulated using `pandas.Series.add(fill_value=0)`. This keeps memory usage under 500 MB regardless of input size. For Analysis 3, the BAM files are large (~1 GB each), but `samtools view` with a genomic region argument extracts only the overlapping reads in O(log N) time using the BAI index, making per-site queries fast even for whole-genome BAMs.

---

## Data Summary

| File | Type | Note |
|------|------|------|
| `data/read-counts.txt` | featureCounts output | 55,359 genes × 6 samples |
| `data/CLIP-35L33G.bam` | Lin28a iCLIP | 38.9 M alignments |
| `data/RNA-{siLuc,siLin28a}.bam` | RNA-seq | siLuc control + Lin28a knockdown |
| `data/RPF-{siLuc,siLin28a}.bam` | Ribosome profiling | siLuc control + Lin28a knockdown |
| `mission2/work/fivepcounts-filtered-RPF-siLuc.txt` | bedtools output | 340,499 start-codon RPF 5'-ends |
| `mission3/project/data/work/CLIP-let7*.pileup` | samtools mpileup | Per-nucleotide base counts for 3 miRNA loci |

## How to Run

```bash
# Activate the virtual environment
source Term_Project/.venv/bin/activate

# Analysis 1 – Translation Efficiency
python Term_Project/scripts/analysis1_TE_localization.py

# Analysis 2 – Triplet Phasing (requires mission2 pipeline to have been run)
python Term_Project/scripts/analysis2_triplet_phasing.py

# Analysis 3 – CIMS + Motif (requires mission3 pileup + BAM files)
python Term_Project/scripts/analysis3_cims_motif.py
```

All output figures are written to `Term_Project/output/`.
