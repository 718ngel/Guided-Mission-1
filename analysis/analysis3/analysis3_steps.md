# Analysis 3 — CLIP-seq CIMS Shannon Entropy and k-mer Motif Enrichment at Lin28a Binding Sites

**Course:** 생물정보학 및 실습 1  
**Institution:** Seoul National University  
**Reference:** `analysis/CLAUDE.md` → Section "Analysis 3 (Extended)"  
**Script:** `analysis/analysis3/analysis3_cims_motif.py`

---

## Overview

This analysis characterises the **molecular binding footprint of Lin28a on its target pre-miRNA transcripts** by combining three orthogonal signals extracted from CLIP-seq pileup data:

1. **Shannon entropy** — a position-wise measure of nucleotide diversity that peaks where reverse transcription errors cluster (i.e., at UV cross-linking sites).
2. **Deletion rate** — the fraction of reads carrying a deletion (`*` in samtools pileup format) at each position; the primary biomarker of CIMS (Cross-linking-Induced Mutation Sites).
3. **k-mer motif enrichment** — comparison of 4-mer nucleotide frequencies in reads overlapping confirmed CIMS sites versus background reads, to recover the RNA sequence motif bound by Lin28a.

The three target loci are the mouse precursor let-7 family members **let7g** (chr9), **let7f1** (chr13), and **let7d** (chr13) — well-characterised direct binding targets of Lin28a.

The analysis extends the Mission 3 single-nucleotide Shannon entropy pipeline by adding the deletion-rate confirmatory filter, a Mann-Whitney U statistical test, BAM-level read extraction, and comparative k-mer frequency analysis.

---

## Biological Background

### Lin28a as a regulator of let-7 biogenesis

Lin28a (Lin28 homolog A) is a highly conserved RNA-binding protein that functions as a master regulator of the let-7 microRNA family. In stem cells and early development, Lin28a binds to the terminal loop of let-7 precursor miRNAs (pre-let-7) and inhibits their processing by two mechanisms:

- **Microprocessor blockade:** Lin28a interacts with the Drosha–DGCR8 complex and sterically prevents cleavage of primary miRNA (pri-let-7) to precursor miRNA (pre-let-7).
- **TUT4/TUT7 recruitment:** Lin28a recruits the poly(U) polymerases TUT4 and TUT7 to oligouridylate the 3′ tail of pre-let-7, targeting it for degradation by the DIS3L2 exoribonuclease.

Both mechanisms converge on the same effect: Lin28a expression suppresses mature let-7 accumulation. Because let-7 miRNAs target oncogenes (KRAS, MYC, HMGA2), Lin28a overexpression is associated with embryonic stemness and tumourigenesis.

The primary Lin28a RNA-binding motif is a **GGAGA** or **GGAG** tetranucleotide within the terminal loop of pre-let-7. The Lin28a cold shock domain (CSD) contacts this motif directly, while the CCHC zinc-knuckles (ZKD) bind the adjacent single-stranded RNA. This makes the GGAG-containing hairpin loop of pre-let-7 both a biologically essential sequence and a strong predicted enrichment signal in CIMS k-mer analysis.

### CLIP-seq and UV cross-linking

CLIP-seq (Cross-linking and ImmunoPrecipitation followed by high-throughput Sequencing) identifies the in vivo RNA-binding sites of an RNA-binding protein (RBP) at nucleotide resolution:

1. Living cells are irradiated with 254 nm UV light. This induces covalent cross-links (specifically pyrimidine-like adducts) between RNA bases and nearby amino acid side chains when the protein is in direct contact with the RNA.
2. The cross-linked RBP–RNA complexes are immunoprecipitated, the RNA is partially fragmented, and the remaining RNA fragment is reverse-transcribed into cDNA.
3. During reverse transcription, the covalent adduct at the cross-linking site causes the reverse transcriptase to either **misincorporate** a base (creating a substitution) or **fall off** the template (creating a deletion or truncation). These are detected as pileup anomalies in the final alignment.

### CIMS (Cross-linking-Induced Mutation Sites)

CIMS are positions in the genome where a statistically elevated frequency of deletions or mutations in CLIP-seq reads marks the precise amino acid–nucleotide contact point. The original CIMS methodology (Zhang and Darnell, 2011) focused on deletions because:

- Deletions are rare in standard RNA-seq libraries (< 0.1% per position), making them a highly specific signal.
- UV cross-linking dramatically elevates the deletion rate at the contact nucleotide (typically 5–20% of reads).
- Other sources of artefact (SNPs, sequencing errors) contribute uniformly across positions, not at specific CLIP-enriched sites.

### Shannon entropy as a complementary CIMS detector

Shannon entropy measures the informational diversity of base calls at a single genomic position. Given counts of A, C, G, T reads at position i:

$$H_i = -\sum_{b \in \{A,C,G,T\}} \frac{n_b}{N} \log_2 \frac{n_b}{N}$$

where $N = \sum n_b$ is the total coverage depth.

- A perfectly conserved position (all reads report the same base) yields $H = 0$ bits.
- A position where the four bases appear equally often yields $H = 2$ bits (maximum for a 4-symbol alphabet).
- Cross-linking-induced misincorporation events increase the proportion of non-reference bases, elevating $H$ above the background level (typically ~0.1–0.3 bits for a position with only sequencing-error-level noise).

Using entropy *alone* is insufficient because genomic SNPs, RNA editing (A-to-I), and repetitive elements can also elevate entropy. Requiring a **co-occurring deletion rate ≥ 5%** provides an orthogonal biochemical confirmation: only UV cross-linking generates the distinctive deletion signature.

The dual criterion (entropy ≥ 90th percentile **AND** deletion rate ≥ 5%) is therefore more specific than either signal alone.

### Target loci

| Gene | Chromosome | Start | End | Strand | Biology |
|------|-----------|-------|-----|--------|---------|
| let7g | chr9 | 106,056,039 | 106,056,126 | + | pre-miR-let-7g; GGAG motif in terminal loop |
| let7f1 | chr13 | 48,691,305 | 48,691,393 | − | pre-miR-let-7f-1; paralogous to let7g |
| let7d | chr13 | 48,689,488 | 48,689,590 | − | pre-miR-let-7d; clustered with let7f1 |

All three encode stem-loop precursor miRNAs that are direct Lin28a substrates in mouse embryonic stem cells.

---

## Input Data

### Pileup files (primary analysis)

| File | Description |
|------|-------------|
| `mission3/project/data/work/CLIP-let7g-gene.pileup` | samtools mpileup output at the let7g locus |
| `mission3/project/data/work/CLIP-let7f1-gene.pileup` | samtools mpileup output at the let7f1 locus |
| `mission3/project/data/work/CLIP-let7d-gene.pileup` | samtools mpileup output at the let7d locus |

**Pileup column format (6 columns, tab-separated, no header):**

| Column | Field | Content |
|--------|-------|---------|
| 1 | `chrom` | Chromosome name |
| 2 | `pos` | 1-based genomic position |
| 3 | `_ref` | Reference base |
| 4 | `depth` | Total read depth at this position |
| 5 | `basereads` | Base call string (`.` = reference match on `+` strand, `,` = reference match on `−` strand, `A/C/G/T` = mismatch, `*` = deletion, `^` = read start, `$` = read end, `+nNNN`/`-nNNN` = indel) |
| 6 | `quals` | Base quality string (ASCII-encoded Phred scores) |

> The `basereads` string must be cleaned before counting: read-start markers `^<ASCII quality>`, read-end `$`, indel notation `+nNNN`/`-nNNN`, and mapping annotations are stripped with a regular expression before base tallying. Deletion tokens (`*`) are counted *before* cleaning to compute the deletion rate.

### BAM files (k-mer analysis)

| File | Description |
|------|-------------|
| `mission3/project/data/work/CLIP-let7g.bam` | Coordinate-sorted BAM for let7g |
| `mission3/project/data/work/CLIP-let7f1.bam` | Coordinate-sorted BAM for let7f1 |
| `mission3/project/data/work/CLIP-let7d.bam` | Coordinate-sorted BAM for let7d |

Corresponding `.bam.bai` index files must also be present for `samtools view` region queries.

> **Prerequisite:** Mission 3 CLIP-seq alignment pipeline must be complete (BAM sorted and indexed, pileup files generated with `samtools mpileup`).

---

## Step-by-Step Pipeline

### Step 1 — Parse pileup files and compute per-position metrics

**What:**  
For each of the three genes, read the corresponding pileup file into a `pandas.DataFrame`. Apply the following transformations per row (i.e., per genomic position):

1. **Strip annotation tokens** from the `basereads` string using the regex pattern  
   ```
   [<>$*#^]|\+\d+[ACGTNacgtn]+|-\d+[ACGTNacgtn]+
   ```
   to obtain `clean_reads` (contains only base-call characters).

2. **Count bases** — build a `Counter` of uppercase `{A, C, G, T, ., ,}` from `clean_reads`. Dots (`.`) and commas (`,`) denote reference-matching bases on the `+` and `−` strands, respectively.

3. **Shannon entropy** — compute $H_i$ using the formula above. Positions with `depth < 5` are discarded to prevent noise-dominated entropy estimates.

4. **Deletion rate** — count raw `*` characters in the original (uncleaned) `basereads` string and divide by `depth`:
   $$\text{del\_rate}_i = \frac{\#(*\text{ in basereads})}{depth_i}$$
   This must be computed *before* string cleaning because `*` tokens are removed by the regex.

5. **Mismatch rate** — the fraction of `clean_reads` positions that carry an explicit base letter (`A/C/G/T`, i.e., a mismatch vs. reference) out of total depth.

**Why depth ≥ 5 filter:**  
At very low coverage (1–4 reads), a single sequencing error or one deletion event inflates both entropy and deletion rate to unreliable extremes. The threshold of 5 reads provides a minimal denominator such that a single deletion represents a maximum deletion rate of 20%, preventing false positives caused by stochastic sampling of very sparse data.

**Output:** A per-gene `DataFrame` with columns `[chrom, pos, depth, entropy, del_rate, mis_rate, gene, strand]`, concatenated into `combined`.

---

### Step 2 — CIMS site calling with dual criteria

**What:**  
Compute the 90th-percentile entropy threshold across all positions in `combined`:

$$\theta_{H} = \text{quantile}_{0.90}(\{H_i : i \in \text{all positions}\})$$

Classify each position into one of three categories:

| Category | Criterion | Interpretation |
|----------|-----------|----------------|
| Background | $H_i < \theta_H$ | Low base diversity — reference-matching region |
| High-entropy only | $H_i \geq \theta_H$ | Elevated diversity but no deletion enrichment — possible SNP, RNA editing, or mapping ambiguity |
| **Confirmed CIMS** | $H_i \geq \theta_H$ **AND** $\text{del\_rate}_i \geq 0.05$ | High diversity + deletion enrichment — UV cross-linking signature |

**Why 90th percentile (not absolute threshold):**  
The absolute entropy scale depends on the local read depth and nucleotide composition of the pre-miRNA hairpin. Using a data-adaptive quantile threshold makes the calling criterion robust to systematic differences in library depth between let7g, let7f1, and let7d. Only the top 10% of entropy values are considered as candidates, of which only a subset will also pass the deletion-rate filter.

**Why del_rate ≥ 5%:**  
The human CLIP-CIMS literature (Zhang et al., 2011; Sugimoto et al., 2012) reports that authentic cross-linking sites carry deletion rates of 5–30% depending on RBP–amino acid chemistry. A 5% threshold provides sensitivity to detect most genuine sites while excluding background noise (typically < 1–2% deletion rate).

---

### Step 3 — Mann-Whitney U test for deletion rate enrichment

**What:**  
Test whether confirmed CIMS sites have significantly higher deletion rates than background positions using a one-tailed Mann-Whitney U test:

$$H_0: \text{del\_rate}_{\text{CIMS}} \leq \text{del\_rate}_{\text{background}}$$
$$H_1: \text{del\_rate}_{\text{CIMS}} > \text{del\_rate}_{\text{background}}$$

```python
u_stat, p_mwu = stats.mannwhitneyu(
    cims_df["del_rate"],
    bg_df["del_rate"],
    alternative="greater"
)
```

**Why Mann-Whitney U (not t-test):**  
Deletion rate distributions are strongly right-skewed and non-normal (most positions have near-zero deletion rates, with a long tail at CIMS sites). The Mann-Whitney U is a non-parametric rank-based test that does not assume normality and is robust to the large imbalance between the small CIMS group and large background group.

**Expected result:**  
U statistic significantly > chance level; p < 0.05 (likely p << 0.001 given the biological effect size of UV cross-linking).

---

### Step 4 — Per-gene entropy and deletion-rate profile plot (Figure 3A)

**What:**  
For each gene, draw a dual-axis line plot (one panel per gene, stacked vertically):

- **Left y-axis (blue fill):** Shannon entropy per position as an area fill (`fill_between`).
- **Right y-axis (red line):** Deletion rate (%) as a line overlay.
- **Red scatter dots:** Positions classified as confirmed CIMS sites are overlaid on the entropy fill.
- **x-axis:** Genomic coordinate in the form `chr:start–end`.

**What to look for:**  
The entropy and deletion-rate signals should be **correlated at CIMS sites** — positions where entropy peaks should coincide with elevated deletion rate. The co-occurrence of both signals distinguishes genuine Lin28a contact sites from background noise.

Within each pre-miRNA, the CIMS sites should localise to the terminal loop region (the hairpin apex), which is the known single-stranded region that accommodates Lin28a CSD binding.

**Output:** `analysis/output/analysis3_entropy_deletion.png`

---

### Step 5 — Entropy vs. deletion-rate scatter plot (Figure 3B)

**What:**  
Draw a single scatter plot pooling all three genes (coloured by gene: let7g = blue, let7f1 = orange, let7d = green):

- **x-axis:** Shannon entropy (bits)
- **y-axis:** Deletion rate (%)
- **Horizontal dashed line:** deletion rate = 5% (CIMS threshold)
- **Vertical dotted line:** 90th-percentile entropy threshold

Points in the **upper-right quadrant** (above both threshold lines) are the confirmed CIMS sites.

**Why this plot:**  
The scatter visualises the relationship between the two complementary signals and makes the dual-threshold criterion transparent. It also allows visual assessment of whether the two signals are genuinely correlated (as expected from biology) or whether one threshold is capturing sites that the other misses (suggesting one signal is noisy in this dataset).

**Output:** `analysis/output/analysis3_cims_scatter.png`

---

### Step 6 — 4-mer k-mer enrichment from BAM reads (Figure 3C)

**What:**  
To identify the RNA sequence motif that Lin28a contacts, compare the nucleotide composition of reads at CIMS sites vs. all reads from the gene locus:

**6a. Extract read sequences from BAM files:**  
For each confirmed CIMS site, run:

```bash
samtools view CLIP-<gene>.bam <chrom>:<pos>-<pos>
```

Parse field 10 (the SEQ column) from each alignment record. These are the **CIMS reads** — reads that directly overlap a confirmed cross-linking site.

For background, also extract all reads from the full gene region (`chr:gstart-gend`). These are the **background reads** — the full population of CLIP reads at each locus.

**6b. Count all 4-mers:**  
For each read sequence, slide a 4-nucleotide window along the read and count every 4-mer that contains only `{A, C, G, T}` (no ambiguous `N`). Accumulate into `Counter` objects separately for CIMS reads and background reads.

**6c. Compute enrichment ratio:**  
Normalise each 4-mer count to frequency (count ÷ total k-mer count):

$$\text{enrichment}(km) = \frac{f_{\text{CIMS}}(km)}{f_{\text{background}}(km)}$$

Discard k-mers with zero background frequency (appear only in CIMS reads — likely artefacts from very small CIMS read sets).

Sort by enrichment (descending) and take the top 20.

**6d. Highlight known Lin28a motifs:**  
The following 4-mers overlap the canonical Lin28a binding motif (GGAG, GAGA, AGGA, GAAG):

| Highlight k-mer | Overlap with motif |
|----------------|--------------------|
| GGAG | Direct match to Lin28a CSD recognition sequence |
| GAGA | Sub-sequence of GGAGA motif |
| AGGA | Reverse complement sub-sequence |
| GAAG | Adjacent context in let-7 terminal loop |

Bars corresponding to these k-mers are coloured **red** in the plot; all others are coloured **blue**.

**Why 4-mers and not 5-mers or 6-mers:**  
The Lin28a CSD recognises a GGAG tetranucleotide (4 nt) directly. Using 4-mers exactly matches the length of the primary motif, maximising sensitivity. Longer k-mers would fragment the signal across too many variants of the same core motif; shorter 3-mers would not resolve the specific arrangement of guanines.

**Why `samtools view` reads instead of pileup sequences:**  
The pileup format reports per-position base calls after alignment, which provides entropy and deletion information but loses information about contiguous sequence context within a read. To recover the local sequence environment of the cross-linking site, the actual read sequences must be retrieved from the BAM file.

**Expected result:**  
GGAG-containing 4-mers should appear in the top-enriched set if Lin28a binding is authentic. Enrichment > 1.5× for GGAG-overlapping k-mers provides positive control validation that the CIMS sites are genuinely Lin28a binding sites.

**Output:** `analysis/output/analysis3_kmer_enrichment.png`

---

### Step 7 — Save CIMS site table

**What:**  
Export all confirmed CIMS sites to a CSV file with columns:

| Column | Description |
|--------|-------------|
| `gene` | let7g / let7f1 / let7d |
| `chrom` | Chromosome |
| `pos` | 1-based genomic position |
| `depth` | Read depth at CIMS site |
| `entropy` | Shannon entropy (bits) |
| `del_rate` | Deletion rate (0–1) |
| `mis_rate` | Mismatch rate (0–1) |

Rows sorted by gene (alphabetically) then by entropy (descending) so that the highest-confidence sites appear first within each gene.

**Output:** `analysis/output/analysis3_cims_sites.csv`

---

## Expected Outputs

| File | Description |
|------|-------------|
| `analysis/output/analysis3_entropy_deletion.png` | Three stacked panels (one per gene). Each panel shows entropy fill (blue) and deletion rate line (red) along the genomic coordinate. Confirmed CIMS sites are overlaid as red scatter points. CIMS positions should cluster at the terminal loop region (~20–30 nt region at the hairpin apex). |
| `analysis/output/analysis3_cims_scatter.png` | Scatter plot of entropy (x) vs. deletion rate (y) for all positions, coloured by gene. Points in the upper-right quadrant above both threshold lines are confirmed CIMS sites. The CIMS cluster should be visually distinct from the main cloud of background points. |
| `analysis/output/analysis3_kmer_enrichment.png` | Horizontal bar chart of the top 20 enriched 4-mers (CIMS / background frequency ratio). GGAG-overlapping k-mers should be among the most enriched and are highlighted in red. |
| `analysis/output/analysis3_cims_sites.csv` | Table of confirmed CIMS sites. Expected: O(1–10) sites per gene given the small (~80–100 bp) target loci. Each site should have del_rate ≥ 0.05 and entropy in the top 10% of all positions. |

---

## How to Run

```bash
# Navigate to the project root
cd /path/to/workspace   # the directory containing analysis/ and mission3/

# Activate the virtual environment
source analysis/.venv/bin/activate    # or conda activate binfo1

# Confirm input data exists
ls mission3/project/data/work/CLIP-let7g-gene.pileup
ls mission3/project/data/work/CLIP-let7g.bam

# Confirm samtools is available (required for BAM k-mer extraction)
samtools --version

# Run Analysis 3
python analysis/analysis3/analysis3_cims_motif.py
```

All four output files are written to `analysis/output/` (one directory above `analysis3/`). The output directory is created automatically if it does not exist.

If BAM files are absent or `samtools` is not in `PATH`, the k-mer enrichment step is skipped gracefully — the three pileup-based outputs (entropy/deletion plot, scatter, CSV) will still be generated.

---

## Parameters Reference

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| `ENTROPY_QUANTILE` | 0.90 | Top 10% of entropy values selected as CIMS candidates; data-adaptive threshold robust to between-locus depth variation |
| `DEL_RATE_MIN` | 0.05 (5%) | Conservative minimum deletion rate consistent with authentic UV cross-linking; avoids SNP and sequencing-error false positives |
| `MIN_DEPTH` | 5 reads | Minimum coverage denominator to prevent stochastic inflation of entropy and deletion rate in sparse regions |
| `K` | 4 (4-mer) | Matches the length of the primary Lin28a CSD recognition motif GGAG; balances specificity with sufficient k-mer count statistics |
| Entropy window | All positions with depth ≥ 5 | No position windowing; all covered positions in the gene pileup are assessed |
| Highlight motifs | GGAG, GAGA, AGGA, GAAG | Canonical Lin28a binding sequences from structural studies (Nam et al., 2011; Loughlin et al., 2012) |

---

## Interpretation Guide

### CIMS site quality criteria

| Metric | Strong CIMS evidence | Weak / inconclusive |
|--------|---------------------|---------------------|
| Deletion rate at CIMS | ≥ 10% | 5–7% (marginal) |
| Entropy at CIMS (bits) | ≥ 1.0 | 0.5–0.8 (borderline) |
| Mann-Whitney U p-value | < 0.001 | > 0.05 (non-significant) |
| GGAG k-mer enrichment | ≥ 1.5× | < 1.2× (no motif recovery) |
| Spatial clustering | CIMS sites within terminal loop | Scattered across full locus |

### Interpreting the entropy-vs-deletion scatter plot

- **Main cloud (low entropy, low deletion rate):** Reference-matching positions — the ribosome or protein is not contacting these bases, or contact does not cause RT errors.
- **Upper-left (high entropy, low deletion rate):** Possible genomic SNPs, RNA editing events, or regions with repetitive mapping. These are *not* confirmed CIMS.
- **Lower-right (low entropy, high deletion rate):** Rare; could reflect systematic deletion artefacts at homopolymer stretches.
- **Upper-right (high entropy, high deletion rate):** Confirmed CIMS — the biologically meaningful cluster where UV cross-linking-induced RT errors concentrate.

### Connecting k-mer results to Lin28a structural biology

The Lin28a cold shock domain (CSD) makes base-specific contacts with the GGAG tetraloop of pre-let-7. If the CIMS sites identified here genuinely mark Lin28a contact points, then:

1. The CIMS positions should map within the terminal loop (≈ 20–30 nt from the hairpin apex of each pre-miRNA).
2. The top-enriched k-mers at CIMS sites should overlap GGAG, confirming that the reads covering these positions are enriched for the known binding motif.
3. The enrichment should be gene-consistent — the same motif family enriched at let7g CIMS sites should also appear at let7f1 and let7d sites.

If these three conditions are met, the analysis provides independent computational evidence that recapitulates known Lin28a binding biochemistry purely from sequencing data, without prior knowledge of the binding motif.

### Limitations and caveats

1. **Single cross-linking condition:** These BAM files represent a single CLIP experiment. Without a no-UV negative control, the baseline deletion rate cannot be calibrated to zero UV background. The 5% threshold is a literature-informed heuristic, not experiment-specific.

2. **Small target loci:** Each pre-miRNA is only ~80–100 bp. The small number of positions limits the statistical power of the Mann-Whitney test and reduces k-mer count statistics (low read coverage → noisy k-mer frequencies). Results should be interpreted with this constraint in mind.

3. **Pileup-to-BAM coordinate concordance:** samtools pileup positions are 1-based; samtools view region coordinates are also 1-based. The script uses the `pos` column from the pileup directly as the region coordinate in `samtools view`, which is correct. Verify this assumption if results appear positionally shifted.

4. **Strand directionality:** let7f1 and let7d are on the `−` strand. Pileup base calls are reported relative to the sequenced read orientation (not the reference strand), so `.`/`,` correctly reflect reference-matching reads on either strand. k-mer enrichment from BAM reads should ideally be reverse-complemented for `−`-strand genes to recover motifs in the RNA (not cDNA) context. The current implementation does not perform this step — enriched k-mers from `−`-strand genes represent the cDNA sequence and must be reverse-complemented to compare with the known Lin28a motif (GGAG on RNA = CTCC on cDNA-sense).

---

## Connection to Other Analyses

| Analysis | Connection to Analysis 3 |
|----------|--------------------------|
| Analysis 1 (TE) | Genes with confirmed Lin28a CIMS sites (pre-let-7 loci) are expected to show indirect TE changes: Lin28a suppresses let-7, which in turn de-represses let-7 target mRNAs. CIMS validation strengthens the mechanistic chain from Lin28a binding → let-7 suppression → TE change. |
| Analysis 2 (Triplet Phasing) | The ribosome-phased RPF data validates that the RNA-seq and RPF libraries used in Analysis 1 are high quality. Similarly, Analysis 3's CIMS validation confirms that the CLIP library has authentic cross-linking chemistry, providing orthogonal quality evidence for the overall dataset. |
| Mission 3 (UCSC / pileup) | Analysis 3 directly extends Mission 3's per-nucleotide entropy calculation by adding deletion-rate analysis, statistical testing, and motif recovery — the three components that elevate the Mission 3 exercise to publication-quality analysis. |
