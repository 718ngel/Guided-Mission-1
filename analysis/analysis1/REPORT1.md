# REPORT1 — Analysis 1: Translation Efficiency (TE) Change and Cellular Localization

**Course:** 생물정보학 및 실습 1  
**Institution:** Seoul National University  
**Date:** 2026-05-21  
**Author:** Angela Kim  

---

## 1. Task List (Pre-execution Checklist)

- [x] Pull latest code from remote: `git pull origin main`
- [ ] Set up `analysis/analysis1/` directory scaffold (scripts/, output/, data/)
- [ ] Copy `read-counts.txt` into `analysis/analysis1/data/`
- [ ] Copy `analysis1_TE_localization.py` into `analysis/analysis1/scripts/`
- [ ] Write and commit this REPORT1.md **before** running any code
- [ ] Run `python analysis/analysis1/scripts/analysis1_TE_localization.py`
- [ ] Verify output: `analysis1_TE_scatter.png`, `analysis1_localization.png`, `analysis1_results.csv`
- [ ] Write `result_analysis1.md` with actual execution results

---

## 2. Biological Objective

### What Translation Efficiency (TE) measures

Translation efficiency quantifies the decoupling between mRNA abundance and protein synthesis. For a given gene, the ribosome density per unit of mRNA — measured by ribosome-protected fragments (RPF-seq) relative to total mRNA (RNA-seq) — gives a genome-wide readout of how efficiently each transcript is being translated at a given moment. Genes with high TE produce many protein molecules per mRNA copy; genes with low TE accumulate mRNA that is not productively translated.

$$\text{TE} = \frac{\text{RPF CPM}}{\text{RNA CPM} + 0.5}$$

The pseudocount of 0.5 CPM is added to the RNA denominator to prevent division by zero and to stabilise log2 ratios for genes expressed at very low levels.

### Why Lin28a knockdown is the perturbation of interest

Lin28a is an RNA-binding protein (RBP) expressed during embryogenesis and in stem cells that directly promotes the translation of specific target mRNAs, including IGF-2 and Hmga2, while also blocking maturation of the let-7 microRNA family. When Lin28a is knocked down using siRNA (siLin28a, compared to a luciferase-targeting siLuc control), we lose the translational boost it normally provides to its target mRNAs. This knockdown experiment therefore reveals — at genome-wide scale — which genes depend on Lin28a for efficient translation.

### What a negative ΔTE means biologically

$$\Delta\text{TE} = \log_2(\text{TE}_{siLin28a}) - \log_2(\text{TE}_{siLuc})$$

A **negative ΔTE** means: upon Lin28a knockdown, the RPF signal dropped more than the RNA signal (relative to the control). In other words, Lin28a normally *promotes* translation of that gene. The translation machinery is still recruited to the mRNA under normal conditions (Lin28a present), but once Lin28a is removed the ribosome engagement falls even if the mRNA itself persists. A **positive ΔTE** would indicate Lin28a normally *suppresses* translation of that gene (a rarer scenario). The global distribution of ΔTE therefore maps Lin28a's translational regulatory landscape at single-gene resolution.

---

## 3. Methodological Rationale: CPM vs. Raw Counts

### Why CPM normalisation is mandatory here

The original Mission 1/2 scripts computed a "ribosome density change" using raw read counts:

```python
# Original approach — incorrect without library-size correction
rden_change = (RPF_siLin28a / RNA_siLin28a) / (RPF_siLuc / RNA_siLuc)
```

This is only mathematically valid if the total library sizes of all four samples are identical. In this dataset they are not:

| Sample | Assigned reads |
|--------|---------------|
| RNA-siLuc | 28.1 M |
| RNA-siLin28a | 35.1 M |

A 25% difference in library size means that every raw TE estimate for siLin28a is systematically inflated by ~25% compared to siLuc. Without correction, we would falsely conclude that Lin28a knockdown *increases* TE on average — which is the opposite of the expected biology.

CPM (Counts Per Million) normalises each sample independently by dividing each gene's count by the total library size of that sample, then multiplying by 10⁶:

$$\text{CPM}_i = \frac{\text{count}_i}{\sum_j \text{count}_j} \times 10^6$$

This eliminates sequencing-depth differences between samples before any ratio is taken.

### Justification of the 0.5 CPM pseudocount

A pseudocount of 0.5 CPM is added only to the RNA denominator:

```python
TE_ctrl = RPF_ctrl_cpm / (RNA_ctrl_cpm + 0.5)
```

- **Why the denominator only:** Adding it to the numerator (RPF) would artificially elevate TE for genes with zero RPF signal, masking true translational silencing.
- **Why 0.5:** Half of 1 CPM is the smallest non-zero value detectable at a typical library depth of ~30 M reads. It is small enough not to distort TE estimates for well-expressed genes (where RNA CPM ≫ 0.5), but large enough to prevent −Inf values when computing log2(TE) for genes with near-zero RNA signal in one condition.

---

## 4. Step-by-Step Pipeline Description

### Step 1 — Load featureCounts matrix

**Input:** `data/read-counts.txt` — a tab-separated featureCounts v2.0.1 output with 55,359 Ensembl gene IDs × 6 samples (CLIP-35L33G, RNA-control, RNA-siLin28a, RNA-siLuc, RPF-siLin28a, RPF-siLuc).

**What is done:** Skip the first comment line (`#`) and load the matrix into a pandas DataFrame indexed by Geneid. Only the four experimental columns (RNA-siLuc, RNA-siLin28a, RPF-siLuc, RPF-siLin28a) are used for TE calculation.

**Why:** featureCounts provides un-normalised integer counts — the direct output of read-to-gene alignment. These raw integers are the correct starting point before any normalisation is applied.

### Step 2 — CPM normalisation

**What is done:** Each of the four sample columns is independently divided by its column sum and multiplied by 10⁶, producing four CPM vectors.

**Why:** As discussed in Section 3, library-size differences would confound the downstream TE ratio. CPM is the minimal correction needed (length normalisation is unnecessary because TE is a within-gene ratio).

### Step 3 — Low-count filtering

**Thresholds applied:**
- RNA CPM ≥ 1.0 in *both* siLuc and siLin28a conditions
- RPF raw read count ≥ 5 in *both* siLuc and siLin28a conditions

**Why:** Genes expressed at very low levels produce highly variable TE estimates because small absolute changes in count produce large log2 fold-changes. The RNA CPM filter removes genes where the mRNA is essentially absent. The RPF raw-count filter (applied before CPM conversion, to preserve integer precision) removes genes with insufficient ribosome coverage to compute a meaningful TE. Using ≥5 raw RPF reads is a standard threshold in ribosome profiling literature (Ingolia et al. 2009).

**Expected outcome:** ~35–40% of the 55,359 genes survive filtering, leaving ~18,000–22,000 genes for downstream analysis.

### Step 4 — Translation Efficiency and ΔTE calculation

**Formulas:**
```
TE_ctrl = RPF_ctrl_CPM / (RNA_ctrl_CPM + 0.5)
TE_kd   = RPF_kd_CPM   / (RNA_kd_CPM   + 0.5)
ΔTE     = log2(TE_kd) − log2(TE_ctrl)
ΔRNA    = log2(RNA_kd_CPM + 0.5) − log2(RNA_ctrl_CPM + 0.5)
ΔRPF    = log2(RPF_kd_CPM + 0.5) − log2(RPF_ctrl_CPM + 0.5)
```

**Biological meaning:** ΔTE isolates the translational component of Lin28a's effect. A gene where both ΔRNA and ΔRPF are equally negative has simply lower mRNA abundance (transcriptional effect); a gene where ΔRPF is more negative than ΔRNA has additionally lost translational efficiency (TE effect). This decomposition is only visible when RNA and RPF are compared at the same CPM scale.

### Step 5 — Cellular localisation annotation and merge

**Input:** `mouselocalization-20210507.txt` downloaded from `https://hyeshik.qbio.io/binfo/mouselocalization-20210507.txt`. Contains gene_id (Ensembl, no version suffix) and a `type` column with values such as `cytoplasm`, `nucleus`, `integral membrane`.

**What is done:** The Ensembl gene ID version suffix (e.g., `.2` in `ENSMUSG00000102693.2`) is stripped from the featureCounts index, and the localisation table is left-joined on `gene_id`.

**Why:** The localisation data was compiled from UniProt GO-slim subcellular location annotations for mouse proteins. Merging on the versionless Ensembl ID is necessary because the localisation file uses the base ID while featureCounts output carries the GENCODE version suffix.

### Step 6 — One-way ANOVA + pairwise t-tests

**What is done:**
1. ΔTE values are grouped into four sets: cytoplasm, nucleus, integral membrane, and unknown (no localisation annotation).
2. One-way ANOVA tests whether the mean ΔTE differs across the three annotated groups.
3. Pairwise Welch t-tests (unequal variances assumed) are performed for all three pairs with Bonferroni correction: α_adjusted = 0.05 / 3 = 0.0167.

**Why ANOVA first:** ANOVA controls the family-wise type I error rate before any pairwise comparison. Running t-tests directly on all pairs without ANOVA first would inflate the false-discovery rate.

**Why Bonferroni:** With only 3 comparisons, Bonferroni is not overly conservative and gives a clear, defensible significance threshold.

---

## 5. Expected Outputs

| File | Description |
|------|-------------|
| `output/analysis1_TE_scatter.png` | Scatter plot with ΔRNA (x-axis) vs. ΔRPF (y-axis), each point coloured by ΔTE (blue = more negative = Lin28a normally promotes translation). The diagonal dashed line represents ΔTE = 0. Points below the diagonal are translationally suppressed upon Lin28a knockdown. |
| `output/analysis1_localization.png` | Violin plot showing the ΔTE distribution for each localisation group. The median is shown as a horizontal line. The ANOVA F-statistic and p-value are annotated in the title. Groups with significantly different ΔTE from cytoplasmic proteins should be visible. |
| `output/analysis1_results.csv` | Per-gene table (filtered genes only) with columns: gene_id, localisation type, rna_ctrl_cpm, rna_kd_cpm, rpf_ctrl_cpm, rpf_kd_cpm, TE_ctrl, TE_kd, ΔTE, ΔRNA, ΔRPF. This CSV is the primary data table for any further downstream analysis. |

**Expected biological signal:**
- Global ΔTE mean should be negative (Lin28a globally promotes translation genome-wide)
- Integral membrane proteins should show a significantly different ΔTE compared to cytoplasmic proteins, consistent with Lin28a's known role in ER-associated translation (Cho et al. 2012, Cell)
- The scatter plot should show a cloud of points slightly below the ΔTE = 0 diagonal
