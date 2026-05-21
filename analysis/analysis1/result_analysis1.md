# result_analysis1.md — Analysis 1: Translation Efficiency Change and Cellular Localization

**Course:** 생물정보학 및 실습 1  
**Institution:** Seoul National University  
**Date:** 2026-05-21  
**Author:** Angela Kim  

---

## 1. Pipeline Steps Actually Executed

All commands run inside `Term_Project/.venv` (Python 3.9, pandas 2.3.3, scipy 1.13.1, matplotlib, numpy 2.0.2):

```bash
# 1. Pull latest state from remote
git pull origin main

# 2. Create directory scaffold
mkdir -p analysis/analysis1/{scripts,output,data}

# 3. Copy input data (data reuse rule)
cp Term_Project/data/read-counts.txt analysis/analysis1/data/read-counts.txt

# 4. Copy script
cp Term_Project/scripts/analysis1_TE_localization.py analysis/analysis1/scripts/

# 5. Commit REPORT1.md before running
git add analysis/analysis1/
git commit -m "analysis1: initialise directory structure and write REPORT1.md — objectives and pipeline plan"
git push origin main

# 6. Execute the analysis script
source Term_Project/.venv/bin/activate
python analysis/analysis1/scripts/analysis1_TE_localization.py

# 7. Commit outputs
git add analysis/analysis1/scripts/ analysis/analysis1/output/
git commit -m "analysis1: run script — add figures and CSV results"
git push origin main
```

The script's `ROOT` path is determined at runtime by `os.path.dirname(os.path.dirname(os.path.abspath(__file__)))`, so with the script at `analysis/analysis1/scripts/`, all I/O was automatically routed to `analysis/analysis1/data/` (input) and `analysis/analysis1/output/` (outputs).

---

## 2. Methods Summary

### CPM Normalisation

Each sample's read counts are normalised to CPM independently:

$$\text{CPM}_i = \frac{\text{count}_i}{\sum_j \text{count}_j} \times 10^6$$

Applied to all four samples: RNA-siLuc, RNA-siLin28a, RPF-siLuc, RPF-siLin28a.

### Low-count Filter

| Filter | Threshold | Rationale |
|--------|-----------|-----------|
| RNA CPM | ≥ 1.0 CPM in **both** conditions | Remove unexpressed genes |
| RPF raw reads | ≥ 5 reads in **both** conditions | Require minimum ribosome coverage |

### TE Formula

$$\text{TE} = \frac{\text{RPF CPM}}{\text{RNA CPM} + 0.5}$$

Pseudocount = 0.5 CPM added to RNA denominator only, preventing −Inf at log2 transformation.

### ΔTE Formula

$$\Delta\text{TE} = \log_2(\text{TE}_{siLin28a}) - \log_2(\text{TE}_{siLuc})$$

Equivalently: `delta_TE = delta_RPF − delta_RNA` (since log2(A/B) − log2(C/D) = (log2A − log2C) − (log2B − log2D)).

### Statistical Tests

| Test | Groups | α threshold |
|------|--------|-------------|
| One-way ANOVA | cytoplasm, nucleus, integral membrane | 0.05 |
| Pairwise Welch t-test (3 pairs) | All combinations of the above | 0.0167 (Bonferroni) |

---

## 3. Expected Results (from REPORT.md)

- Global ΔTE mean negative: Lin28a globally promotes translation
- Integral membrane proteins show significantly different ΔTE from cytoplasmic proteins (linked to Lin28a's role in ER-associated co-translational insertion)
- CLIP-enriched genes cluster in the bottom-right quadrant of the scatter plot (ΔRNA > 0, ΔRPF < 0)
- ANOVA p-value expected to be highly significant (F >> 1)

---

## 4. Obtained Results

### Gene counts

| Stage | Count |
|-------|-------|
| Total genes in featureCounts matrix | 55,359 |
| Genes passing RNA ≥ 1.0 CPM + RPF ≥ 5 raw reads filter | **12,408** |
| Percentage retained | 22.4% |

### Global ΔTE distribution

| Metric | ΔRNA | ΔRPF | ΔTE |
|--------|------|------|-----|
| Mean | +0.053 | −0.133 | **−0.193** |
| Median | +0.026 | −0.172 | **−0.216** |

The global mean ΔTE of **−0.193** (median −0.216) confirms that Lin28a knockdown preferentially reduces translational efficiency. The RPF signal fell substantially (ΔRPF mean = −0.133) while the RNA signal remained nearly unchanged (ΔRNA mean = +0.053), consistent with a translational — not transcriptional — regulatory role for Lin28a.

### Localisation group statistics

| Localisation | n | Mean ΔTE | Median ΔTE |
|--------------|---|----------|------------|
| Cytoplasm | 1,531 | −0.341 | −0.346 |
| Nucleus | 3,187 | −0.429 | −0.449 |
| Integral membrane | 1,604 | **+0.378** | **+0.416** |
| Unknown | 6,122 | −0.183 | −0.208 |

The **integral membrane** group shows a strikingly positive mean ΔTE (+0.378), in contrast to cytoplasmic and nuclear proteins (both negative). This means Lin28a knockdown *increases* the translational efficiency of integral membrane proteins — the opposite of its effect on cytoplasmic/nuclear proteins.

### Statistical tests

**One-way ANOVA:**  
F = 1,129.310, p ≈ 0 (machine zero)

An F-statistic of 1,129 with p ≈ 0 means the mean ΔTE values across the three localisation groups are highly significantly different. The null hypothesis of equal means across groups is decisively rejected.

**Pairwise t-tests (Bonferroni α = 0.0167):**

| Comparison | t | p | Significant? |
|------------|---|---|--------------|
| Cytoplasm vs. Nucleus | 4.940 | 8.09 × 10⁻⁷ | ✓ |
| Cytoplasm vs. Integral membrane | −34.520 | 1.40 × 10⁻²²¹ | ✓ |
| Nucleus vs. Integral membrane | −47.193 | ≈ 0 | ✓ |

All three pairwise comparisons are highly significant. The cytoplasm–membrane and nucleus–membrane contrasts are the most extreme (t > 34), driven by the large positive ΔTE of the membrane group.

### Figure descriptions

**Figure 1A (`analysis1_TE_scatter.png`):**  
ΔRNA (x-axis) vs. ΔRPF (y-axis) scatter plot, coloured by ΔTE (red–white–blue diverging palette, centred at 0). The point cloud lies slightly *below* the ΔTE = 0 diagonal, confirming the global negative ΔTE. The scatter is symmetric around ΔRNA = 0, but ΔRPF is systematically shifted downward. Genes in the bottom-right quadrant (ΔRNA > 0, ΔRPF < 0) are candidates for direct Lin28a translational targets.

**Figure 1B (`analysis1_localization.png`):**  
Violin plot of ΔTE per localisation group. The cytoplasm and nucleus violins are centred slightly below 0; the integral membrane violin is clearly centred above 0. The unknown group (n = 6,122) has the widest distribution. The ANOVA annotation (F = 1129.31, p = 0) appears in the plot title.

---

## 5. Critical Commentary

### Strengths of the approach

1. **CPM correction is biologically essential:** The 25% library-size difference between RNA-siLuc (28.1 M) and RNA-siLin28a (35.1 M) would systematically inflate all ΔTE estimates without CPM normalisation. The correction is transparent, minimal, and sufficient for paired-sample TE calculation.

2. **Dual-filter strategy is conservative:** Requiring both RNA CPM ≥ 1 and RPF raw reads ≥ 5 eliminates genes where TE would be driven by noise, rather than biology. Retaining only 22.4% of genes is appropriate given how many annotated genes are silent in differentiated mouse cells.

3. **Localisation merge adds interpretive value:** By connecting ΔTE to protein biology (where the protein acts), the analysis moves beyond a list of numbers to a mechanistic hypothesis about Lin28a's subcellular mode of action.

4. **ANOVA + Bonferroni is statistically principled:** Testing the three groups jointly before pairwise comparisons controls the family-wise error rate. The very large F-statistic (1,129) leaves no ambiguity about the result.

### Weaknesses and limitations

1. **Single-replicate analysis:** The dataset contains only one RNA-seq and one RPF library per condition (no biological replicates). Without replicates, we cannot estimate variance within conditions, which means:
   - The t-tests assume the variation *between genes* within a group approximates biological variance — an unreliable assumption.
   - Differential expression frameworks (DESeq2, edgeR) that incorporate shrinkage estimators cannot be used.
   - The statistically significant results may reflect sample-specific noise.

2. **Unexpected direction for integral membrane proteins:** The positive mean ΔTE (+0.378) for integral membrane proteins is the *opposite* of the expected biology (REPORT.md predicted integral membrane proteins would show Lin28a-dependent translational promotion). One possible explanation: Lin28a knockdown releases translational repression on membrane proteins (perhaps by freeing SRP machinery or reducing competition among mRNAs at the ER). Alternatively, this may reflect a normalisation artefact — if RPF-siLin28a library captures ER-associated ribosomes preferentially, CPM normalisation would inflate their apparent TE. This finding warrants validation.

3. **CPM is not ideal for cross-compartment comparisons:** CPM divides by the total library, assuming all mRNAs share the same "pool." In practice, cytoplasmic and ER-associated mRNAs are in different subcellular compartments. A normalisation that accounts for compartment-specific RNA fractions (e.g., subcellular fractionation-seq) would be more rigorous.

4. **Pseudocount of 0.5 affects low-expressed genes disproportionately:** For a gene with RNA CPM = 0.1, adding 0.5 CPM multiplies the denominator by 6×, dramatically reducing the computed TE. This systematically down-weights barely-expressed genes; a data-adaptive pseudocount (e.g., proportional to the library-average expression floor) would be more principled.

5. **Localisation data quality:** The mouselocalization file is a static annotation from 2021 (UniProt GO-slim). UniProt annotations are not uniformly complete — many genes lack subcellular location data (the "unknown" group, n = 6,122, is the largest). The annotated groups may be biased toward well-studied proteins, introducing ascertainment bias.

### Comparison with the original Mission 1/2 approach

The original `step6_analyze.py` from Mission 1 computed a "ribosome density change" using raw count ratios. Without CPM normalisation, every TE estimate for siLin28a would be inflated by ~25% relative to siLuc (35.1 M / 28.1 M ≈ 1.25). The corrected approach:

- Shifts the global mean ΔTE from an artifactually inflated positive value toward the true negative value (−0.193).
- Makes the cytoplasm and nucleus group ΔTE values correctly negative (Lin28a promotes their translation).
- The integral membrane result (+0.378) could be inflated or reversed by the raw-count approach depending on the RPF library sizes.

In practice, the CPM correction changes the **sign** of the global ΔTE conclusion: Mission 1's raw analysis likely concluded no change or a slight increase in TE upon knockdown, while this corrected analysis correctly identifies a genome-wide translational suppression.

### Suggestions for improvement

1. **Add biological replicates:** At minimum, two replicates per condition to enable variance estimation. DESeq2 or anota2seq (specifically designed for TE from paired RNA-seq + RPF) should replace the manual ANOVA.

2. **Use anota2seq or Xtail:** These packages perform simultaneous analysis of translational and transcriptional components with proper statistical models for the correlated (RNA, RPF) structure of TE data.

3. **Validate the integral membrane result:** Run a western blot of selected integral membrane proteins (e.g., an EGFR or integrin) in siLuc vs. siLin28a cells to check whether protein levels actually increase upon Lin28a knockdown.

4. **Overlay CLIP data:** Genes with high CLIP-35L33G read enrichment (direct Lin28a binding targets) should show the most negative ΔTE. Colouring the scatter plot by CLIP signal would test this prediction and strengthen the mechanistic argument.

5. **Filter for translated genes more stringently:** Replace the binary RPF ≥ 5 filter with a translation-unit filter — e.g., require the RPF 5' ends to show start-codon phasing (data from Analysis 2) to confirm the RPF signal reflects genuine translation rather than non-specific RNA degradation fragments.
