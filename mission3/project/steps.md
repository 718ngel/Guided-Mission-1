# Mission 3 — CIMS Analysis (Crosslinking-Induced Mutation Sites)

> **Goal:** Identify protein-RNA binding sites by analyzing reverse transcriptase error rates in CLIP-seq data. Fully run locally in VS Code.

---

## Project Structure

```
project/
├── steps.md              ← this file
├── data/                 ← raw data (never committed)
│   ├── binfo1-datapack1/ ← BAM files and GTF
│   └── work/             ← derived working files
├── notebooks/            ← main Jupyter notebook
├── scripts/              ← Python and Bash scripts
└── results/              ← bedgraph files and UCSC screenshots
```

---

## Step 0 — VS Code Prerequisites

### 0.1 VS Code Extensions to Install
- **Python** (Microsoft)
- **Jupyter** (Microsoft)

### 0.2 Conda Environment

Create a dedicated environment with all required bioinformatics tools:

```bash
conda create -n binfo-mission3 python=3.11 -y
conda activate binfo-mission3
conda install -c bioconda -c conda-forge samtools bedtools bioawk -y
pip install pandas matplotlib jupyter
```

Verify installations:

```bash
samtools --version
bedtools --version
python -c "import pandas; print(pandas.__version__)"
```

### 0.3 Select Kernel in VS Code

Open the `.ipynb` notebook in VS Code, click **"Select Kernel"** in the top right, and choose the `binfo-mission3` environment.

---

## Step 1 — Data Download and Preparation

### 1.1 Download the Datapack

```bash
cd project/data
wget --no-check-certificate https://hyeshik.qbio.io/binfo/binfo1-datapack1.tar
tar -xf binfo1-datapack1.tar -C binfo1-datapack1/
```

> If `wget` fails on the first attempt, run the command again.

### 1.2 Verify MD5 Checksums

```bash
md5sum binfo1-datapack1/*
```

Expected values:
```
140aaf30bcb9276cc716f8699f04ddd6  CLIP-35L33G.bam
f1b3336ed7e2f97d562dcc71641251bd  CLIP-35L33G.bam.bai
328883a73d507eafbf5b60bd6b906201  RNA-control.bam
02073818e2f398a73c3b76e5169de1ca  RNA-control.bam.bai
...
```

### 1.3 Index BAM Files (if index is missing)

```bash
samtools index binfo1-datapack1/CLIP-35L33G.bam
```

---

## Step 2 — Find Genomic Coordinates

### 2.1 Search for Mirlet7g in the GTF File

```bash
grep -i mirlet7g data/binfo1-datapack1/gencode.gtf
```

> The region of interest is: **chr9:106056039-106056126** (mm39 genome).

Do the same for the two other genes:

```bash
grep -i mirlet7f-1 data/binfo1-datapack1/gencode.gtf
grep -i mirlet7d  data/binfo1-datapack1/gencode.gtf
```

> Record the coordinates for all three genes.

---

## Step 3 — Extract CLIP Reads by Region

For each gene, extract reads from the global BAM file into a local BAM file:

```bash
# Mirlet7g
samtools view -b -o data/work/CLIP-let7g.bam \
    data/binfo1-datapack1/CLIP-35L33G.bam \
    chr9:106056039-106056126
samtools view data/work/CLIP-let7g.bam | wc -l

# Mirlet7f-1 (use coordinates found in Step 2)
samtools view -b -o data/work/CLIP-let7f1.bam \
    data/binfo1-datapack1/CLIP-35L33G.bam \
    <chrom>:<start>-<end>

# Mirlet7d (use coordinates found in Step 2)
samtools view -b -o data/work/CLIP-let7d.bam \
    data/binfo1-datapack1/CLIP-35L33G.bam \
    <chrom>:<start>-<end>
```

---

## Step 4 — Generate the Pileup File

`samtools mpileup` summarizes, for each position, the bases read by all overlapping reads (coverage, substitutions, indels…).

```bash
# Mirlet7g
samtools mpileup data/work/CLIP-let7g.bam > data/work/CLIP-let7g.pileup
wc -l data/work/CLIP-let7g.pileup

# Mirlet7f-1
samtools mpileup data/work/CLIP-let7f1.bam > data/work/CLIP-let7f1.pileup

# Mirlet7d
samtools mpileup data/work/CLIP-let7d.bam  > data/work/CLIP-let7d.pileup
```

> The pileup also covers neighboring regions — this is expected since reads extend beyond the extracted region.

### 4.1 Filter Pileup to the Gene Region Only

```bash
awk '$2 >= 106056039 && $2 <= 106056126 { print $0; }' \
    data/work/CLIP-let7g.pileup > data/work/CLIP-let7g-gene.pileup
```

> Adapt the positions for each gene.

---

## Step 5 — Python Analysis in the Jupyter Notebook

Open `notebooks/mission3_analysis.ipynb` in VS Code with the `binfo-mission3` kernel.

### 5.1 Load Pileup into pandas

```python
import pandas as pd

cols = ['chrom', 'pos', '_ref', 'count', 'basereads', 'quals']
pileup = pd.read_csv('data/work/CLIP-let7g-gene.pileup', sep='\t', names=cols)
pileup.tail()
```

### 5.2 Clean Basereads

Remove special characters from the pileup that do not correspond to bases (insertion/deletion markers, read start/end tags):

```python
import re

to_remove = re.compile(r'[<>$*#^]|\+\d+[ACGTNacgtn]+|-\d+[ACGTNacgtn]+')

def clean_basereads(s):
    return to_remove.sub('', s)

pileup['matches'] = pileup['basereads'].apply(clean_basereads)
```

### 5.3 Count Bases per Position

```python
from collections import Counter

def count_bases(s):
    return Counter(c.upper() for c in s if c.upper() in 'ACGTN.,')

pileup['base_counts'] = pileup['matches'].apply(count_bases)
```

### 5.4 Calculate Shannon Entropy per Position

Shannon entropy measures base diversity at each position. A high entropy value indicates strong heterogeneity (potential mutation site).

```python
import math

def shannon_entropy(counts_dict):
    total = sum(counts_dict.values())
    if total == 0:
        return 0.0
    return -sum(
        (n / total) * math.log2(n / total)
        for n in counts_dict.values()
        if n > 0
    )

pileup['entropy'] = pileup['base_counts'].apply(shannon_entropy)
```

### 5.5 Export to bedgraph Format

The bedgraph format is: `chrom  start  end  value` (0-based coordinates).

```python
bedgraph_rows = []
for _, row in pileup.iterrows():
    start = int(row['pos']) - 1  # pileup is 1-based → bedgraph is 0-based
    end   = int(row['pos'])
    bedgraph_rows.append(f"{row['chrom']}\t{start}\t{end}\t{row['entropy']:.6f}")

header = 'track type=bedGraph name="CLIP-let7g Shannon entropy" visibility=full'
with open('results/CLIP-let7g-entropy.bedgraph', 'w') as f:
    f.write(header + '\n')
    f.write('\n'.join(bedgraph_rows) + '\n')
```

### 5.6 Repeat for Mirlet7f-1 and Mirlet7d

Re-run steps 5.1 to 5.5 with the corresponding pileup files and name output files accordingly:
- `results/CLIP-let7f1-entropy.bedgraph`
- `results/CLIP-let7d-entropy.bedgraph`

---

## Step 6 — Visualization on UCSC Genome Browser

1. Go to [UCSC Genome Browser — Mirlet7g region (mm39)](http://genome.ucsc.edu/cgi-bin/hgTracks?db=mm39&position=chr9%3A106056039-106056126)
2. Confirm the genome is set to **mm39**
3. At the bottom of the page, click **"add custom tracks"**
4. Click **"Choose File"** and upload `results/CLIP-let7g-entropy.bedgraph`
5. Click **"Submit"** then **"go to genome browser"**
6. Explore the track — zoom in, adjust track height
7. Go to **View → PDF/PS** to export a screenshot
8. Save the screenshot to `results/`
9. Repeat for Mirlet7f-1 and Mirlet7d

---

## Step 7 — Commit and Push to GitHub

```bash
git add mission3/project/notebooks/mission3_analysis.ipynb
git add mission3/project/scripts/
git add mission3/project/results/*.bedgraph mission3/project/results/*.pdf
git commit -m "Add Mission 3: CIMS Shannon entropy analysis for Mirlet7g, 7f-1, 7d

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>"
git push
```

> Never commit `.bam`, `.pileup`, `.fastq`, or `.tar` files — verify `.gitignore` excludes them.

---

## Summary of Output Files

| File | Description |
|---|---|
| `data/work/CLIP-let7g.bam` | CLIP reads extracted for Mirlet7g |
| `data/work/CLIP-let7g.pileup` | Raw pileup (full region) |
| `data/work/CLIP-let7g-gene.pileup` | Pileup filtered to gene coordinates |
| `results/CLIP-let7g-entropy.bedgraph` | Per-position entropy — UCSC format |
| `results/CLIP-let7g-ucsc.pdf` | UCSC Genome Browser screenshot |
| *(same for let7f-1 and let7d)* | |

---

## Useful Links

- [UCSC Genome Browser — Mirlet7g mm39](http://genome.ucsc.edu/cgi-bin/hgTracks?db=mm39&position=chr9%3A106056039-106056126)
- [bedgraph format — UCSC documentation](https://genome.ucsc.edu/goldenPath/help/bedgraph.html)
- [samtools mpileup documentation](http://www.htslib.org/doc/samtools-mpileup.html)
