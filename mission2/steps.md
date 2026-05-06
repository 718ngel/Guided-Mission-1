# Mission 2 – Steps to Reproduce in VS Code

> **Objectif biologique** : Observer la distribution des 5'-ends de fragments protégés par les ribosomes (RPF) autour des start codons, afin de reproduire la Figure S5A du papier.

---

## Structure du projet

```
mission2/
├── data/          ← fichiers d'entrée (BAM, GTF) — non trackés par git
├── work/          ← fichiers intermédiaires générés par le pipeline
├── output/        ← figure finale (figure_S5A.png)
├── scripts/
│   └── analysis.py   ← script principal à lancer depuis VS Code
└── steps.md          ← ce fichier
```

---

## Étape 0 — Prérequis : installer les outils

Ouvre un terminal dans VS Code (`Ctrl+`` ` `` ` ou `Terminal > New Terminal`).

### Outils bioinformatiques (via Conda / Bioconda)

```bash
conda install -y -c bioconda -c conda-forge bedtools bioawk samtools
```

Vérifie que tout fonctionne :

```bash
bedtools --version
samtools --version
bioawk --version
```

### Bibliothèques Python

```bash
pip install pandas matplotlib
```

---

## Étape 1 — Télécharger les données

Place les fichiers suivants dans le dossier `data/` :

### 1a. Fichiers BAM (données RPF)

Télécharge depuis la source du cours ou copie-les :

```bash
# Option A — depuis le datapack du cours (en local)
cp /chemin/vers/binfo1-datapack1/RPF-siLuc.bam data/

# Option B — télécharger le datapack complet
wget -O - --no-check-certificate https://hyeshik.qbio.io/binfo/binfo1-datapack1.tar \
  | tar -C data/ -xf -
```

Fichiers requis dans `data/` :
- `RPF-siLuc.bam`
- `RPF-siLuc.bam.bai`  _(index — pas strictement requis pour ce script)_

### 1b. Annotation GENCODE (GTF de la souris vM27)

```bash
wget --no-check-certificate \
  -O data/gencode.gtf.gz \
  http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M27/gencode.vM27.annotation.gtf.gz

gunzip data/gencode.gtf.gz
```

### 1c. Vérification des checksums MD5 (optionnel mais recommandé)

```bash
md5sum data/RPF-siLuc.bam
# attendu : f2eebf50943024d0116c9cd3e744c707
```

---

## Étape 2 — Extraire les start codons (+ strand, TSL 1)

> **Pourquoi** : on veut les positions exactes des start codons sur le brin +, uniquement pour les transcrits bien supportés (transcript support level 1) pour éviter de surestimer les comptages.

```bash
grep '	start_codon	.*	+	.*transcript_support_level "1"' data/gencode.gtf \
  | sed -e 's/	[^	]*transcript_id "\([^"]*\)".*$/	\1/g' \
  > work/gencode-start.gtf
```

⚠️ Les espaces dans la commande `grep` sont des **tabulations** (Tab), pas des espaces normaux.

Résultat : `work/gencode-start.gtf` — colonnes 1–8 standard GTF + col. 9 = transcript_id

---

## Étape 3 — Extraire les exons du brin +

> **Pourquoi** : on n'a besoin que des exons sur le brin + pour simplifier l'analyse (le brin − nécessite une transformation de coordonnées supplémentaire).

```bash
grep '	exon	.*	+	' data/gencode.gtf \
  | sed -e 's/	[^	]*transcript_id "\([^"]*\)".*$/	\1/g' \
  > work/gencode-plusexon.gtf
```

Résultat : `work/gencode-plusexon.gtf` — même format, mais uniquement les features `exon` du brin +.

---

## Étape 4 — Trouver les exons contenant un start codon

> **Pourquoi** : on filtre les exons en ne gardant que ceux qui chevauchent le start codon du même transcrit, puis on convertit en BED 0-based pour faciliter les calculs.

```bash
bedtools intersect \
  -a work/gencode-start.gtf \
  -b work/gencode-plusexon.gtf \
  -wa -wb \
  | awk -F'	' -v OFS='	' '$9 == $18 { print $10, $13-1, $14, $18, $4-1, $16; }' \
  | sort -k1,1 -k2,3n -k4,4 \
  > work/gencode-exons-containing-startcodon.bed
```

Format de sortie (colonnes) :

| col | contenu |
|-----|---------|
| 1 | chromosome |
| 2 | début exon (0-based) |
| 3 | fin exon |
| 4 | transcript_id |
| 5 | position start codon (0-based) |
| 6 | strand (+) |

---

## Étape 5 — Filtrer le BAM RPF (brin +, longueur ≥ 25 nt)

> **Pourquoi** : les RPF font normalement 25–30 nt. Les reads plus courts sont probablement des contaminations rRNA/tRNA. `-F20` garde uniquement le brin + (élimine les reads reverse-strand et unmapped).

```bash
(samtools view -H data/RPF-siLuc.bam; \
 samtools view -F20 data/RPF-siLuc.bam \
   | bioawk -c sam '{ if (length($seq) >= 25) print $0; }') \
| samtools view -b -o work/filtered-RPF-siLuc.bam
```

Vérifie la taille obtenue :

```bash
ls -lh data/RPF-siLuc.bam work/filtered-RPF-siLuc.bam
```

---

## Étape 6 — Calculer la couverture des 5'-ends (bedgraph)

> **Pourquoi** : `-5` dit à bedtools de ne compter que l'extrémité 5' de chaque read (= position d'entrée du ribosome), pas toute la longueur du read.

```bash
bedtools genomecov \
  -ibam work/filtered-RPF-siLuc.bam \
  -bg -5 \
  > work/fivepcounts-RPF-siLuc.bed
```

Format : `chr  start  end  count` (1 ligne par position avec au moins 1 read)

---

## Étape 7 — Intersecter avec les exons contenant un start codon

> **Pourquoi** : on ne veut garder que les 5'-ends qui tombent dans un exon portant un start codon, pour ensuite calculer la position relative de chaque read par rapport au start codon.

```bash
bedtools intersect \
  -a work/fivepcounts-RPF-siLuc.bed \
  -b work/gencode-exons-containing-startcodon.bed \
  -wa -wb -nonamecheck \
  > work/fivepcounts-filtered-RPF-siLuc.txt
```

Format de sortie (10 colonnes) :

| cols 1–4 | 5'-end position + count (depuis fivepcounts) |
| cols 5–10 | exon + start codon info (depuis le BED) |

---

## Étape 8 — Tracer la Figure S5A (Python)

Ouvre [scripts/analysis.py](scripts/analysis.py) dans VS Code.

**Option A — lancer tout le pipeline d'un coup**

```bash
python scripts/analysis.py
```

Ce script exécute toutes les étapes 2 à 7 ci-dessus via `subprocess`, puis génère le graphique.

**Option B — tracer uniquement le graphique** (si les fichiers intermédiaires existent déjà)

Colle ce code dans un fichier `scripts/plot_only.py` et lance-le :

```python
import pandas as pd
import matplotlib.pyplot as plt

COLS = [
    "read_chr", "read_start", "read_end", "count",
    "exon_chr", "exon_start", "exon_end",
    "transcript_id", "startcodon_pos", "strand",
]

WINDOW  = (-50, 100)
CHUNK   = 500_000
profile = pd.Series(dtype=float)

for chunk in pd.read_csv(
    "work/fivepcounts-filtered-RPF-siLuc.txt",
    sep="\t", header=None, names=COLS, chunksize=CHUNK,
):
    chunk["rel_pos"] = chunk["read_start"] - chunk["startcodon_pos"]
    chunk = chunk[
        (chunk["rel_pos"] >= WINDOW[0]) & (chunk["rel_pos"] <= WINDOW[1])
    ]
    profile = profile.add(chunk.groupby("rel_pos")["count"].sum(), fill_value=0)

profile = profile.sort_index()

fig, ax = plt.subplots(figsize=(12, 4))
ax.bar(profile.index, profile.values, width=1, color="steelblue", alpha=0.85)
ax.axvline(0, color="red", linestyle="--", linewidth=1.5, label="Start codon (A of AUG)")
ax.set_xlabel("Position relative to start codon (nt)")
ax.set_ylabel("RPF 5'-end count")
ax.set_title("RPF 5'-end distribution around start codons\n(+ strand, RPF-siLuc, TSL-1 transcripts)")
ax.legend()
plt.tight_layout()
plt.savefig("output/figure_S5A.png", dpi=150)
print("Figure saved → output/figure_S5A.png")
plt.show()
```

---

## Ce que représente la figure

- **Axe X** : position en nt par rapport au A du codon AUG (= 0)
- **Axe Y** : nombre de 5'-ends de reads RPF à cette position
- **Signal attendu** : un pic fort autour de −15 nt (empreinte du ribosome placé sur l'AUG), puis un patron en « dents de scie » à période 3 le long de la CDS (codons triplets)
- **Différence avec Fig. S5A du papier** : ici on ne traite que le brin + et les TSL 1 ; les valeurs absolues diffèrent, mais le profil qualitatif doit être similaire.

---

## Résumé rapide des commandes

```bash
# 0. Installation (une seule fois)
conda install -y -c bioconda bedtools bioawk samtools
pip install pandas matplotlib

# 1. Données → data/
#    RPF-siLuc.bam  +  gencode.gtf

# 2–7. Pipeline complet
python scripts/analysis.py

# Résultat → output/figure_S5A.png
```
