# Étapes pour exécuter le projet de bioinformatique (Guide Mission 1)

Ce guide explique comment exécuter le notebook de projet de bioinformatique sur VS Code au lieu de Google Colab.

---

## Étape 1: Préparer l'environnement Python

**Objectif:** Installer les bibliothèques nécessaires pour le projet.

**Action requise:** Confirmer pour installer les packages Python (pandas, matplotlib, numpy, etc.)

---

## Étape 2: Installer les outils Bioinformatiques (Bioconda)

**Objectif:** Installer featureCounts (du package subread) pour compter les reads.

**Action requise:** Confirmer pour installer Bioconda et subread.

---

## Étape 3: Télécharger les données de séquençage

**Objectif:** Télécharger les fichiers BAM et les annotations GENCODE depuis le serveur.

**Action requise:** Confirmer pour télécharger les données (~2GB).

---

## Étape 4: Vérifier l'intégrité des données

**Objectif:** Vérifier que les fichiers téléchargés sont intacts avec MD5.

**Action requise:** Cette étape est automatique après le téléchargement.

---

## Étape 5: Exécuter featureCounts

**Objectif:** Compter les reads par gène à partir des fichiers BAM et GTF.

**Action requise:** Confirmer pour exécuter featureCounts.

---

## Step 6: Charger et analyser les données

**Objectif:** Charger les résultats dans pandas et créer les graphiques.

**Action requise:** Confirmer pour exécuter l'analyse et générer les figures.

---

## Étape 7: Créer le graphique final avec localisation protéique

**Objectif:** Créer un scatter plot coloré selon la localisation cellulaire des protéines.

**Action requise:** Confirmer pour générer le graphique final.

---

## Commandes rapides (si vous voulez exécuter manuellement)

```bash
# Installer les dépendances Python
pip install pandas matplotlib numpy

# Installer Bioconda et subread
conda install -y subread

# Télécharger les données
wget -O - --no-check-certificate https://hyeshik.qbio.io/binfo/binfo1-datapack1.tar | tar -xf -

# Télécharger GENCODE
wget --no-check-certificate -O gencode.gtf.gz http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M27/gencode.vM27.annotation.gtf.gz

# Exécuter featureCounts
featureCounts -a gencode.gtf -o read-counts.txt *.bam
```