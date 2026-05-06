# Instructions for Claude Code Agent

## Commit & Push automatique

Après chaque modification de fichier, effectue automatiquement les étapes suivantes :

```bash
git add .
git commit -m "description courte de la modification"
git push
```

## Fichiers à ne JAMAIS committer

Les fichiers suivants sont trop lourds ou inutiles pour GitHub. Ne les ajoute jamais au commit :

```
*.bam
*.fastq
*.fastq.gz
*.tar
*.tar.gz
*.zip
colab-biolab/
*.sam
*.bed
*.bigwig
*.bw
```

Ces règles doivent aussi être présentes dans le fichier `.gitignore` à la racine du projet.

## .gitignore à maintenir

Vérifie que le fichier `.gitignore` contient toujours ces lignes :

```
# Données biologiques lourdes
*.bam
*.fastq
*.fastq.gz
*.sam
*.tar
*.tar.gz

# Setup Colab
colab-biolab/

# Fichiers système
.ipynb_checkpoints/
__pycache__/
*.pyc
.DS_Store
```

Si le fichier `.gitignore` n'existe pas, crée-le avant tout commit.

## Règles générales

- Ne committe que le **code** et les **notebooks** (`.ipynb`, `.py`, `.R`, `.md`)
- Ne committe jamais de fichiers de plus de **50 MB**
- Message de commit : décris brièvement ce qui a changé (ex: `"ajout analyse différentielle"`)
- Si `git push` échoue à cause d'un fichier trop lourd, ajoute-le au `.gitignore` et réessaie


TOKEN pour commit et push = [voir mémoire locale Claude — ne pas committer ici]

