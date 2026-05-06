# Instructions for Claude Code Agent

## Commit & Push automatique

After every file modification, automatically run:

```bash
git add <modified files>   # never use git add -A or git add . blindly
git commit -m "short description in English of what changed"
git push
```

**Remote:** `https://github.com/718ngel/Guided-Mission2.git`
**Token:** stored in Claude memory — retrieve with:
`cat ~/.claude/projects/-Users-angela-Desktop----------Term-Project/memory/project_github_token.md`

Set the remote URL with token when needed:
```bash
git remote set-url origin "https://TOKEN@github.com/718ngel/Guided-Mission2.git"
```

Commit message rules:
- Written in **English**
- Describe what changed and why, not just the file name
- Always append: `Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>`

If `git push` fails with HTTP 400 or authentication error, retrieve the token from memory and update the remote URL above.

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

