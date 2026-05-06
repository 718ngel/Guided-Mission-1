# Mission Instructions for Claude Code Agent

This file captures all standing orders given by the user for this project.
Any future agent working on a mission of this type should follow these rules without the user having to repeat them.

---

## 1. Execution Protocol

- Read `docs/steps.md` at the start of every mission to understand the pipeline steps.
- **Ask the user before executing each step.** Wait for explicit confirmation ("oui", "yes", "ok", etc.) before proceeding.
- If a step requires installing libraries or tools, install them automatically — no need to ask permission.
- Execute steps in order. Do not skip steps.

---

## 2. Git: Commit & Push Rules

After **every file modification**, run:

```bash
git add <modified files>   # never use git add -A or git add . blindly
git commit -m "short description of what changed"
git push
```

- Commit messages must be written **in English**.
- Describe **what changed and why**, not just what the file is.
- Always use a heredoc for multi-line messages:
  ```bash
  git commit -m "$(cat <<'EOF'
  Short title

  - Bullet detail 1
  - Bullet detail 2

  Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
  EOF
  )"
  ```
- Always add `Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>` to commit messages.

### If `git push` fails with HTTP 400 / chunked encoding error

The remote URL token may be expired. Retrieve the current token from Claude memory
(`~/.claude/projects/.../memory/project_github_token.md`) and update the remote:

```bash
git remote set-url origin "https://NEW_TOKEN@github.com/718ngel/Guided-Mission-1.git"
git push
```

### If GitHub blocks the push with "secret scanning" / GH013 error

A secret (API key, token) is present in a tracked file. Remove it before pushing:
1. Replace the secret in the file with a placeholder.
2. `git add` the file and `git commit --amend --no-edit`.
3. Retry `git push`.
4. Store the secret in Claude memory instead (never in version-controlled files).

---

## 3. Files Never to Commit

The following file types are too large or irrelevant for GitHub. **Never stage them.**

```
*.bam          # aligned read files (can be several GB)
*.bai          # BAM index files
*.fastq        # raw sequencing reads
*.fastq.gz
*.sam
*.gtf          # genome annotation (can be 800 MB+)
*.gtf.gz
*.tar
*.tar.gz
*.zip
*.bigwig / *.bw
colab-biolab/
```

Also avoid any single file larger than **50 MB**.
If a large file was accidentally staged, run `git rm --cached <file>` and add it to `.gitignore`.

---

## 4. .gitignore Maintenance

Always keep `.gitignore` up to date with the patterns above.
If the `.gitignore` does not exist, create it before the first commit.

---

## 5. Project File Structure

Organize files into these directories:

```
project-root/
├── data/          ← input data files (BAM, GTF, etc.) — NOT tracked by git
├── docs/          ← documentation: steps.md, REPORT.md, mission.md
├── notebooks/     ← Jupyter notebooks (.ipynb)
├── output/        ← generated figures and result files (.png, .pdf, .csv)
└── scripts/       ← analysis scripts (.py, .R, .sh)
```

When scripts reference `data/` or `output/` using `os.path.dirname(__file__)`, update the path to go one level up after moving them into `scripts/`:

```python
ROOT = os.path.join(os.path.dirname(__file__), "..")
DATA_DIR = os.path.join(ROOT, "data")
OUTPUT_DIR = os.path.join(ROOT, "output")
```

---

## 6. Final Deliverable

At the end of every mission, produce a `docs/REPORT.md` file in **English** that includes:

1. **Project overview** — biological question and data types used
2. **Pipeline summary** — diagram or bullet list of the full workflow
3. **Data integrity** — MD5 / checksum results for input files
4. **Tool output summary** — e.g., featureCounts read assignment rates per sample, with interpretation
5. **Quantitative results** — summary statistics (mean, median, std, min, max) for key metrics
6. **Biological interpretation** — what the numbers mean in context
7. **Output files table** — list of all generated files and what they contain
8. **Reproduction instructions** — exact commands to rerun the full pipeline

Commit and push `REPORT.md` after writing it.

---

## 7. Environment

- Python virtual environment: `.venv/` at project root (activate with `source .venv/bin/activate` or use `python3` directly if already active)
- Conda environment for bioinformatics tools: `subread_env` (use `conda run -n subread_env <command>`)
- featureCounts: available via `conda run -n subread_env featureCounts`
- GitHub remote: `https://github.com/718ngel/Guided-Mission-1.git`
- Token: stored in Claude memory — see `~/.claude/projects/.../memory/project_github_token.md`
