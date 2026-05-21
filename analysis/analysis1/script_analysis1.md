# Agent Prompt — Analysis 1: Translation Efficiency (TE) & Cellular Localization

## Context

You are working in the directory:
`/Users/angela/Desktop/서울대학교/생물정보학 실습/`

The Git remote repository is:
`https://github.com/718ngel/Guided-Mission-1.git` (remote name: `origin`)

The reference specification for this analysis is:
`analysis/docs/REPORT.md` → section **"Analysis 1 — Translation Efficiency (TE) Change and Cellular Localization"**

The existing analysis scripts and notebooks are located in `analysis/scripts/`.

---

## Git Workflow — Follow These Rules for Every Interaction with the Repository

Before doing any work:
```bash
git pull origin main
```

After completing each major step (REPORT1.md written, script executed, result_analysis1.md written), stage and push:
```bash
git add analysis/analysis1/
git commit -m "analysis1: <short description of what was just done>"
git push origin main
```

All files created during this analysis must be placed inside `analysis/analysis1/` so they are tracked under the correct folder in the remote repository. Never commit files outside that directory as part of this task.

---

## Your Objectives (in order)

### Step 0 — Sync with the remote repository

```bash
git pull origin main
```

Verify that the `analysis/analysis1/` directory exists locally. If it does not, create it now. This directory is your working scope for the entire task.

---

### Step 1 — Read the specification

Read `analysis/docs/REPORT.md` in full. Extract all information relevant to Analysis 1:
- biological background and goal
- methodological corrections (CPM vs raw count)
- input files, filters, and parameters
- expected output files and figure descriptions
- script path: `analysis/scripts/analysis1_TE_localization.py`

---

### Step 2 — Set up the working directory

All work for this analysis must live in `analysis/analysis1/`. Create the following sub-structure if it does not already exist:

```
analysis/analysis1/
├── REPORT1.md          ← written by you before executing anything
├── scripts/            ← copy of the script you will run
├── output/             ← all generated figures and CSV files
└── result_analysis1.md ← written by you after execution
```

Then commit and push the empty scaffold:
```bash
git add analysis/analysis1/
git commit -m "analysis1: initialise directory structure"
git push origin main
```

---

### Step 3 — Write REPORT1.md (before running any code)

Write `analysis/analysis1/REPORT1.md` containing:

1. **Task list** — a numbered checklist of every concrete action you will take
2. **Biological objective** — in 2–3 paragraphs: what TE measures, why Lin28a knockdown is the perturbation of interest, and what a negative ΔTE means biologically
3. **Methodological rationale** — explain why CPM normalisation is used instead of raw counts, and justify the pseudocount of 0.5 CPM
4. **Step-by-step pipeline description** — for each step in the Methods table of REPORT.md, explain *what* is done and *why* (input, transformation, output, biological meaning)
5. **Expected outputs** — list every figure and file that should be produced, with a one-sentence description of what each should show

Then commit and push:
```bash
git add analysis/analysis1/REPORT1.md
git commit -m "analysis1: write REPORT1.md — objectives and pipeline plan"
git push origin main
```

---

### Step 4 — Execute the analysis

Run the script:
```bash
source analysis/.venv/bin/activate   # or Term_Project/.venv/bin/activate if that is the correct venv
python analysis/scripts/analysis1_TE_localization.py
```

If the script does not exist or fails, diagnose the error, fix it, and re-run. Do not skip execution — results must be real, not hypothetical.

Copy all output files (`.png`, `.csv`) into `analysis/analysis1/output/`.

Then commit and push:
```bash
git add analysis/analysis1/scripts/ analysis/analysis1/output/
git commit -m "analysis1: run script — add figures and CSV results"
git push origin main
```

---

### Step 5 — Write result_analysis1.md (after execution)

Write `analysis/analysis1/result_analysis1.md` containing:

1. **Pipeline steps actually executed** — what commands were run, in what order, with what parameters
2. **Methods summary** — CPM normalisation formula, TE formula, ΔTE formula, filter thresholds, statistical tests used (ANOVA + Bonferroni pairwise t-tests)
3. **Expected results** (from REPORT.md) — what the analysis was supposed to show
4. **Obtained results** — describe what the actual figures and CSV show: distribution of ΔTE, which localisation group showed the largest effect, ANOVA p-value, pairwise test results
5. **Critical commentary** — your own scientific assessment:
   - Strengths of the approach (what makes this analysis rigorous)
   - Weaknesses or limitations (what biases remain, what data is missing, what assumptions are not verified)
   - Comparison with the original Mission 1/2 approach (what the CPM correction changes in practice)
   - Suggestions for improvement if this were a published paper

Then commit and push the final report:
```bash
git add analysis/analysis1/result_analysis1.md
git commit -m "analysis1: write result_analysis1.md — results, methods, critical commentary"
git push origin main
```

---

## Constraints

- Remote repository: `https://github.com/718ngel/Guided-Mission-1.git` (remote: `origin`, branch: `main`)
- Always `git pull origin main` before starting any step.
- Every completed step must end with a `git add` + `git commit` + `git push origin main`.
- Do not modify `analysis/docs/REPORT.md` or any file outside `analysis/analysis1/`.
- All output figures must be saved to `analysis/analysis1/output/`, not elsewhere.
- `REPORT1.md` must be written and pushed **before** running the script.
- `result_analysis1.md` must reflect **actual** execution results, not hypothetical ones.
- If a dependency (Python package, data file) is missing, report it explicitly and explain how to resolve it before proceeding.
