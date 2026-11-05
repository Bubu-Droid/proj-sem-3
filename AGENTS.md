# AGENTS.md

Statistics course project: proportional-odds and EM (Dawid–Skene) modelling of
LIDC-IDRI lung-nodule malignancy ratings. `PLAN.md` is the brief, `STEPS.md`
the step-by-step walkthrough (`STEPS.html` is generated from it).

Keep this file current: when you change the pipeline, dependencies or data
conventions, update it in the same commit.

## Environment

- Python is managed with **uv** only (`uv add`, `uv run`); never `pip install`.
- Python is pinned to 3.11 with `numpy<1.24`, `sqlalchemy<2.0`,
  `setuptools<80` (and matching pandas/scipy/scikit-image) because pylidc
  0.2.3 needs `pkg_resources` and the old NumPy/SQLAlchemy APIs. Don't lift
  these pins.
- Python is used only for extraction (`python/`); all statistics are in R
  (`R/`). Run R scripts from the repo root: `Rscript R/<script>.R`.

## Pipeline

```bash
uv run python python/extract.py   # -> data/raw/annotations.csv (slow; run once)
Rscript R/01_clean.R              # -> data/processed/*.csv
Rscript R/check_data.R            # sanity checks; exits non-zero on failure
Rscript R/eda.R                   # -> figures/*.pdf, prints numbers used in main.tex
```

`data/raw/annotations.csv` is the master file and is never edited by hand.
Re-run `R/check_data.R` after any change to extraction or cleaning.

Installing R packages from source (`ordinal`, `car`, and their dependencies)
needs a Fortran compiler: `sudo pacman -S gcc-fortran`.

## Report

All LaTeX lives in `core/`: `main.tex` uses `statproj.sty` (title page,
minted) and `refs.bib` with plain BibTeX, and includes figures from
`../figures/`. Build from inside `core/`:

```bash
cd core && latexmk -pdf -shell-escape main.tex
```

The user's global latexmkrc puts aux files in `core/buildtmp/` (ignored).
`core/main.pdf` is committed. Numbers in `main.tex` are typed in by hand
from `R/eda.R` and `R/01_clean.R` output, so update them when the data or scripts change. Do not double-space:
keep the style's default one-and-a-half spacing. Author, roll number
(`\rollno`) and date are placeholders in `main.tex`; never hard-code personal
details in `statproj.sty`.

## Data conventions

- Unit of analysis is the **nodule** (`nodule_id` = `<scan_id>_<cluster>`),
  not the annotation. Train/test split is by **patient**.
- `reader_slot` / `R1..R4` order is arbitrary: LIDC reader identities are not
  consistent across scans.
- `calcification` and `internalStructure` are nominal and excluded from
  modelling; they stay in `annotations.csv` and `nodules.csv` only.
- `mal_median` rounds half-grade medians **up**; `mal_median_tie` flags them.
- Processed CSVs hold unstandardised values; standardise in the modelling
  script using training-set means/SDs.

## Git

Conventional Commits (`feat`, `fix`, `docs`, `build`, `chore`, `refactor`,
with an optional scope such as `feat(clean):`). Keep commits small and
focused: one logical change per commit, never one large catch-all commit.

After editing `STEPS.md`, regenerate the HTML:

```bash
pandoc STEPS.md -f gfm+tex_math_dollars -s --toc --toc-depth=2 \
  --math-method=mathjax --metadata pagetitle="Project Walkthrough" -o STEPS.html
```
