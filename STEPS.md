# Project Walkthrough: Ordinal Modelling of LIDC-IDRI Nodule Malignancy

A step-by-step guide from an empty folder to the final defense. Read alongside `PLAN.md`, which states *what* the project is. This file says *how* to do it and in what order.

Conventions used below:

- **Python** is used only for data extraction (`pylidc` is Python-only).
- **R** is used for all statistics (`MASS::polr`, `ordinal`, `brant` are mature and accept case weights, which the EM M-step needs).
- `☐` marks a concrete deliverable that should exist before you move on.

---

## Step 0 — Set up the project (Day 1)

### 0.1 Folder layout

```text
proj-sem-3/
├── PLAN.md
├── STEPS.md
├── data/
│   ├── raw/            # CSVs straight out of pylidc (never edited by hand)
│   └── processed/      # nodule-level analysis tables
├── python/
│   └── extract.py
├── main.tex
├── main.pdf
├── refs.bib
├── statproj.sty
├── isilogo.png
├── codes/
│   ├── compileplots.sh
│   ├── r-codes/
│   │   ├── clean-data.R
│   │   ├── check-data.R
│   │   ├── summary-cohort.R
│   │   ├── summary-all.R
│   │   ├── summary-grades.R
│   │   ├── eda-summary.R
│   │   ├── polr-fit.R
│   │   ├── em-fit.R
│   │   ├── hypo-test.R
│   │   └── simulation.R
│   └── plots/
│       └── eda/        # one .R and one .pdf per figure
└── slides/
    └── slides.tex
```

### 0.2 Python environment

The environment is managed with **uv** and is already set up in `pyproject.toml` / `uv.lock`. `pylidc` 0.2.3 is old: it needs `pkg_resources` and the pre-2.0 NumPy/SQLAlchemy APIs, so Python is pinned to 3.11 with `numpy<1.24`, `sqlalchemy<2.0`, `setuptools<80` (and matching pandas/scipy/scikit-image).

```bash
uv sync
uv run python -c "import pylidc as pl; print(pl.query(pl.Scan).count(), pl.query(pl.Annotation).count())"
```

This prints **1018 6859** (CT scans, annotations). The annotation database ships inside the package, so no DICOM download is needed for anything in this project.

### 0.3 R environment

Steps 2–3 (`clean-data.R`, `check-data.R`, `eda-summary.R`) need only base R. Install the rest before Step 4. Some dependencies of `ordinal` compile Fortran, so install a Fortran compiler first:

```bash
sudo pacman -S gcc-fortran
```


```r
install.packages(c("ordinal", "brant",
                   "ggplot2", "corrplot", "xtable"))
```

☐ `pylidc` returns 1018 scans and 6859 annotations; all R packages load.

---

## Step 1 — Background reading (Week 1, in parallel with Step 2)

Do this before writing any modelling code. Each item maps to a part of the report.

| Topic                   | What you must be able to do                                                                                                                                                              | Source                                                                          |
| ----------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------- |
| Proportional odds model | Write the cumulative logit model, explain it as a latent $Y^* = x^\top\beta + \varepsilon$ with logistic $\varepsilon$ and thresholds $\alpha_k$, interpret $e^{\beta}$ as an odds ratio | Agresti, *Categorical Data Analysis*, ch. on ordinal models (Topic 2 notes)     |
| Brant test              | It is a **Wald** test comparing separate binary-logit fits at each threshold; know the alternative LR test of partial proportional odds (`ordinal::nominal_test`)                        | Brant (1990), *Biometrics* 46(4)                                                |
| Wilks' theorem / LRT    | $-2\log\Lambda \to \chi^2_{df}$ under $H_0$; regularity conditions; nested models only                                                                                                   | Casella & Berger, ch. 10                                                        |
| EM algorithm            | Derive the Q-function, prove monotonicity via Jensen/KL, know local-maximum and slow-convergence issues                                                                                  | McLachlan & Krishnan ch. 1–4; Wasserman, *All of Statistics*                    |
| Dawid–Skene             | Latent true class + confusion matrices; E- and M-step updates; label switching                                                                                                           | Dawid & Skene (1979), JRSS-C 28(1) 20–28, doi:10.2307/2346806; arXiv:1812.02736 |

**Warm-up exercise (strongly recommended):** derive and code EM for a 2-component Gaussian mixture on simulated data. It has the same E-step/M-step structure as Dawid–Skene and takes an evening.

☐ One page of handwritten notes per row of the table; Gaussian-mixture EM working.

---

## Step 2 — Data extraction with pylidc (Week 1)

### 2.1 Understand the data structure

- **Scan**: one CT study. Has `patient_id`, `pixel_spacing`, `slice_thickness`, `slice_spacing`.
- **Annotation**: one radiologist's outline of one nodule (≥ 3 mm) in one scan. Has:
  - geometric: `diameter`, `surface_area`, `volume` (computed from the contours, no images needed);
  - semantic ratings (1–5 or 1–6): `subtlety`, `internalStructure`, `calcification`, `sphericity`, `margin`, `lobulation`, `spiculation`, `texture`, **`malignancy`**.
- **Nodule** (not a table): `scan.cluster_annotations()` groups annotations that refer to the same physical nodule. Each cluster has 1–4 annotations. **The nodule is your unit of analysis.**

Hierarchy: patient → scan(s) → nodule(s) → 1–4 reader annotations.

### 2.2 Extraction script (`python/extract.py`)

```bash
uv run python python/extract.py
```

Writes `data/raw/annotations.csv`, the **master file**, with one row per **annotation**. It carries patient, scan and nodule IDs so you can aggregate later. Run it once (it is slow: pylidc's surface-area and volume computations dominate) and never edit the output by hand.

| Group | Columns |
|---|---|
| IDs | `patient_id`, `scan_id`, `nodule_id` (`<scan_id>_<cluster>`), `ann_id`, `reader_slot` (**not** a stable reader identity), `n_readers` |
| Scan geometry | `pixel_spacing`, `slice_thickness`, `slice_spacing` |
| pylidc geometry | `diameter` (mm), `surface_area` (mm²), `volume` (mm³) |
| Mask shape | `n_voxels`, `n_slices`, `bbox_x/y/z` (mm), `elongation`, `flatness` (principal-axis ratios, (0, 1]), `convex_volume` (mm³, convex hull counted in voxels), `solidity` (`n_voxels` / hull voxels, (0, 1]) |
| Ratings | `subtlety`, `internalStructure`, `calcification`, `sphericity`, `margin`, `lobulation`, `spiculation`, `texture`, **`malignancy`** |

Nodules drawn on a single slice are flat in 3D, so their convex hull is taken in 2D.

### 2.3 Sanity checks

- Number of annotations, nodules, scans, patients. Expect 6,859 annotations and roughly 2,600 nodule clusters.
- Distribution of `n_readers` (how many nodules have 1, 2, 3, 4 readers).
- No missing `malignancy`; ratings all within their documented ranges.
- Some patients have more than one scan. Note this for the train/test split.

### 2.4 Important limitation to acknowledge now

Without the DICOM images you **cannot** compute intensity-based radiomic texture features (GLCM, GLRLM, etc.). Your feature pool is:

1. contour-derived geometry (≈ 10–15 features from 2.2), and
2. the radiologists' semantic ratings (8 features).

That is not "dozens of texture features". State this honestly in the report. (Optional extension: download the DICOMs for a subset and run `pyradiomics`. Skip unless you have spare time.)

☐ `data/raw/annotations.csv` exists; a short table of counts for the report's Data section.

---

## Step 3 — Cleaning and building the nodule-level table (Week 2)

```bash
Rscript codes/r-codes/clean-data.R
```

Base R only. Writes to `data/processed/`:

| File | Contents |
|---|---|
| `nodules.csv` | Master nodule table: every column below plus `calcification`, `internalStructure` and `split` |
| `nodules_train.csv`, `nodules_test.csv` | Modelling tables (no nominal ratings, no `split`) |
| `median_ties.csv` | Nodules whose median malignancy fell between two grades |

### 3.1 Decide which features are predictors

The semantic ratings (spiculation, margin, …) come from the **same radiologists** who rate malignancy, so using them as predictors is partly circular. Recommended design:

- **Primary model:** contour-derived geometric features only.
- **Secondary model:** geometric + semantic features (reported as a comparison, with the circularity caveat).

`calcification` and `internalStructure` are **nominal** (category labels, not rankings), and almost every nodule has the same value (no calcification, soft tissue). They are **dropped from modelling**; their per-nodule mode is kept in `nodules.csv` for reference.

### 3.2 Clusters with more than 4 annotations

LIDC had four readers per scan, so a cluster of 5+ annotations means `cluster_annotations()` merged neighbouring lesions. These clusters are dropped; report how many.

### 3.3 Aggregate to one row per nodule

For each `nodule_id`:

- geometric features: **mean** across readers (each reader drew a slightly different contour);
- ordinal semantic features: **median** across readers (half-grades such as 2.5 are kept, since these are numeric predictors);
- malignancy: keep **all** reader ratings in wide form `R1..R4` (with `NA` for missing readers) for the EM model, plus `mal_median` for the naive models.

Also store the count matrix `C` (n × 5) as columns `C1..C5`: `C[i, r]` = number of readers who gave nodule *i* rating *r*. The EM code uses this.

**Median ties.** With an even number of readers (2 or 4) the median falls between two grades whenever the two middle ratings differ, e.g. (3, 4) → 3.5 or (2, 3, 4, 4) → 3.5. `mal_median` rounds these **up** (conservative for triage), `mal_median_raw` keeps the unrounded value and `mal_median_tie` flags them. State the rule in the report. The EM model does not need this rule, since it uses all ratings.

### 3.4 Transformations

- `volume`, `surface_area`, `convex_volume` are heavily right-skewed: `log_volume`, `log_surface_area`, `log_convex_volume` are added.
- Derived: `compactness = volume^(2/3) / surface_area`. (`solidity` already comes from extraction.)
- The CSVs store **unstandardised** values. Standardise (mean 0, sd 1) inside each modelling script using the **training-set** means/SDs, so test-set information never leaks into the saved files.

### 3.5 Train/test split by patient

Split at the **patient** level (70/30, seed `20261105`), not the nodule level, so nodules from one patient never appear in both sets. The assignment is the `split` column of `nodules.csv`.

**All model fitting uses the training set only.** The test set is touched once, in Step 8.

### 3.6 Exploratory analysis

- Bar chart of `mal_median`; bar chart of all individual ratings.
- Boxplots of each geometric feature by `mal_median`.
- Correlation heatmap of predictors (expect diameter, log-volume, log-surface-area to be ≈ 0.9+ correlated).
- **Inter-reader agreement:** proportion of nodules where readers disagree; Fleiss' κ or the distribution of (max − min) rating per nodule. This motivates the EM model.

☐ `data/processed/nodules_train.csv`, `nodules_test.csv`; 3–4 EDA figures in `codes/plots/eda/`.

---

## Step 4 — Naive proportional odds model (Week 3)

Script: `codes/r-codes/polr-fit.R`. Topic 2. This is the baseline the EM model will be compared against.

### 4.1 Fit

```r
library(MASS)
m_naive <- polr(mal_median ~ x1 + x2 + x3, data = train, Hess = TRUE, method = "logistic")
summary(m_naive)
exp(coef(m_naive))           # odds ratios
exp(confint(m_naive))        # profile-likelihood CIs
```

Model: $\operatorname{logit} P(Y \le k \mid x) = \alpha_k - x^\top\beta$, $k = 1,\dots,4$.

### 4.2 Interpret

For a one-SD increase in $x_j$, the odds of being in a **higher** malignancy category (vs. at or below any given threshold) multiply by $e^{\beta_j}$, the same factor at every threshold. Note R's `polr` sign convention ($\alpha_k - x^\top\beta$), so positive $\beta$ means higher malignancy.

### 4.3 Figures

- Predicted category probabilities $P(Y = k \mid x)$ vs. a size variable, with other predictors at their mean: five curves, one per category.
- Odds-ratio forest plot with 95% CIs.

### 4.4 Secondary model

Repeat with semantic features added; report the comparison in one table.

☐ Coefficient table, odds-ratio plot, probability-curve plot.

---

## Step 5 — EM latent true-grade model (Weeks 4–5)

Script: `codes/r-codes/em-fit.R`. This is the core of the extension (Phase 4 of `PLAN.md`).

### 5.1 The model

- Latent true grade $Z_i \in \{1,\dots,5\}$ with $P(Z_i \le k \mid x_i) = \operatorname{logit}^{-1}(\alpha_k - x_i^\top\beta)$.
- Reader ratings: $P(R_{ij} = r \mid Z_i = k) = \pi_{kr}$, one shared $5\times5$ matrix, rows summing to 1.
- Ratings conditionally independent given $Z_i$.
- Observed-data log-likelihood:
  $$\ell(\theta) = \sum_i \log \sum_{k=1}^5 p_{ik}(\alpha,\beta) \prod_{r=1}^5 \pi_{kr}^{C_{ir}}$$
  where $p_{ik} = P(Z_i = k \mid x_i)$ and $C_{ir}$ is the count matrix from Step 3.3.

### 5.2 Write the derivation first (on paper, then into the appendix)

1. Complete-data log-likelihood: $\ell_c = \sum_i \sum_k \mathbb{1}\{Z_i=k\}\left[\log p_{ik} + \sum_r C_{ir}\log\pi_{kr}\right]$.
2. Q-function: replace $\mathbb{1}\{Z_i=k\}$ by $w_{ik} = P(Z_i = k \mid R_i, x_i; \theta^{(t)})$.
3. Q separates into a $\pi$ part (closed-form maximiser) and an $(\alpha,\beta)$ part (a weighted ordinal likelihood).
4. Monotonicity proof: $\ell(\theta) = Q(\theta\mid\theta^{(t)}) - H(\theta\mid\theta^{(t)})$, and $H$ is maximised at $\theta^{(t)}$ by Jensen / Gibbs' inequality.

### 5.3 Implementation

```r
library(MASS)

logsumexp <- function(M) {             # row-wise log-sum-exp
  m <- apply(M, 1, max)
  m + log(rowSums(exp(M - m)))
}

em_dawid_skene_polr <- function(X, C, K = 5, w_init, max_iter = 500, tol = 1e-8) {
  n <- nrow(X); w <- w_init; ll_trace <- numeric(0)
  X <- as.data.frame(X)
  long <- X[rep(seq_len(n), times = K), , drop = FALSE]
  long$z <- factor(rep(seq_len(K), each = n), levels = seq_len(K), ordered = TRUE)

  for (t in seq_len(max_iter)) {
    ## ---- M-step: confusion matrix (closed form) ----
    num <- t(w) %*% C                              # K x K
    Pi  <- num / rowSums(num)
    Pi  <- pmax(Pi, 1e-10); Pi <- Pi / rowSums(Pi)

    ## ---- M-step: weighted proportional-odds fit ----
    long$wt <- as.vector(w)
    keep <- long$wt > 1e-8
    fit <- polr(z ~ . - wt, data = long[keep, ], weights = wt, Hess = FALSE)

    ## ---- E-step ----
    P    <- predict(fit, newdata = X, type = "probs")          # n x K prior
    logW <- log(pmax(P, 1e-300)) + C %*% t(log(Pi))            # n x K
    lse  <- logsumexp(logW)
    w    <- exp(logW - lse)

    ## ---- convergence ----
    ll_trace <- c(ll_trace, sum(lse))              # observed-data log-lik
    if (t > 1 && abs(diff(tail(ll_trace, 2))) < tol * abs(tail(ll_trace, 1))) break
  }
  list(fit = fit, Pi = Pi, w = w, loglik = tail(ll_trace, 1), trace = ll_trace)
}
```

Notes:

- **Initialisation:** `w_init[i, k] = 1` if `k == mal_median[i]` (optionally smoothed, e.g. 0.8 on the median and 0.05 elsewhere). Also run 10–20 random Dirichlet initialisations and keep the fit with the largest log-likelihood.
- **Check the trace is non-decreasing.** A decrease means a bug (usually in the E-step or in how weights are passed to `polr`).
- **Label switching:** after fitting, check that $\hat\pi$ is diagonally dominant ($\hat\pi_{kk}$ is the largest entry in row $k$). If not, the run converged to a permuted solution; discard it.
- **If $\hat\pi$ is unstable** (rows for rare categories estimated from very few nodules): constrain it to be banded ($\pi_{kr} = 0$ for $|k - r| > 2$) or use the one-coin model ($\pi_{kk} = p$, off-diagonals equal).
- `polr` may warn about non-integer weights; that is harmless here. You can also use `ordinal::clm` with `weights =`.

### 5.4 Standard errors

Parametric or nonparametric bootstrap at the **patient** level:

1. Resample patients with replacement (B = 200 is enough for a course project).
2. Re-run EM, initialised at the full-data MLE (fast convergence).
3. Report bootstrap SEs and percentile CIs for $\beta$ and the diagonal of $\pi$.

(Louis' identity is the analytical alternative; mention it in the theory section even if you use the bootstrap.)

### 5.5 Check with Crowd-Kit (optional but quick)

Fit plain Dawid–Skene without covariates (set the structural part to just class proportions) in your code and in Python's `crowd-kit` (`crowdkit.aggregation.DawidSkene`). The posterior class probabilities should agree closely.

### 5.6 Outputs

- Table: $\hat\beta$ (naive, Step 4) vs. $\hat\beta$ (EM), with SEs.
- Heatmap of $\hat\pi$ (rows = true grade, columns = reported rating).
- Log-likelihood trace plot.
- Distribution of posterior entropy of $w_i$: which nodules are "ambiguous".

☐ Converged EM fit; four outputs above.

---

## Step 6 — Simulation study (Week 5)

Script: `codes/r-codes/simulation.R`. This shows your EM implementation is correct, and it is a strong appendix item.

1. Fix "true" values: $\beta$ and $\alpha$ near your EM estimates, and a realistic $\pi$ (e.g. 0.6 on the diagonal, 0.15 on adjacent cells, the rest small).
2. Simulate $X$ by resampling rows of the real design matrix, $Z_i$ from the PO model, and 1–4 ratings per nodule from $\pi$ (matching the real distribution of `n_readers`).
3. Fit the EM model and the naive median-rating model.
4. Repeat 200 times. Report bias and RMSE of $\hat\beta$ for both methods, and coverage of bootstrap CIs if time allows.

Expected finding: the naive model's $\hat\beta$ is attenuated (biased toward 0) because the median rating is a noisy version of $Z$; EM is approximately unbiased. Whether or not that happens, report it.

☐ Simulation table (bias / RMSE, naive vs EM).

---

## Step 7 — Diagnostics and hypothesis tests (Week 6)

Script: `codes/r-codes/hypo-test.R`. Phase 5 of `PLAN.md`.

### 7.1 Proportional odds assumption

- **Brant test** (a Wald test) on the naive model:

  ```r
  library(brant)
  brant(m_naive)          # omnibus + per-variable chi-square tests
  ```

- **Likelihood-ratio alternative:** fit with `ordinal::clm` and run `nominal_test(model)`, which compares the PO model with a partial PO model that frees each predictor's coefficient across thresholds.
- If a variable violates PO, either report it as a limitation or fit a partial PO model for that variable (`clm(..., nominal = ~ x_bad)`).
- Note: these tests are applied to the naive model; in the EM model PO is assumed for the latent grade.

### 7.2 Likelihood-ratio tests (Wilks)

For nested models $M_0 \subset M_1$: $\Lambda = -2(\ell_0 - \ell_1) \xrightarrow{d} \chi^2_{df}$, $df$ = difference in number of parameters.

- **Joint significance of all predictors:** full model vs. intercept-only (thresholds only).
- **Each predictor:** drop one at a time (`drop1(m_naive, test = "Chisq")`).
- **Semantic features add anything?** Primary vs. secondary model.
- **For the EM model:** refit EM under the reduced model and use the **observed-data** log-likelihoods from the `loglik` output. Do not use the complete-data likelihood.

Write a short paragraph linking this to Neyman–Pearson: the LRT is the most powerful test in the simple-vs-simple case, and Wilks gives its asymptotic null distribution in the composite case.

### 7.3 Model fit

- For the EM model, compare against the plain Dawid–Skene model without covariates (LRT: are the features informative about the true grade beyond the readers' ratings?).

☐ Brant/nominal-test table; LRT table.

---

## Step 8 — Evaluation on the test set (Week 6)

Touch the test set **once**, after all modelling decisions are frozen.

1. Apply the training-set standardisation to the test predictors.
2. For each model (naive primary, naive secondary, EM), compute:
   - accuracy and mean absolute error of the predicted category (mode of $P(Z = k \mid x)$) vs. `mal_median`;
   - quadratic-weighted Cohen's κ;
   - AUC for the binary triage rule "malignancy ≥ 4" using $P(Z \ge 4 \mid x)$;
   - test log-likelihood of the observed **reader ratings** (for the EM model this is the observed-data log-likelihood; it is the fairest comparison).
3. Calibration plot for $P(Z \ge 4 \mid x)$ (binned predicted vs. observed proportion).

☐ Test-set performance table; calibration plot.

---

## Step 9 — Policy interpretation (Week 7)

Tie the results to the triage recommendation in `PLAN.md`:

- Propose a rule such as: flag a nodule for priority review if $P(Z \ge 4 \mid x) > \tau$. Choose $\tau$ from the test-set ROC (e.g. sensitivity ≥ 0.9) and report the resulting workload (fraction flagged).
- Use $\hat\pi$ to say something concrete about reader reliability (e.g. "a true grade-4 nodule is rated 3 or lower X% of the time"). This supports using a model as a second reader.

**Limitations you must state:**

1. LIDC "malignancy" is radiologists' **suspicion**, not pathology-confirmed diagnosis. The model predicts expert opinion, not cancer. (TCIA hosts a diagnosis spreadsheet for a small subset of LIDC patients; mention it as possible external validation.)
2. No intensity-based radiomics (no DICOM).
3. Circularity of semantic features.
4. Conditional-independence assumption of readers in the EM model.
5. LIDC is a research dataset from a few US centres; generalisation to a public-health setting such as Indian hospitals is untested.

---

## Step 10 — Write the report (Weeks 7–8)

`main.tex`. **Max 5 pages, double-spaced, A4.** Appendix does not count.

### 10.1 Preamble

```latex
\documentclass[11pt]{article}
\usepackage[a4paper, margin=2.5cm]{geometry}
\usepackage{setspace}\doublespacing
\usepackage{amsmath, amssymb, booktabs, graphicx, hyperref}
\usepackage[style=authoryear]{biblatex}\addbibresource{refs.bib}
```

### 10.2 Main text (≈ page budget)

1. **Introduction and research questions** (0.5 p). Example RQs:
   - RQ1: Which geometric nodule features are associated with radiologist-assessed malignancy?
   - RQ2: Does the proportional-odds model adequately describe that relationship?
   - RQ3: Does accounting for inter-reader disagreement (EM) change the estimated effects?
2. **Data** (0.75 p): LIDC-IDRI, pylidc, unit of analysis, aggregation, train/test split, counts table.
3. **Methodology** (1.5 p): PO model → EM latent-grade model → Brant/LRT. Equations, not code.
4. **Results** (1.5 p): one coefficient table (naive vs EM), one key figure (probability curves or $\hat\pi$ heatmap), test results in text.
5. **Discussion, policy recommendations, limitations** (0.75 p).

### 10.3 Appendix (no page limit)

- A: Full code (use `listings` or `minted`, or `\input` the scripts).
- B: Full coefficient tables, odds-ratio plot, probability curves for each category.
- C: EM derivation, monotonicity proof, log-likelihood trace.
- D: Simulation study.
- E: Brant and LRT tables; calibration plot.

### 10.4 Bibliography (minimum)

Armato et al. (2011, *Medical Physics*) for LIDC-IDRI; Hancock & Magnan (2016) for pylidc; McCullagh (1980) for PO models; Brant (1990); Dempster, Laird & Rubin (1977); Dawid & Skene (1979); Agresti.

Compile, check the page count with the appendix excluded, and ask someone outside the group to read it once.

☐ `report.pdf`, ≤ 5 pages of main text.

---

## Step 11 — Slides and defense (Week 8)

`slides/slides.tex` (Beamer). **Fewer than 16 slides, 20–25 minutes, 5–7 minutes per member.**

Suggested 13-slide structure and split (assuming 4 members):

| #   | Slide                                                | Speaker |
| --- | ---------------------------------------------------- | ------- |
| 1   | Title                                                | A       |
| 2   | Motivation: nodule triage, reader disagreement       | A       |
| 3   | Data: LIDC-IDRI, unit of analysis, counts            | A       |
| 4   | EDA: malignancy distribution, reader disagreement    | B       |
| 5   | Proportional odds model                              | B       |
| 6   | Naive model results: odds ratios, probability curves | B       |
| 7   | Why EM: latent true grade                            | C       |
| 8   | E-step / M-step and monotonicity                     | C       |
| 9   | EM results: $\hat\pi$ heatmap, $\beta$ comparison    | C       |
| 10  | Simulation study                                     | C       |
| 11  | Brant test and LRTs                                  | D       |
| 12  | Test-set performance and triage rule                 | D       |
| 13  | Limitations and conclusions                          | D       |

Preparation:

- Every member should be able to answer questions on every slide, not only their own.
- Likely questions: Why not plain logistic regression? Why is Brant a Wald test? Why does EM increase the likelihood? Are readers really conditionally independent? What does $\hat\pi$ mean clinically? Why is malignancy not ground truth?
- Do two full timed rehearsals.

☐ `slides.pdf`; two rehearsals done.

---

## Suggested timeline (8 weeks; compress if your deadline is sooner)

| Week | Steps      | Milestone                                                       |
| ---- | ---------- | --------------------------------------------------------------- |
| 1    | 0, 1, 2    | Environment works; `annotations.csv` extracted; reading started |
| 2    | 3          | Nodule-level train/test tables; EDA                             |
| 3    | 4          | Naive PO model fitted                                           |
| 4    | 5.1–5.3    | EM derivation written; EM code converging                       |
| 5    | 5.4–5.6, 6 | Bootstrap SEs; simulation study                                 |
| 6    | 7, 8       | Diagnostics, LRTs, test-set evaluation                          |
| 7    | 9, 10      | Policy section; report draft                                    |
| 8    | 10, 11     | Report final; slides; rehearsals                                |

**Fallback if time runs short:** the EM model (Steps 5–6) is the extension. Steps 2–4 and 7–11 alone form a complete project that matches the original plan. Get those working first; do not let EM debugging block the rest.
