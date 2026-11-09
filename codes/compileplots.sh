#!/usr/bin/env bash

set -e
cd "$(dirname "$0")/.."

Rscript codes/r-codes/clean-data.R
Rscript codes/r-codes/check-data.R

Rscript codes/r-codes/summary-cohort.R
Rscript codes/r-codes/summary-all.R
Rscript codes/r-codes/summary-grades.R
Rscript codes/r-codes/eda-summary.R
Rscript codes/r-codes/polr-fit.R
Rscript codes/r-codes/em-fit.R
Rscript codes/r-codes/hypo-test.R
Rscript codes/r-codes/hypo-test-em.R
Rscript codes/r-codes/test-eval.R

Rscript codes/plots/data-description/mal-distribution.R
Rscript codes/plots/data-description/features-by-grade.R
Rscript codes/plots/data-description/correlation.R
Rscript codes/plots/data-description/reader-agreement.R

Rscript codes/plots/model-formulation/probability-curves.R
Rscript codes/plots/model-formulation/odds-ratio-forest.R
Rscript codes/plots/model-formulation/em-convergence.R
Rscript codes/plots/model-formulation/em-rating-probs.R

Rscript codes/plots/model-analysis/threshold-odds-ratios.R

Rscript codes/plots/model-testing/predicted-grade-counts.R
Rscript codes/plots/model-testing/triage-curve.R
