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

Rscript codes/plots/eda/mal-distribution.R
Rscript codes/plots/eda/features-by-grade.R
Rscript codes/plots/eda/correlation.R
Rscript codes/plots/eda/reader-agreement.R
