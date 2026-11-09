# Hypothesis tests for the naive proportional odds models (STEPS.md Step 6).
# Predictors are standardised with training-set means and SDs.
#
# Run from the repository root:
#   Rscript codes/r-codes/hypo-test.R

library(MASS)
library(ordinal)

data <- read.csv("data/processed/nodules_train_trimmed.csv", stringsAsFactors = FALSE)

PRIMARY <- c("log_volume", "elongation", "flatness", "solidity", "compactness")
SECONDARY <- c("spiculation", "lobulation", "margin", "texture", "sphericity")

# Same nodules as polr-fit.R
data <- data[complete.cases(data[c(PRIMARY, SECONDARY, "mal_median")]), ]
cat("Nodules:", nrow(data), "\n\n")

data[c(PRIMARY, SECONDARY)] <- scale(data[c(PRIMARY, SECONDARY)])
data$y <- factor(data$mal_median, levels = 1:5, ordered = TRUE)

fit_polr <- function(vars) {
  f <- if (length(vars)) reformulate(vars, response = "y") else y ~ 1
  polr(f, data = data, Hess = TRUE, method = "logistic")
}

null_fit <- fit_polr(character(0))
primary <- fit_polr(PRIMARY)
secondary <- fit_polr(c(PRIMARY, SECONDARY))

lr_line <- function(label, small, big) {
  lr <- 2 * as.numeric(logLik(big) - logLik(small))
  df <- length(coef(big)) - length(coef(small))
  cat(sprintf("%s: loglik %.3f vs %.3f, LR %.3f, df %d, p %.3g\n", label,
              as.numeric(logLik(big)), as.numeric(logLik(small)), lr, df,
              pchisq(lr, df, lower.tail = FALSE)))
}

# ======================================
#   JOINT LIKELIHOOD RATIO TESTS
# ======================================

lr_line("Primary vs thresholds only", null_fit, primary)
lr_line("Secondary vs thresholds only", null_fit, secondary)
lr_line("Secondary vs primary", primary, secondary)

# ======================================
#   DROP-ONE LIKELIHOOD RATIO TESTS
# ======================================

drop_table <- function(fit, vars, label) {
  cat("\nDrop-one tests,", label, "\n")
  for (v in vars) {
    reduced <- fit_polr(setdiff(vars, v))
    lr_line(v, reduced, fit)
  }
}

drop_table(primary, PRIMARY, "primary model")
drop_table(secondary, c(PRIMARY, SECONDARY), "secondary model")

# ======================================
#   PROPORTIONAL ODDS ASSUMPTION
# ======================================

po_tests <- function(vars, label) {
  cat("\n======================================\n")
  cat(" ", label, "\n")
  cat("======================================\n\n")

  # clm() re-evaluates its call, so the formulas go into it directly
  f <- reformulate(vars, response = "y")
  nom <- reformulate(vars)
  po <- eval(bquote(clm(.(f), data = data, link = "logit")))
  free <- eval(bquote(clm(y ~ 1, nominal = .(nom), data = data, link = "logit")))

  lr <- 2 * as.numeric(logLik(free) - logLik(po))
  df <- length(coef(free)) - length(coef(po))
  cat(sprintf("All predictors: loglik %.3f vs %.3f, LR %.3f, df %d, p %.3g\n\n",
              as.numeric(logLik(free)), as.numeric(logLik(po)), lr, df,
              pchisq(lr, df, lower.tail = FALSE)))

  cat("One predictor free: LR, df, p\n")
  nt <- nominal_test(po)
  for (v in vars)
    cat(sprintf("%s: %.3f, %d, %.3g\n", v, nt[v, "LRT"], nt[v, "Df"], nt[v, "Pr(>Chi)"]))

  # clm's coefficients act on P(Y <= k), so the sign is flipped for higher grades
  cat("\nOdds ratio of a higher grade at each threshold (all predictors free):\n")
  b <- coef(free)
  for (v in vars)
    cat(sprintf("%s: %s\n", v, paste(sprintf("%.3f", exp(-b[paste0(1:4, "|", 2:5, ".", v)])),
                                     collapse = " & ")))
}

po_tests(PRIMARY, "PRIMARY MODEL")
po_tests(c(PRIMARY, SECONDARY), "SECONDARY MODEL")
