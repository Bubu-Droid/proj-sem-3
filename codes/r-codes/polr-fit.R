# Naive proportional odds models for median malignancy (STEPS.md Step 4).
# Predictors are standardised with training-set means and SDs.
#
# Run from the repository root:
#   Rscript codes/r-codes/polr-fit.R

library(MASS)

data <- read.csv("data/processed/nodules_train.csv", stringsAsFactors = FALSE)

PRIMARY <- c("log_volume", "elongation", "flatness", "solidity", "compactness")
SECONDARY <- c("spiculation", "lobulation", "margin", "texture", "sphericity")

# Same nodules in both models, so that they can be compared later
data <- data[complete.cases(data[c(PRIMARY, SECONDARY)]), ]
cat("Nodules:", nrow(data), "\n\n")

data[c(PRIMARY, SECONDARY)] <- scale(data[c(PRIMARY, SECONDARY)])
data$mal_median <- factor(data$mal_median, levels = 1:5, ordered = TRUE)

report <- function(vars, label) {
  f <- reformulate(vars, response = "mal_median")
  fit <- polr(f, data = data, Hess = TRUE, method = "logistic")
  est <- coef(summary(fit))
  ci <- suppressMessages(confint(fit))

  cat("======================================\n")
  cat(" ", label, "\n")
  cat("======================================\n\n")

  cat("Coefficients: beta & SE & OR & OR 2.5% & OR 97.5%\n")
  for (v in vars)
    cat(sprintf("%s: %.3f & %.3f & %.3f & %.3f & %.3f\n", v, est[v, "Value"],
                est[v, "Std. Error"], exp(est[v, "Value"]),
                exp(ci[v, 1]), exp(ci[v, 2])))

  cat("\nThresholds: alpha & SE\n")
  for (a in names(fit$zeta))
    cat(sprintf("%s: %.3f & %.3f\n", a, est[a, "Value"], est[a, "Std. Error"]))

  cat(sprintf("\nConverged: %s\n", fit$convergence == 0))
  cat(sprintf("Log-likelihood: %.3f\n\n", as.numeric(logLik(fit))))
  invisible(fit)
}

report(PRIMARY, "PRIMARY MODEL")
report(c(PRIMARY, SECONDARY), "SECONDARY MODEL")
