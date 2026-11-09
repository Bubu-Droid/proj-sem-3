# EM fit of the latent true-grade model (STEPS.md Step 5) on the training nodules.
# Predictors are the primary features, standardised with training-set means and SDs.
# Ratings come from nodules_train_trimmed.csv (see clean-data.R).
#
# Run from the repository root:
#   Rscript codes/r-codes/em-fit.R

library(MASS)
library(parallel)

data <- read.csv("data/processed/nodules_train.csv", stringsAsFactors = FALSE)
trimmed <- read.csv("data/processed/nodules_train_trimmed.csv", stringsAsFactors = FALSE)

PRIMARY <- c("log_volume", "elongation", "flatness", "solidity", "compactness")
SECONDARY <- c("spiculation", "lobulation", "margin", "texture", "sphericity")
K <- 5
EPS <- 1e-6
MAX_ITER <- 3000
N_RANDOM <- 5
N_CORES <- 3
set.seed(20261105)

# Same nodules as polr-fit.R
keep <- complete.cases(data[c(PRIMARY, SECONDARY)])
data <- data[keep, ]
trimmed <- trimmed[keep, ]
X <- scale(data[PRIMARY])
C_all <- as.matrix(data[paste0("C", 1:K)])
n <- nrow(X)
cat("Nodules:", n, "\n")

# Ratings used for fitting come from the trimmed table
C <- as.matrix(trimmed[paste0("C", 1:K)])
cat("Ratings in the fit:", sum(C), "of", sum(C_all), "\n\n")

# One row per nodule and grade, weighted by w in the M-step
long <- as.data.frame(X[rep(seq_len(n), times = K), ])
long$z <- factor(rep(seq_len(K), each = n), levels = seq_len(K), ordered = TRUE)
form <- reformulate(PRIMARY, response = "z")

# Row-wise log-sum-exp
logsumexp <- function(M) {
  m <- apply(M, 1, max)
  m + log(rowSums(exp(M - m)))
}

# P(Z_i = k | x_i) under the proportional odds model
grade_probs <- function(alpha, beta) {
  eta <- drop(X %*% beta)
  cum <- cbind(0, sapply(alpha, function(a) plogis(a - eta)), 1)
  cum[, -1] - cum[, -(K + 1)]
}

update_pi <- function(w) {
  num <- t(w) %*% C
  Pi <- pmax(num / rowSums(num), 1e-10)
  Pi / rowSums(Pi)
}

em_step <- function(logprior, Pi) {
  logW <- logprior + C %*% t(log(Pi))
  lse <- logsumexp(logW)
  list(w = exp(logW - lse), loglik = sum(lse))
}

# Latent-grade model with covariates
em_fit <- function(w) {
  start <- NULL
  trace <- numeric(0)
  for (t in seq_len(MAX_ITER)) {
    # M-step
    Pi <- update_pi(w)
    long$wt <- as.vector(w)
    keep <- long$wt > 1e-8
    ctrl <- list(reltol = 1e-10, maxit = 500)
    fit <- suppressWarnings(
      if (is.null(start))
        polr(form, data = long[keep, ], weights = wt, method = "logistic",
             control = ctrl)
      else
        polr(form, data = long[keep, ], weights = wt, start = start,
             method = "logistic", control = ctrl))
    beta <- coef(fit); alpha <- fit$zeta
    start <- c(beta, alpha)

    # E-step
    e <- em_step(log(pmax(grade_probs(alpha, beta), 1e-300)), Pi)
    w <- e$w
    trace <- c(trace, e$loglik)
    if (t > 1 && abs(trace[t] - trace[t - 1]) < EPS) break
  }
  list(beta = beta, alpha = alpha, Pi = Pi, w = w, loglik = e$loglik,
       iter = t, trace = trace, converged = t < MAX_ITER)
}

# Same model with beta = 0: one set of grade probabilities for every nodule
em_fit_nocov <- function(w) {
  trace <- numeric(0)
  for (t in seq_len(MAX_ITER)) {
    Pi <- update_pi(w)
    p <- colMeans(w)
    e <- em_step(matrix(log(pmax(p, 1e-300)), n, K, byrow = TRUE), Pi)
    w <- e$w
    trace <- c(trace, e$loglik)
    if (t > 1 && abs(trace[t] - trace[t - 1]) < EPS) break
  }
  list(p = p, Pi = Pi, loglik = e$loglik, iter = t, trace = trace,
       converged = t < MAX_ITER)
}

diag_dominant <- function(Pi) all(max.col(Pi) == seq_len(K))

# ======================================
#   STARTING VALUES
# ======================================

# Median start: 0.8 on the median grade, 0.05 elsewhere
w_median <- matrix(0.05, n, K)
w_median[cbind(seq_len(n), data$mal_median)] <- 0.8

# Random starts: the median start mixed with a random Dirichlet(1) draw
w_random <- function() {
  g <- matrix(rexp(n * K), n, K)
  0.5 * w_median + 0.5 * g / rowSums(g)
}

# ======================================
#   FIT WITH COVARIATES
# ======================================

starts <- c(list(w_median), replicate(N_RANDOM, w_random(), simplify = FALSE))
fits <- mclapply(starts, em_fit, mc.cores = N_CORES, mc.preschedule = FALSE)

for (s in seq_along(fits)) {
  f <- fits[[s]]
  cat(sprintf("start %2d: loglik %.4f, iterations %d, converged %s, min step %.2e\n",
              s - 1, f$loglik, f$iter, f$converged, min(diff(f$trace))))
}
best <- fits[[which.max(sapply(fits, `[[`, "loglik"))]]

cat("\n======================================\n")
cat("  EM MODEL WITH COVARIATES\n")
cat("======================================\n\n")

cat("Coefficients: beta & OR\n")
for (v in PRIMARY)
  cat(sprintf("%s: %.3f & %.3f\n", v, best$beta[v], exp(best$beta[v])))

cat("\nThresholds:", paste(sprintf("%.3f", best$alpha), collapse = ", "), "\n")
cat(sprintf("Observed-data log-likelihood: %.3f\n", best$loglik))
cat(sprintf("Iterations: %d, converged: %s\n\n", best$iter, best$converged))

cat("Confusion matrix (rows = true grade, columns = rating):\n")
dimnames(best$Pi) <- list(1:K, 1:K)
print(round(best$Pi, 3))
cat("Likeliest rating per grade:", max.col(best$Pi), "\n")

cat("\nMean posterior grade probability (weights):",
    paste(sprintf("%.3f", colMeans(best$w)), collapse = " & "), "\n")
cat("Modal posterior grade counts:",
    paste(tabulate(max.col(best$w), K), collapse = " & "), "\n")
cat("Median grade counts:        ",
    paste(tabulate(data$mal_median, K), collapse = " & "), "\n")
cat(sprintf("Modal posterior grade equals median grade: %.1f%%\n",
            100 * mean(max.col(best$w) == data$mal_median)))

cl <- max.col(best$w)
cat("Mean rating (all ratings) per modal grade:",
    paste(sprintf("%.2f", sapply(seq_len(K), function(k)
      sum(sweep(C_all[cl == k, , drop = FALSE], 2, seq_len(K), "*")) /
        sum(C_all[cl == k, ]))), collapse = " & "), "\n")

# ======================================
#   FIT WITHOUT COVARIATES
# ======================================

base <- em_fit_nocov(w_median)

cat("\n======================================\n")
cat("  EM MODEL WITHOUT COVARIATES\n")
cat("======================================\n\n")
cat("Grade probabilities:", paste(sprintf("%.3f", base$p), collapse = " & "), "\n")
cat(sprintf("Observed-data log-likelihood: %.3f\n", base$loglik))
cat(sprintf("Iterations: %d, converged: %s, diagonal dominant: %s\n",
            base$iter, base$converged, diag_dominant(base$Pi)))
