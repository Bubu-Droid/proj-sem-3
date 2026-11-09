# Likelihood ratio tests for the EM model (STEPS.md Step 6), using observed-data
# log-likelihoods. Each reduced model is refitted from two starts: the median
# start and the full model's posterior weights.
#
# Run from the repository root:
#   Rscript codes/r-codes/hypo-test-em.R

library(MASS)
library(parallel)

data <- read.csv("data/processed/nodules_train.csv", stringsAsFactors = FALSE)
trimmed <- read.csv("data/processed/nodules_train_trimmed.csv", stringsAsFactors = FALSE)

PRIMARY <- c("log_volume", "elongation", "flatness", "solidity", "compactness")
SECONDARY <- c("spiculation", "lobulation", "margin", "texture", "sphericity")
K <- 5
EPS <- 1e-6
MAX_ITER <- 3000
N_CORES <- 12

# Same nodules as em-fit.R
keep <- complete.cases(data[c(PRIMARY, SECONDARY)])
data <- data[keep, ]
C <- as.matrix(trimmed[keep, paste0("C", 1:K)])
X_all <- scale(data[PRIMARY])
n <- nrow(X_all)

logsumexp <- function(M) {
  m <- apply(M, 1, max)
  m + log(rowSums(exp(M - m)))
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

em_fit <- function(vars, w) {
  X <- X_all[, vars, drop = FALSE]
  long <- as.data.frame(X[rep(seq_len(n), times = K), , drop = FALSE])
  long$z <- factor(rep(seq_len(K), each = n), levels = seq_len(K), ordered = TRUE)
  form <- reformulate(vars, response = "z")
  ctrl <- list(reltol = 1e-10, maxit = 500)
  start <- NULL
  last <- -Inf
  for (t in seq_len(MAX_ITER)) {
    Pi <- update_pi(w)
    long$wt <- as.vector(w)
    use <- long$wt > 1e-8
    fit <- suppressWarnings(
      if (is.null(start))
        polr(form, data = long[use, ], weights = wt, method = "logistic",
             control = ctrl)
      else
        polr(form, data = long[use, ], weights = wt, start = start,
             method = "logistic", control = ctrl))
    start <- c(coef(fit), fit$zeta)
    eta <- drop(X %*% coef(fit))
    cum <- cbind(0, sapply(fit$zeta, function(a) plogis(a - eta)), 1)
    e <- em_step(log(pmax(cum[, -1] - cum[, -(K + 1)], 1e-300)), Pi)
    w <- e$w
    if (t > 1 && abs(e$loglik - last) < EPS) break
    last <- e$loglik
  }
  list(loglik = e$loglik, iter = t, converged = t < MAX_ITER, w = w)
}

em_fit_nocov <- function(w) {
  last <- -Inf
  for (t in seq_len(MAX_ITER)) {
    Pi <- update_pi(w)
    p <- colMeans(w)
    e <- em_step(matrix(log(pmax(p, 1e-300)), n, K, byrow = TRUE), Pi)
    w <- e$w
    if (t > 1 && abs(e$loglik - last) < EPS) break
    last <- e$loglik
  }
  list(loglik = e$loglik, iter = t, converged = t < MAX_ITER)
}

w_median <- matrix(0.05, n, K)
w_median[cbind(seq_len(n), data$mal_median)] <- 0.8

# ======================================
#   FULL AND NO-COVARIATE MODELS
# ======================================

full <- em_fit(PRIMARY, w_median)
base <- em_fit_nocov(w_median)
cat(sprintf("Full model: loglik %.3f, iterations %d, converged %s\n",
            full$loglik, full$iter, full$converged))
cat(sprintf("No covariates: loglik %.3f, iterations %d, converged %s\n\n",
            base$loglik, base$iter, base$converged))

report <- function(label, ll_reduced, df) {
  lr <- 2 * (full$loglik - ll_reduced)
  cat(sprintf("%s: reduced loglik %.3f, LR %.3f, df %d, p %.3g\n", label,
              ll_reduced, lr, df, pchisq(lr, df, lower.tail = FALSE)))
}

report("Features vs no features", base$loglik, length(PRIMARY))

# ======================================
#   DROP-ONE TESTS
# ======================================

jobs <- expand.grid(v = PRIMARY, s = 1:2, stringsAsFactors = FALSE)
fits <- mclapply(seq_len(nrow(jobs)), function(j)
  em_fit(setdiff(PRIMARY, jobs$v[j]), if (jobs$s[j] == 1) w_median else full$w),
  mc.cores = N_CORES, mc.preschedule = FALSE)

cat("\nDrop-one tests\n")
for (v in PRIMARY) {
  f <- fits[jobs$v == v]
  lls <- sapply(f, `[[`, "loglik")
  report(v, max(lls), 1)
  cat(sprintf("  starts agree to %.4f; converged %s\n", abs(diff(lls)),
              all(sapply(f, `[[`, "converged"))))
}
