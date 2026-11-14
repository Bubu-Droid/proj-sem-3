# Test-set errors and triage rule of the frozen models (STEPS.md Steps 7-8). Models are fitted on the
# trimmed training table and scored on nodules_test.csv, which keeps all ratings.
#
# Run from the repository root:
#   Rscript codes/r-codes/test-eval.R

library(MASS)

train_trim <- read.csv("data/processed/nodules_train_trimmed.csv", stringsAsFactors = FALSE)
train_full <- read.csv("data/processed/nodules_train.csv", stringsAsFactors = FALSE)
test <- read.csv("data/processed/nodules_test.csv", stringsAsFactors = FALSE)

PRIMARY <- c("log_volume", "elongation", "flatness", "solidity", "compactness")
SECONDARY <- c("spiculation", "lobulation", "margin", "texture", "sphericity")
ALL <- c(PRIMARY, SECONDARY)
K <- 5
EPS <- 1e-6
MAX_ITER <- 3000

test <- test[complete.cases(test[ALL]), ]
y <- test$mal_median
cat("Test nodules:", nrow(test), " patients:", length(unique(test$patient_id)), "\n")
cat("Test grade counts:", paste(tabulate(y, K), collapse = " & "), "\n\n")

errors <- function(pred) {
  c(misclassification = mean(pred != y), MAE = mean(abs(pred - y)),
    RMSE = sqrt(mean((pred - y)^2)))
}

# ======================================
#   NAIVE MODELS
# ======================================

naive <- train_trim[complete.cases(train_trim[c(ALL, "mal_median")]), ]
center <- colMeans(naive[ALL])
spread <- apply(naive[ALL], 2, sd)
naive[ALL] <- scale(naive[ALL], center, spread)
naive$y <- factor(naive$mal_median, levels = 1:K, ordered = TRUE)
test_naive <- test
test_naive[ALL] <- scale(test[ALL], center, spread)

fit_primary <- polr(reformulate(PRIMARY, "y"), data = naive, method = "logistic")
fit_secondary <- polr(reformulate(ALL, "y"), data = naive, method = "logistic")
cat(sprintf("Naive log-likelihoods: %.3f, %.3f\n",
            as.numeric(logLik(fit_primary)), as.numeric(logLik(fit_secondary))))

pred_primary <- max.col(predict(fit_primary, test_naive, type = "probs"), "first")
pred_secondary <- max.col(predict(fit_secondary, test_naive, type = "probs"), "first")

# ======================================
#   EM MODEL
# ======================================

keep <- complete.cases(train_full[ALL])
X <- scale(train_full[keep, PRIMARY])
C <- as.matrix(train_trim[keep, paste0("C", 1:K)])
n <- nrow(X)
long <- as.data.frame(X[rep(seq_len(n), times = K), ])
long$z <- factor(rep(seq_len(K), each = n), levels = seq_len(K), ordered = TRUE)
form <- reformulate(PRIMARY, response = "z")

grade_probs <- function(M, alpha, beta) {
  eta <- drop(M %*% beta)
  cum <- cbind(0, sapply(alpha, function(a) plogis(a - eta)), 1)
  cum[, -1] - cum[, -(K + 1)]
}

logsumexp <- function(M) {
  m <- apply(M, 1, max)
  m + log(rowSums(exp(M - m)))
}

w <- matrix(0.05, n, K)
w[cbind(seq_len(n), train_full$mal_median[keep])] <- 0.8
start <- NULL
last <- -Inf
ctrl <- list(reltol = 1e-10, maxit = 500)
for (t in seq_len(MAX_ITER)) {
  num <- t(w) %*% C
  Pi <- pmax(num / rowSums(num), 1e-10)
  Pi <- Pi / rowSums(Pi)
  long$wt <- as.vector(w)
  use <- long$wt > 1e-8
  fit <- suppressWarnings(
    if (is.null(start))
      polr(form, data = long[use, ], weights = wt, method = "logistic", control = ctrl)
    else
      polr(form, data = long[use, ], weights = wt, start = start, method = "logistic",
           control = ctrl))
  beta <- coef(fit); alpha <- fit$zeta
  start <- c(beta, alpha)
  logW <- log(pmax(grade_probs(X, alpha, beta), 1e-300)) + C %*% t(log(Pi))
  lse <- logsumexp(logW)
  w <- exp(logW - lse)
  ll <- sum(lse)
  if (t > 1 && abs(ll - last) < EPS) break
  last <- ll
}
cat(sprintf("EM log-likelihood: %.3f, iterations %d\n\n", ll, t))

# Standardize with the EM fit's means and SDs, then convert to the rating scale
X_test <- scale(test[PRIMARY], attr(X, "scaled:center"), attr(X, "scaled:scale"))
rating_probs <- grade_probs(X_test, alpha, beta) %*% Pi
pred_em <- max.col(rating_probs, "first")

# ======================================
#   BASELINE AND ERRORS
# ======================================

mode_grade <- which.max(tabulate(naive$mal_median, K))
cat("Baseline grade:", mode_grade, "\n\n")

preds <- list(Baseline = rep(mode_grade, length(y)), `Naive primary` = pred_primary,
              `Naive secondary` = pred_secondary, EM = pred_em)
tab <- t(sapply(preds, errors))
cat("Errors:\n"); print(round(tab, 4))

cat("\nReduction relative to the baseline (%):\n")
print(round(100 * sweep(-sweep(tab, 2, tab["Baseline", ]), 2, tab["Baseline", ], "/"), 2))

cat("\nPredicted grade counts:\n")
print(t(sapply(preds, tabulate, nbins = K)))

# ======================================
#   READER BENCHMARK
# ======================================

# Each reader and each model is scored against the median of the other readers
# on the same nodule, with half-grades rounded up as in mal_median.
R <- as.matrix(test[paste0("R", 1:4)])
bench <- do.call(rbind, lapply(seq_len(nrow(R)), function(i) {
  obs <- which(!is.na(R[i, ]))
  if (length(obs) < 2) return(NULL)
  do.call(rbind, lapply(obs, function(j)
    data.frame(i = i, rating = R[i, j],
               target = ceiling(median(R[i, setdiff(obs, j)])))))
}))
cat("\nMulti-reader test nodules:", length(unique(bench$i)), " rating comparisons:",
    nrow(bench), "\n")
cat("Readers disagreeing on a nodule (%):",
    round(100 * mean(tapply(bench$rating, bench$i, function(r) length(unique(r)) > 1)), 1), "\n")
mae_target <- function(pred) mean(abs(pred[bench$i] - bench$target))
cat("MAE against the median of the other readers:\n")
cat(sprintf("  Single reader: %.4f\n", mean(abs(bench$rating - bench$target))))
for (m in names(preds)) cat(sprintf("  %s: %.4f\n", m, mae_target(preds[[m]])))

# ======================================
#   ERRORS BY CALCIFICATION (EXPLORATORY)
# ======================================

all <- read.csv("data/processed/nodules.csv", stringsAsFactors = FALSE)
calc <- all$calcification[match(test$nodule_id, all$nodule_id)] != 6
cat("\nCalcified test nodules:", sum(calc), "of", length(calc), "\n")
cat("Their mean median grade:", round(mean(y[calc]), 2), " others:", round(mean(y[!calc]), 2), "\n")
for (m in names(preds)) {
  e <- abs(preds[[m]] - y)
  cat(sprintf("%s: MAE calcified %.4f, not calcified %.4f; misclassification calcified %.4f, not calcified %.4f; predicted grade 1 among calcified %d\n",
              m, mean(e[calc]), mean(e[!calc]), mean(e[calc] > 0), mean(e[!calc] > 0),
              sum(preds[[m]][calc] == 1)))
}

# ======================================
#   TRIAGE RULE
# ======================================

# Flag a nodule if its predicted probability of grade >= 4 exceeds tau. The
# threshold is the largest one reaching 90% sensitivity on the training set.
high <- function(P) rowSums(P[, 4:5, drop = FALSE])
choose_tau <- function(p, truth, target = 0.9) {
  p_pos <- sort(p[truth], decreasing = TRUE)
  p_pos[ceiling(target * length(p_pos))]
}
truth_train <- naive$mal_median >= 4
truth_test <- y >= 4

em_keep <- !is.na(train_trim$mal_median[keep])
train_prob <- list(
  `Naive primary` = high(predict(fit_primary, naive, type = "probs")),
  `Naive secondary` = high(predict(fit_secondary, naive, type = "probs")),
  EM = high((grade_probs(X, alpha, beta) %*% Pi)[em_keep, ]))
test_prob <- list(
  `Naive primary` = high(predict(fit_primary, test_naive, type = "probs")),
  `Naive secondary` = high(predict(fit_secondary, test_naive, type = "probs")),
  EM = high(rating_probs))

# Results for the plots in codes/plots/
write.csv(data.frame(nodule_id = test$nodule_id, mal_median = y,
                     pred_primary = pred_primary, pred_secondary = pred_secondary,
                     pred_em = pred_em, prob_primary = test_prob[["Naive primary"]],
                     prob_secondary = test_prob[["Naive secondary"]],
                     prob_em = test_prob[["EM"]]),
          "data/processed/test_predictions.csv", row.names = FALSE)
write.csv(data.frame(model = names(test_prob),
                     tau = sapply(names(test_prob), function(m)
                       choose_tau(train_prob[[m]], truth_train))),
          "data/processed/test_tau.csv", row.names = FALSE)

cat("\nTest nodules with grade >= 4:", sum(truth_test), "of", length(truth_test), "\n")
cat("Model: tau; training sensitivity, workload; test sensitivity, workload, precision, missed\n")
for (m in names(test_prob)) {
  tau <- choose_tau(train_prob[[m]], truth_train)
  flag_train <- train_prob[[m]] >= tau
  flag <- test_prob[[m]] >= tau
  cat(sprintf("%s: %.4f; %.3f, %.3f; %.3f, %.3f, %.3f, %d\n", m, tau,
              mean(flag_train[truth_train]), mean(flag_train),
              mean(flag[truth_test]), mean(flag), mean(truth_test[flag]),
              sum(truth_test & !flag)))
}

# ======================================
#   RISK SCORE RANGE OF EACH PREDICTED GRADE
# ======================================

# Naive primary model: where each grade is the most likely one
eta <- seq(-8, 8, by = 0.001)
cum <- cbind(0, sapply(fit_primary$zeta, function(a) plogis(a - eta)), 1)
best_grade <- max.col(cum[, -1] - cum[, -(K + 1)], "first")
cat("\nRisk score range per predicted grade (naive primary):\n")
for (k in sort(unique(best_grade)))
  cat(sprintf("Grade %d: %.3f to %.3f\n", k, min(eta[best_grade == k]), max(eta[best_grade == k])))
cat(sprintf("Largest P(grade 2) at any risk score: %.3f\n",
            max((cum[, -1] - cum[, -(K + 1)])[, 2])))
