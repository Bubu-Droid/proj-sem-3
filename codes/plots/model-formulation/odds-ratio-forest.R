# Run from the repository root:
#   Rscript codes/plots/model-formulation/odds-ratio-forest.R

library(MASS)

data <- read.csv("data/processed/nodules_train_trimmed.csv", stringsAsFactors = FALSE)
em <- read.csv("data/processed/em_coef.csv", stringsAsFactors = FALSE)
PRIMARY <- c("log_volume", "elongation", "flatness", "solidity", "compactness")
SECONDARY <- c("spiculation", "lobulation", "margin", "texture", "sphericity")
ALL <- c(PRIMARY, SECONDARY)
data <- data[complete.cases(data[c(ALL, "mal_median")]), ]
data[ALL] <- scale(data[ALL])
data$y <- factor(data$mal_median, levels = 1:5, ordered = TRUE)

ors <- function(vars) {
  fit <- polr(reformulate(vars, "y"), data = data, Hess = TRUE, method = "logistic")
  ci <- suppressMessages(confint(fit))
  list(est = exp(coef(fit)), low = exp(ci[, 1]), high = exp(ci[, 2]))
}
prim <- ors(PRIMARY)
sec <- ors(ALL)
em_or <- setNames(exp(em$beta), em$predictor)

BLUE <- "#4c72b0"; ORANGE <- "#dd8452"; GREEN <- "#1a8060"
y <- rev(seq_along(ALL))
names(y) <- ALL

pdf("codes/plots/model-formulation/odds-ratio-forest.pdf", width = 7.5, height = 5)
par(mar = c(4.2, 7, 0.8, 0.8), las = 1, col.axis = "grey25", col.lab = "grey25", fg = "grey60")
plot(NA, xlim = c(0.5, 6), ylim = c(0.4, length(ALL) + 0.6), log = "x", yaxt = "n",
     xlab = "Odds ratio per 1 SD increase (log scale)", ylab = "", bty = "n")
axis(2, y, ALL, tick = FALSE, family = "mono", cex.axis = 0.8)
abline(v = 1, col = "grey55", lty = 2)
segments(prim$low, y[PRIMARY] + 0.2, prim$high, y[PRIMARY] + 0.2, col = BLUE, lwd = 2)
points(prim$est, y[PRIMARY] + 0.2, pch = 16, col = BLUE, cex = 1.1)
segments(sec$low, y[ALL], sec$high, y[ALL], col = ORANGE, lwd = 2)
points(sec$est, y[ALL], pch = 15, col = ORANGE, cex = 1.1)
points(em_or[PRIMARY], y[PRIMARY] - 0.2, pch = 17, col = GREEN, cex = 1.1)
legend("bottomright", c("Naive primary (95% CI)", "Naive secondary (95% CI)", "EM (no CI)"),
       col = c(BLUE, ORANGE, GREEN), pch = c(16, 15, 17), lty = c(1, 1, NA), lwd = 2,
       bty = "n", cex = 0.85, text.col = "grey20")
invisible(dev.off())
