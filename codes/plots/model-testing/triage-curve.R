# Run from the repository root:
#   Rscript codes/plots/model-testing/triage-curve.R

d <- read.csv("data/processed/test_predictions.csv")
tau <- read.csv("data/processed/test_tau.csv")
truth <- d$mal_median >= 4
models <- list(list(p = d$prob_primary, name = "Naive primary", col = "#4c72b0", lty = 1, pch = 16),
               list(p = d$prob_secondary, name = "Naive secondary", col = "#dd8452", lty = 2, pch = 15),
               list(p = d$prob_em, name = "EM", col = "#1a8060", lty = 4, pch = 17))
curve_at <- function(p, cuts) t(sapply(cuts, function(k) c(mean(p >= k), mean(p[truth] >= k))))

pdf("codes/plots/model-testing/triage-curve.pdf", width = 6, height = 4.6)
par(mar = c(4.2, 4.2, 0.8, 0.8), las = 1, col.axis = "gray25", col.lab = "gray25", fg = "gray60")
plot(c(0, 1), c(0, 1), type = "n", xlab = "Workload (share of nodules flagged)",
     ylab = "Sensitivity (share of grade 4+ flagged)", bty = "n")
abline(h = seq(0, 1, 0.2), v = seq(0, 1, 0.2), col = "gray92")
abline(0, 1, col = "gray60", lty = 3)
text(0.72, 0.65, "Random flagging", col = "gray35", cex = 0.8, srt = 37)
for (m in models) {
  xy <- curve_at(m$p, c(Inf, sort(unique(m$p), decreasing = TRUE)))
  lines(xy[, 1], xy[, 2], col = m$col, lwd = 2, lty = m$lty)
}
for (m in models) {
  xy <- curve_at(m$p, tau$tau[tau$model == m$name])
  points(xy[1], xy[2], pch = m$pch, col = m$col, cex = 1.3)
}
legend("bottomright", sapply(models, `[[`, "name"), col = sapply(models, `[[`, "col"),
       lty = sapply(models, `[[`, "lty"), pch = sapply(models, `[[`, "pch"), lwd = 2, bty = "n",
       cex = 0.85, text.col = "gray20", title = "Points mark the triage rule")
invisible(dev.off())
