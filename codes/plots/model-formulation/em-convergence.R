# Run from the repository root:
#   Rscript codes/plots/model-formulation/em-convergence.R

tr <- read.csv("data/processed/em_trace.csv")

pdf("codes/plots/model-formulation/em-convergence.pdf", width = 7, height = 3.6)
par(mar = c(4.2, 4.8, 0.8, 0.8), las = 1, col.axis = "gray25", col.lab = "gray25", fg = "gray60")
plot(tr$iteration, tr$loglik, type = "l", col = "#4c72b0", lwd = 2, log = "x", bty = "n",
     xlab = "Iteration (log scale)", ylab = "")
mtext("Observed-data log-likelihood", side = 2, line = 3.6, las = 0, col = "gray25")
abline(h = max(tr$loglik), col = "gray75", lty = 2)
invisible(dev.off())
