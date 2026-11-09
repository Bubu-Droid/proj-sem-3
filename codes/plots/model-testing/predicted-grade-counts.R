# Run from the repository root:
#   Rscript codes/plots/model-testing/predicted-grade-counts.R

d <- read.csv("data/processed/test_predictions.csv")
tab <- sapply(list(Actual = d$mal_median, Primary = d$pred_primary,
                   Secondary = d$pred_secondary, EM = d$pred_em), tabulate, nbins = 5)
cols <- c("grey60", "#4c72b0", "#dd8452", "#1a8060")

pdf("codes/plots/model-testing/predicted-grade-counts.pdf", width = 7.5, height = 4)
par(mar = c(4.2, 4.2, 0.8, 0.8), las = 1, col.axis = "grey25", col.lab = "grey25", fg = "grey60")
bp <- barplot(t(tab), beside = TRUE, col = cols, border = "white", names.arg = 1:5,
              xlab = "Grade", ylab = "Test nodules", ylim = c(0, 760))
abline(h = seq(200, 600, 200), col = "grey90")
barplot(t(tab), beside = TRUE, col = cols, border = "white", names.arg = 1:5, add = TRUE,
        axes = FALSE)
text(bp, t(tab) + 22, t(tab), cex = 0.62, col = "grey20")
legend("topright", c("Actual median", "Naive primary", "Naive secondary", "EM"), fill = cols,
       border = NA, bty = "n", cex = 0.85, text.col = "grey20")
invisible(dev.off())
