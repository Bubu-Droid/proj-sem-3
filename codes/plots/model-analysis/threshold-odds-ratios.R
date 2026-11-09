# Run from the repository root:
#   Rscript codes/plots/model-analysis/threshold-odds-ratios.R

library(ordinal)

data <- read.csv("data/processed/nodules_train_trimmed.csv", stringsAsFactors = FALSE)
PRIMARY <- c("log_volume", "elongation", "flatness", "solidity", "compactness")
SECONDARY <- c("spiculation", "lobulation", "margin", "texture", "sphericity")
data <- data[complete.cases(data[c(PRIMARY, SECONDARY, "mal_median")]), ]
data[PRIMARY] <- scale(data[PRIMARY])
data$y <- factor(data$mal_median, levels = 1:5, ordered = TRUE)

free <- clm(y ~ 1, nominal = ~ log_volume + elongation + flatness + solidity + compactness,
            data = data, link = "logit")
po <- clm(y ~ log_volume + elongation + flatness + solidity + compactness, data = data,
          link = "logit")
b <- coef(free)

pdf("codes/plots/model-analysis/threshold-odds-ratios.pdf", width = 9, height = 3)
par(mfrow = c(1, 5), mar = c(4.2, 2.2, 2, 0.4), oma = c(0, 2.4, 0, 0), las = 1,
    col.axis = "grey25", col.lab = "grey25", fg = "grey60")
for (v in PRIMARY) {
  or <- exp(-b[paste0(1:4, "|", 2:5, ".", v)])
  plot(1:4, or, type = "o", pch = 16, lwd = 2, col = "#4c72b0", log = "y", ylim = c(0.5, 7),
       xaxt = "n", xlab = "Threshold", ylab = "", bty = "n", cex.axis = 0.85)
  axis(1, 1:4, c(">1", ">2", ">3", ">4"), cex.axis = 0.85)
  abline(h = 1, col = "grey55", lty = 2)
  segments(1, exp(coef(po)[v]), 4, exp(coef(po)[v]), col = "grey35", lty = 3, lwd = 1.5)
  title(v, font.main = 1, cex.main = 1, family = "mono", col.main = "grey20")
}
mtext("Odds ratio per 1 SD (log scale)", side = 2, outer = TRUE, line = 0.8, las = 0,
      col = "grey25", cex = 0.75)
invisible(dev.off())
