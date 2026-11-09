# Run from the repository root:
#   Rscript codes/plots/data-description/features-by-grade.R

nod <- read.csv("data/processed/nodules_train_trimmed.csv", stringsAsFactors = FALSE)

GEOM <- c("log_volume", "elongation", "flatness", "solidity", "compactness")
LABELS <- c(log_volume = "log Volume", elongation = "Elongation",
            flatness = "Flatness", solidity = "Solidity",
            compactness = "Compactness")
BLUE <- "#4c72b0"

pdf("codes/plots/data-description/features-by-grade.pdf", width = 10, height = 5.5)
par(mfrow = c(2, 3), mar = c(4, 4, 2.5, 1), las = 1, col.axis = "grey25", col.lab = "grey25", font.main = 1, cex.main = 1, col.main = "grey20")
for (v in GEOM)
  boxplot(nod[[v]] ~ nod$mal_median, col = BLUE, border = "grey25",
          outcex = 0.4, main = LABELS[[v]], xlab = "Median malignancy", ylab = "")
invisible(dev.off())
