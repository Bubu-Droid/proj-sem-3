# Run from the repository root:
#   Rscript codes/plots/data-description/mal-distribution.R

nod <- read.csv("data/processed/nodules_train_trimmed.csv", stringsAsFactors = FALSE)

BLUE <- "#4c72b0"; ORANGE <- "#dd8452"
GRADES <- c("1", "2", "3", "4", "5")

med_tab <- table(factor(nod$mal_median, levels = 1:5))
all_tab <- table(factor(unlist(nod[paste0("R", 1:4)]), levels = 1:5))

pdf("codes/plots/data-description/mal-distribution.pdf", width = 9, height = 3.8)
par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1), las = 1, col.axis = "gray25", col.lab = "gray25", font.main = 1, cex.main = 1, col.main = "gray20")
barplot(med_tab, col = BLUE, border = NA, names.arg = GRADES,
        main = "Median malignancy per nodule", xlab = "Grade", ylab = "Nodules")
barplot(all_tab, col = ORANGE, border = NA, names.arg = GRADES,
        main = "Individual ratings", xlab = "Grade", ylab = "Ratings")
invisible(dev.off())
