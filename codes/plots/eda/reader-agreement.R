# Run from the repository root:
#   Rscript codes/plots/eda/reader-agreement.R

nod <- read.csv("data/processed/nodules_train.csv", stringsAsFactors = FALSE)

BLUE <- "#4c72b0"

n <- rowSums(nod[paste0("C", 1:5)])
multi <- n >= 2
rng <- apply(nod[paste0("R", 1:4)], 1, function(r) diff(range(r, na.rm = TRUE)))
rng_tab <- table(n_readers = n[multi], range = factor(rng[multi], levels = 0:4))

pdf("codes/plots/eda/reader-agreement.pdf", width = 7, height = 4)
par(mar = c(4, 4, 1, 1))
barplot(t(rng_tab), beside = TRUE, col = colorRampPalette(c("#cfd8ea", BLUE))(5),
        border = NA, xlab = "Number of readers", ylab = "Nodules",
        legend.text = paste("range", 0:4),
        args.legend = list(x = "topleft", bty = "n", cex = 0.8))
invisible(dev.off())
