# Run from the repository root:
#   Rscript codes/plots/data-description/correlation.R

nod <- read.csv("data/processed/nodules_train.csv", stringsAsFactors = FALSE)

GEOM <- c("log_volume", "elongation", "flatness", "solidity", "compactness",
          "spiculation", "lobulation", "margin", "texture", "sphericity")
LABELS <- c(log_volume = "log Volume", elongation = "Elongation",
            flatness = "Flatness", solidity = "Solidity",
            compactness = "Compactness", spiculation = "Spiculation",
            lobulation = "Lobulation", margin = "Margin", texture = "Texture",
            sphericity = "Sphericity")

cm <- cor(nod[GEOM], use = "pairwise.complete.obs", method = "spearman")

pdf("codes/plots/data-description/correlation.pdf", width = 7.5, height = 7)
par(mar = c(9, 9, 1, 1))
p <- length(GEOM)
pal <- colorRampPalette(c("#b2182b", "#f7f7f7", "#2166ac"))(101)
image(1:p, 1:p, cm[, p:1], zlim = c(-1, 1), col = pal, axes = FALSE,
      xlab = "", ylab = "")
axis(1, 1:p, LABELS[GEOM], las = 2, tick = FALSE)
axis(2, 1:p, rev(LABELS[GEOM]), las = 1, tick = FALSE)
text(rep(1:p, p), rep(p:1, each = p), sprintf("%.2f", cm), cex = 0.7)
invisible(dev.off())
