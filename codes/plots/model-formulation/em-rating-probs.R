# Run from the repository root:
#   Rscript codes/plots/model-formulation/em-rating-probs.R

Pi <- as.matrix(read.csv("data/processed/em_pi.csv"))
K <- nrow(Pi)
pal <- colorRampPalette(c("#f5f7fb", "#4c72b0", "#1c3768"))(101)

pdf("codes/plots/model-formulation/em-rating-probs.pdf", width = 5.2, height = 4.4)
par(mar = c(4.2, 4.2, 0.8, 0.8), las = 1, col.axis = "grey25", col.lab = "grey25")
image(1:K, 1:K, t(Pi[K:1, ]), zlim = c(0, 1), col = pal, axes = FALSE,
      xlab = "Rating", ylab = "Latent grade")
axis(1, 1:K, tick = FALSE)
axis(2, 1:K, K:1, tick = FALSE)
grid_x <- rep(1:K, each = K); grid_y <- rep(K:1, times = K)
val <- as.vector(Pi)
text(grid_x, grid_y, sprintf("%.2f", val), cex = 0.85,
     col = ifelse(val > 0.45, "white", "grey20"))
invisible(dev.off())
