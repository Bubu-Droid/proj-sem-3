# Run from the repository root:
#   Rscript codes/plots/model-formulation/probability-curves.R

library(MASS)

data <- read.csv("data/processed/nodules_train_trimmed.csv", stringsAsFactors = FALSE)
PRIMARY <- c("log_volume", "elongation", "flatness", "solidity", "compactness")
SECONDARY <- c("spiculation", "lobulation", "margin", "texture", "sphericity")
data <- data[complete.cases(data[c(PRIMARY, SECONDARY, "mal_median")]), ]
centre <- mean(data$log_volume); spread <- sd(data$log_volume)
data[PRIMARY] <- scale(data[PRIMARY])
data$y <- factor(data$mal_median, levels = 1:5, ordered = TRUE)
fit <- polr(reformulate(PRIMARY, "y"), data = data, method = "logistic")

# Other predictors at their means, which are 0 after standardising
x <- seq(3, 9.5, length.out = 300)
grid <- data.frame(log_volume = (x - centre) / spread, elongation = 0, flatness = 0,
                   solidity = 0, compactness = 0)
P <- predict(fit, grid, type = "probs")

# Ordered grades: one hue, light to dark
ramp <- c("#a9bddc", "#7f9bc8", "#4c72b0", "#2f5088", "#16294d")

pdf("codes/plots/model-formulation/probability-curves.pdf", width = 7.5, height = 4.2)
par(mar = c(4.2, 4.2, 0.8, 0.8), las = 1, col.axis = "grey25", col.lab = "grey25", fg = "grey60")
plot(range(x), c(0, 1), type = "n", xlab = "log volume (mm\u00b3)",
     ylab = "Predicted probability", bty = "n")
abline(h = seq(0, 1, 0.2), col = "grey90")
for (k in 1:5) lines(x, P[, k], col = ramp[k], lwd = 2.2)
for (k in 1:5) {
  i <- c(55, 70, which.max(P[, 3]), which.max(P[, 4]), which.max(P[, 5]))[k]
  text(x[i], P[i, k] + c(-0.04, 0.04, 0.035, 0.035, 0.035)[k], paste("Grade", k),
       col = "grey20", cex = 0.85)
}
invisible(dev.off())
