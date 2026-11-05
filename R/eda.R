# Exploratory analysis of the training nodules (STEPS.md 3.6).
#
# Outputs: figures/*.pdf and the numbers quoted in main.tex (printed).
#
# Run from the repository root:
#   Rscript R/eda.R

nod <- read.csv("data/processed/nodules_train.csv", stringsAsFactors = FALSE)
raw <- read.csv("data/raw/annotations.csv", stringsAsFactors = FALSE)
dir.create("figures", showWarnings = FALSE)

GEOM <- c("diameter", "log_volume", "log_surface_area", "log_convex_volume",
          "bbox_z", "n_slices", "elongation", "flatness", "solidity",
          "compactness")
SEM <- c("subtlety", "sphericity", "margin", "lobulation", "spiculation",
         "texture")

BLUE <- "#4c72b0"; ORANGE <- "#dd8452"
GRADES <- c("1", "2", "3", "4", "5")

# ---- Summary statistics of the geometric features -------------------------
cat("== geometric features (training nodules)\n")
summ <- t(sapply(nod[GEOM], function(x)
  c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE),
    min = min(x, na.rm = TRUE), max = max(x, na.rm = TRUE),
    n_na = sum(is.na(x)))))
print(round(summ, 3))

# ---- (a) Response distribution --------------------------------------------
cat("\n== malignancy\n")
med_tab <- table(factor(nod$mal_median, levels = 1:5))
all_tab <- table(factor(unlist(nod[paste0("R", 1:4)]), levels = 1:5))
print(rbind(median = med_tab, all_ratings = all_tab))
print(round(rbind(median = prop.table(med_tab), all_ratings = prop.table(all_tab)), 3))

pdf("figures/mal-distribution.pdf", width = 9, height = 3.8)
par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1))
barplot(med_tab, col = BLUE, border = NA, names.arg = GRADES,
        main = "Median malignancy per nodule", xlab = "Grade", ylab = "Nodules")
barplot(all_tab, col = ORANGE, border = NA, names.arg = GRADES,
        main = "All individual ratings", xlab = "Grade", ylab = "Ratings")
invisible(dev.off())

# ---- (b) Features by grade -------------------------------------------------
cat("\n== median of each feature by mal_median\n")
print(round(sapply(nod[GEOM], function(x) tapply(x, nod$mal_median, median, na.rm = TRUE)), 3))

LABELS <- c(diameter = "Diameter (mm)", log_volume = "log Volume",
            log_surface_area = "log Surface area",
            log_convex_volume = "log Convex volume", bbox_z = "Height (mm)",
            n_slices = "Slices", elongation = "Elongation",
            flatness = "Flatness", solidity = "Solidity",
            compactness = "Compactness")
pdf("figures/features-by-grade.pdf", width = 10, height = 8)
par(mfrow = c(3, 3), mar = c(4, 4, 2.5, 1))
for (v in setdiff(GEOM, "n_slices"))
  boxplot(nod[[v]] ~ nod$mal_median, col = BLUE, border = "grey25",
          outcex = 0.4, main = LABELS[[v]], xlab = "Median malignancy", ylab = "")
invisible(dev.off())

# ---- (c) Correlation between predictors ------------------------------------
cm <- cor(nod[GEOM], use = "pairwise.complete.obs", method = "spearman")
cat("\n== Spearman correlation of geometric features\n")
print(round(cm, 2))

pdf("figures/correlation.pdf", width = 7.5, height = 7)
par(mar = c(9, 9, 1, 1))
p <- length(GEOM)
pal <- colorRampPalette(c("#b2182b", "#f7f7f7", "#2166ac"))(101)
image(1:p, 1:p, cm[, p:1], zlim = c(-1, 1), col = pal, axes = FALSE,
      xlab = "", ylab = "")
axis(1, 1:p, LABELS[GEOM], las = 2, tick = FALSE)
axis(2, 1:p, rev(LABELS[GEOM]), las = 1, tick = FALSE)
text(rep(1:p, p), rep(p:1, each = p), sprintf("%.2f", cm), cex = 0.7)
invisible(dev.off())

# n_slices depends on the scan protocol as well as the nodule
cat(sprintf("\nSpearman(n_slices, slice_spacing) over annotations: %.2f\n",
            cor(raw$n_slices, raw$slice_spacing, method = "spearman")))
cat(sprintf("Spearman(bbox_z, slice_spacing) over annotations:   %.2f\n",
            cor(raw$bbox_z, raw$slice_spacing, method = "spearman")))

# ---- (d) Reader agreement --------------------------------------------------
C <- as.matrix(nod[paste0("C", 1:5)])
n <- rowSums(C)
multi <- n >= 2
rng <- apply(nod[paste0("R", 1:4)], 1, function(r) diff(range(r, na.rm = TRUE)))

cat("\n== reader agreement (nodules with >= 2 readers)\n")
cat(sprintf("multi-reader nodules: %d of %d\n", sum(multi), nrow(nod)))
rng_tab <- table(n_readers = n[multi], range = factor(rng[multi], levels = 0:4))
print(addmargins(rng_tab))

# Fleiss' kappa for a varying number of raters per subject
P_i <- rowSums(C * (C - 1)) / (n * (n - 1))
Pbar <- mean(P_i[multi])
p_r <- colSums(C[multi, ]) / sum(C[multi, ])
Pe <- sum(p_r^2)
kappa <- (Pbar - Pe) / (1 - Pe)
cat(sprintf("Fleiss kappa: Pbar = %.3f, Pe = %.3f, kappa = %.3f\n", Pbar, Pe, kappa))

pdf("figures/reader-agreement.pdf", width = 7, height = 4)
par(mar = c(4, 4, 1, 1))
barplot(t(rng_tab), beside = TRUE, col = colorRampPalette(c("#cfd8ea", BLUE))(5),
        border = NA, xlab = "Number of readers", ylab = "Nodules",
        legend.text = paste("range", 0:4),
        args.legend = list(x = "topleft", bty = "n", cex = 0.8))
invisible(dev.off())

# ---- Semantic features vs malignancy (circularity check) ------------------
cat("\n== Spearman correlation with mal_median_raw\n")
print(round(cor(nod[c("diameter", "log_volume", SEM)], nod$mal_median_raw,
                method = "spearman"), 2))
