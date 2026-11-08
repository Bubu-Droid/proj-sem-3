# Exploratory analysis of the training nodules (STEPS.md 3.6): the numbers quoted
# in main.tex. The figures are made by codes/plots/eda/*.R.
#
# Run from the repository root:
#   Rscript codes/r-codes/eda-summary.R

nod <- read.csv("data/processed/nodules_train.csv", stringsAsFactors = FALSE)
raw <- read.csv("data/raw/annotations.csv", stringsAsFactors = FALSE)

GEOM <- c("diameter", "log_volume", "log_surface_area", "log_convex_volume",
          "bbox_z", "n_slices", "elongation", "flatness", "solidity",
          "compactness")
SEM <- c("subtlety", "sphericity", "margin", "lobulation", "spiculation",
         "texture")

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

# ---- (b) Features by grade -------------------------------------------------
cat("\n== median of each feature by mal_median\n")
print(round(sapply(nod[GEOM], function(x) tapply(x, nod$mal_median, median, na.rm = TRUE)), 3))

# ---- (c) Correlation between predictors ------------------------------------
cm <- cor(nod[GEOM], use = "pairwise.complete.obs", method = "spearman")
cat("\n== Spearman correlation of geometric features\n")
print(round(cm, 2))

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

# ---- Semantic features vs malignancy (circularity check) ------------------
cat("\n== Spearman correlation with mal_median_raw\n")
print(round(cor(nod[c("diameter", "log_volume", SEM)], nod$mal_median_raw,
                method = "spearman"), 2))
