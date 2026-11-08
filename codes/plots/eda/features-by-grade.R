# Run from the repository root:
#   Rscript codes/plots/eda/features-by-grade.R

nod <- read.csv("data/processed/nodules_train.csv", stringsAsFactors = FALSE)

GEOM <- c("diameter", "log_volume", "log_surface_area", "log_convex_volume",
          "bbox_z", "n_slices", "elongation", "flatness", "solidity",
          "compactness")
LABELS <- c(diameter = "Diameter (mm)", log_volume = "log Volume",
            log_surface_area = "log Surface area",
            log_convex_volume = "log Convex volume", bbox_z = "Height (mm)",
            n_slices = "Slices", elongation = "Elongation",
            flatness = "Flatness", solidity = "Solidity",
            compactness = "Compactness")
BLUE <- "#4c72b0"

pdf("codes/plots/eda/features-by-grade.pdf", width = 10, height = 8)
par(mfrow = c(3, 3), mar = c(4, 4, 2.5, 1))
for (v in setdiff(GEOM, "n_slices"))
  boxplot(nod[[v]] ~ nod$mal_median, col = BLUE, border = "grey25",
          outcex = 0.4, main = LABELS[[v]], xlab = "Median malignancy", ylab = "")
invisible(dev.off())
