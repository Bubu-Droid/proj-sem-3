# Overall summary statistics of the model variables (training nodules).
#
# Run from the repository root:
#   Rscript codes/r-codes/summary-all.R

data <- read.csv("data/processed/nodules_train.csv", stringsAsFactors = FALSE)

PRIMARY <- c("log_volume", "elongation", "flatness", "solidity", "compactness")
SECONDARY <- c("spiculation", "lobulation", "margin", "texture", "sphericity")

# ======================================
#   AGGREGATE DATA
# ======================================

for (v in c(PRIMARY, SECONDARY)) {
  x <- data[[v]]
  cat(sprintf("%s: %.3f & %.3f & %.3f & %.3f (missing: %d)\n", v,
              mean(x, na.rm = TRUE), min(x, na.rm = TRUE),
              max(x, na.rm = TRUE), sd(x, na.rm = TRUE), sum(is.na(x))))
}

cat("\nNodules with flatness exactly 0:", sum(data$flatness == 0), "\n")
cat("Nodules with every variable observed:",
    sum(complete.cases(data[c(PRIMARY, SECONDARY)])), "of", nrow(data), "\n")
