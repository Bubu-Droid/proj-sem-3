# Mean (SD) of the model variables by median malignancy grade (training nodules).
#
# Run from the repository root:
#   Rscript codes/r-codes/summary-grades.R

data <- read.csv("data/processed/nodules_train.csv", stringsAsFactors = FALSE)

PRIMARY <- c("log_volume", "elongation", "flatness", "solidity", "compactness")
SECONDARY <- c("spiculation", "lobulation", "margin", "texture", "sphericity")
GRADES <- 1:5

# ======================================
#   GRADE COUNTS
# ======================================

counts <- table(factor(data$mal_median, levels = GRADES))
cat("Nodules:", paste(counts, collapse = " & "), "\n")
cat("Share (%):", paste(sprintf("%.1f", 100 * prop.table(counts)),
                        collapse = " & "), "\n\n")

# ======================================
#   GRADE DATA
# ======================================

for (v in c(PRIMARY, SECONDARY)) {
  cells <- sapply(GRADES, function(g) {
    x <- data[[v]][data$mal_median == g]
    sprintf("%.3f (%.3f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
  })
  cat(sprintf("%s: %s\n", v, paste(cells, collapse = " & ")))
}
