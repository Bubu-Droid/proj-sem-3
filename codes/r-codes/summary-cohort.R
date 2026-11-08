# Counts quoted in the Data Description section of main.tex.
#
# Run from the repository root:
#   Rscript codes/r-codes/summary-cohort.R

raw <- read.csv("data/raw/annotations.csv", stringsAsFactors = FALSE)
nod <- read.csv("data/processed/nodules.csv", stringsAsFactors = FALSE)

# ======================================
#   EXTRACTED ANNOTATIONS
# ======================================

cat("Annotations:", nrow(raw), "\n")
cat("Nodule clusters:", length(unique(raw$nodule_id)), "\n")
cat("Scans:", length(unique(raw$scan_id)), "\n")
cat("Patients:", length(unique(raw$patient_id)), "\n\n")

too_big <- unique(raw$nodule_id[raw$n_readers > 4])
cat("Clusters with more than 4 annotations:", length(too_big), "\n")
cat("Annotations in those clusters:", sum(raw$nodule_id %in% too_big), "\n\n")

# ======================================
#   ANALYSIS NODULES
# ======================================

cat("Nodules:", nrow(nod), "\n")
cat("Patients:", length(unique(nod$patient_id)), "\n")
cat("Annotations kept:", sum(nod$n_readers), "\n\n")

cat("Nodules by number of readers:\n")
print(table(nod$n_readers))
cat("\n")

cat("Median ties:", sum(nod$mal_median_tie),
    sprintf("(%.1f%%)", 100 * mean(nod$mal_median_tie)), "\n\n")

# ======================================
#   TRAIN/TEST SPLIT
# ======================================

for (s in c("train", "test")) {
  d <- nod[nod$split == s, ]
  cat(s, "- nodules:", nrow(d), "patients:", length(unique(d$patient_id)), "\n")
}
