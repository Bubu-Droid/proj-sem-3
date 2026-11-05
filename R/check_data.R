# Sanity checks for the extracted and processed data (STEPS.md 2.3 and 3).
# Exits non-zero if any check fails; WARN lines are known, documented quirks.
#
# Run from the repository root:
#   Rscript R/check_data.R

raw <- read.csv("data/raw/annotations.csv", stringsAsFactors = FALSE)
nod <- read.csv("data/processed/nodules.csv", stringsAsFactors = FALSE)
trn <- read.csv("data/processed/nodules_train.csv", stringsAsFactors = FALSE)
tst <- read.csv("data/processed/nodules_test.csv", stringsAsFactors = FALSE)

n_fail <- 0
check <- function(ok, msg) {
  cat(if (isTRUE(ok)) "PASS " else "FAIL ", msg, "\n", sep = "")
  if (!isTRUE(ok)) n_fail <<- n_fail + 1
}
warn <- function(msg) cat("WARN ", msg, "\n", sep = "")

in_range <- function(x, lo, hi) all(!is.na(x) & x >= lo & x <= hi)

# ---- Raw annotations (Step 2.3) --------------------------------------------
cat("== data/raw/annotations.csv\n")
check(nrow(raw) == 6859, sprintf("6859 annotations (got %d)", nrow(raw)))
check(!anyDuplicated(raw$ann_id), "ann_id unique")
cat(sprintf("     patients %d, scans %d, nodule clusters %d\n",
            length(unique(raw$patient_id)), length(unique(raw$scan_id)),
            length(unique(raw$nodule_id))))

sizes <- table(raw$nodule_id)
check(all(raw$n_readers == sizes[raw$nodule_id]),
      "n_readers equals the cluster size")
check(!anyNA(raw$malignancy), "no missing malignancy")

ranges <- list(subtlety = 5, sphericity = 5, margin = 5, lobulation = 5,
               spiculation = 5, texture = 5, malignancy = 5, calcification = 6)
for (v in names(ranges))
  check(in_range(raw[[v]], 1, ranges[[v]]),
        sprintf("%s in 1..%d", v, ranges[[v]]))
bad_is <- sum(!raw$internalStructure %in% 1:4)
if (bad_is > 0) {
  warn(sprintf("internalStructure outside 1..4 in %d annotation(s); column is not modelled",
               bad_is))
} else check(TRUE, "internalStructure in 1..4")

check(all(raw[c("diameter", "surface_area", "volume", "n_voxels")] > 0),
      "diameter, surface_area, volume, n_voxels positive")
check(in_range(na.omit(raw$solidity), 1e-12, 1), "solidity in (0, 1]")
check(in_range(na.omit(raw$elongation), 0, 1) && in_range(na.omit(raw$flatness), 0, 1),
      "elongation, flatness in [0, 1]")
check(!any(sapply(raw, function(x) is.numeric(x) && any(is.infinite(x)))),
      "no infinite values")
warn(sprintf("%d annotation(s) with no convex hull (flat mask): convex_volume/solidity NA",
             sum(is.na(raw$solidity))))

scans_per_patient <- tapply(raw$scan_id, raw$patient_id, function(s) length(unique(s)))
cat(sprintf("     patients with more than one scan: %d\n", sum(scans_per_patient > 1)))

# ---- Processed nodule table (Step 3) ---------------------------------------
cat("\n== data/processed/nodules.csv\n")
big <- names(sizes)[sizes > 4]
check(nrow(nod) == length(sizes) - length(big),
      sprintf("one row per nodule after dropping %d clusters with > 4 annotations",
              length(big)))
check(!anyDuplicated(nod$nodule_id), "nodule_id unique")
check(!any(nod$nodule_id %in% big), "no over-merged clusters remain")
check(in_range(nod$n_readers, 1, 4), "n_readers in 1..4")

Rm <- as.matrix(nod[paste0("R", 1:4)])
Cm <- as.matrix(nod[paste0("C", 1:5)])
check(all(rowSums(!is.na(Rm)) == nod$n_readers), "non-missing R1..R4 equals n_readers")
check(all(rowSums(Cm) == nod$n_readers), "C1..C5 sum to n_readers")
check(all(t(apply(Rm, 1, function(r) tabulate(r, nbins = 5))) == Cm),
      "C1..C5 are the counts of R1..R4")

med <- apply(Rm, 1, median, na.rm = TRUE)
check(all(nod$mal_median_raw == med), "mal_median_raw is the median of R1..R4")
check(all(nod$mal_median == ceiling(med)), "mal_median rounds half-grades up")
check(all(nod$mal_median_tie == (med != round(med))), "mal_median_tie flags half-grades")
check(all(!nod$mal_median_tie | nod$n_readers %in% c(2, 4)),
      "ties only occur with 2 or 4 readers")

# Aggregation against the raw file
agg_d <- tapply(raw$diameter, raw$nodule_id, mean)
check(isTRUE(all.equal(as.vector(agg_d[nod$nodule_id]), nod$diameter)),
      "diameter is the reader mean")
agg_s <- tapply(raw$spiculation, raw$nodule_id, median)
check(isTRUE(all.equal(as.vector(agg_s[nod$nodule_id]), nod$spiculation)),
      "spiculation is the reader median")
check(isTRUE(all.equal(nod$log_volume, log(nod$volume))), "log_volume = log(volume)")

# ---- Train/test split ------------------------------------------------------
cat("\n== train/test split\n")
check(!length(intersect(trn$patient_id, tst$patient_id)), "no patient in both sets")
check(nrow(trn) + nrow(tst) == nrow(nod), "train + test = all nodules")
check(setequal(trn$nodule_id, nod$nodule_id[nod$split == "train"]),
      "train file matches split column")
check(!any(c("calcification", "internalStructure", "split") %in% names(trn)),
      "nominal ratings and split excluded from modelling tables")
check(all(c("calcification", "internalStructure") %in% names(nod)),
      "nominal ratings kept in master table")
frac <- length(unique(trn$patient_id)) / length(unique(nod$patient_id))
check(abs(frac - 0.7) < 0.01, sprintf("70%% of patients in train (got %.3f)", frac))

cat(sprintf("\n%d check(s) failed\n", n_fail))
quit(status = as.integer(n_fail > 0))
