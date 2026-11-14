# Build the nodule-level analysis tables from data/raw/annotations.csv.
#
# Outputs (data/processed/):
#   nodules.csv        master table, one row per nodule, all columns + split
#   nodules_train.csv  modeling table, training patients
#   nodules_test.csv   modeling table, test patients
#   nodules_train_trimmed.csv  training table without outlying ratings
#   median_ties.csv    nodules whose median malignancy fell between two grades
#
# Run from the repository root:
#   Rscript codes/r-codes/clean-data.R

raw <- read.csv("data/raw/annotations.csv", stringsAsFactors = FALSE)

GEOM <- c("diameter", "surface_area", "volume", "convex_volume", "n_voxels",
          "n_slices", "bbox_x", "bbox_y", "bbox_z", "elongation", "flatness",
          "solidity")
SEMANTIC <- c("subtlety", "sphericity", "margin", "lobulation", "spiculation",
              "texture")
NOMINAL <- c("calcification", "internalStructure")  # kept in master only
K <- 5

# ---- Clusters with more than 4 annotations ---------------------------------
# LIDC had 4 readers per scan, so a cluster of 5+ means cluster_annotations()
# merged neighboring lesions. Drop them.
too_big <- unique(raw$nodule_id[raw$n_readers > 4])
cat(sprintf("dropping %d clusters with > 4 annotations (%d annotations)\n",
            length(too_big), sum(raw$nodule_id %in% too_big)))
raw <- raw[!raw$nodule_id %in% too_big, ]

# ---- Aggregate to one row per nodule ---------------------------------------
mode_int <- function(x) {
  tab <- table(x)
  as.integer(names(tab)[which.max(tab)])            # ties -> smallest value
}

by_nod <- split(raw, raw$nodule_id)

nod <- do.call(rbind, lapply(by_nod, function(d) {
  r <- d$malignancy
  med <- median(r)
  data.frame(
    nodule_id  = d$nodule_id[1],
    patient_id = d$patient_id[1],
    scan_id    = d$scan_id[1],
    n_readers  = nrow(d),
    as.list(colMeans(d[GEOM], na.rm = TRUE)),
    as.list(sapply(d[SEMANTIC], median)),
    calcification     = mode_int(d$calcification),
    internalStructure = mode_int(d$internalStructure),
    R1 = r[1], R2 = r[2], R3 = r[3], R4 = r[4],      # NA past n_readers
    as.list(setNames(tabulate(r, nbins = K), paste0("C", 1:K))),
    mal_median_raw = med,
    mal_median_tie = med != round(med),
    mal_median     = ceiling(med)                   # ties rounded up
  )
}))
rownames(nod) <- NULL

# colMeans of an all-NA column gives NaN; make it NA
nod[GEOM] <- lapply(nod[GEOM], function(x) replace(x, is.nan(x), NA))

# ---- Derived features ------------------------------------------------------
nod$log_volume        <- log(nod$volume)
nod$log_surface_area  <- log(nod$surface_area)
nod$log_convex_volume <- log(nod$convex_volume)
nod$compactness       <- nod$volume^(2/3) / nod$surface_area

# ---- Train/test split by patient (70/30) -----------------------------------
set.seed(20261105)
patients <- sort(unique(nod$patient_id))
train_pat <- sample(patients, round(0.7 * length(patients)))
nod$split <- ifelse(nod$patient_id %in% train_pat, "train", "test")

# ---- Write -----------------------------------------------------------------
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
write.csv(nod, "data/processed/nodules.csv", row.names = FALSE)

model_cols <- setdiff(names(nod), c(NOMINAL, "split"))
write.csv(nod[nod$split == "train", model_cols],
          "data/processed/nodules_train.csv", row.names = FALSE)
write.csv(nod[nod$split == "test", model_cols],
          "data/processed/nodules_test.csv", row.names = FALSE)

# ---- Trimmed training table for the EM model --------------------------------
# Removes ratings at least 2 grades from the nodule's unrounded median, then
# recomputes the median columns from the ratings that remain. A nodule with no
# ratings left has NA in them.
trim <- nod[nod$split == "train", model_cols]
R_cols <- paste0("R", 1:4)
C_cols <- paste0("C", 1:K)
total <- rowSums(trim[C_cols])
R <- as.matrix(trim[R_cols])
R[abs(sweep(R, 1, trim$mal_median_raw)) >= 2] <- NA
trim[R_cols] <- R
trim[C_cols] <- t(apply(R, 1, function(r) tabulate(r[!is.na(r)], K)))
trim$n_ratings_removed <- total - rowSums(trim[C_cols])
trim$mal_median_raw <- apply(R, 1, function(r)
  if (all(is.na(r))) NA_real_ else median(r, na.rm = TRUE))
trim$mal_median_tie <- trim$mal_median_raw != round(trim$mal_median_raw)
trim$mal_median <- ceiling(trim$mal_median_raw)
write.csv(trim, "data/processed/nodules_train_trimmed.csv", row.names = FALSE)
cat(sprintf("trimmed: %d of %d training ratings removed (%.1f%%), %d nodules affected\n",
            sum(trim$n_ratings_removed), sum(total),
            100 * sum(trim$n_ratings_removed) / sum(total),
            sum(trim$n_ratings_removed > 0)))

ties <- nod[nod$mal_median_tie,
            c("nodule_id", "patient_id", "n_readers", paste0("R", 1:4),
              "mal_median_raw", "mal_median")]
write.csv(ties, "data/processed/median_ties.csv", row.names = FALSE)

# ---- Summary ---------------------------------------------------------------
cat(sprintf("nodules %d, patients %d (train %d / test %d nodules)\n",
            nrow(nod), length(patients),
            sum(nod$split == "train"), sum(nod$split == "test")))
cat("readers per nodule:\n"); print(table(nod$n_readers))
cat("mal_median:\n"); print(table(nod$mal_median))
cat(sprintf("median ties: %d nodules\n", nrow(ties)))
print(table(n_readers = ties$n_readers, raw_median = ties$mal_median_raw))
cat("missing values per column:\n")
na <- colSums(is.na(nod[c(GEOM, SEMANTIC)])); print(na[na > 0])
