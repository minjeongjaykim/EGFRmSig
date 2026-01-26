# EGFR_mSig
LUAD EGFR mSig Score Calculation

## - Centroid: All_subtypes_centroids (1020 genes x 10), z-score space
## - Output: pred.ind (1 WT-like, 2 MT-like), pred.score (distance)


# Set your working directory
setwd("xxx")

# Load "All_subtypes_centroids"
load("xxx")
dim(All_subtypes_centroids)
head(All_subtypes_centroids)

# Load expression data
load("xxx")
expr <- luad.rsem.exp
dim(expr) #row = gene, column = sample

# Check structure
dim(expr)
head(expr[, 1:5],20)
rownames(expr)
colnames(expr)

# log2 transform
expr_log2 <- log2(expr + 1)

# gene-wise z-score
expr_z <- t(scale(t(expr_log2)))
dim(expr_z) 
head(expr_z)

## Choose which centroid columns to use (keep same style: column 1/2 or named columns)
wt_col <- colnames(All_subtypes_centroids)[3]  # "mean_all_wt"
mt_col <- colnames(All_subtypes_centroids)[4]  # "mean_all_wt"


## Match genes between centroids and expression (keep your "sub.centr.idx" style)
match.idx <- unlist(sapply(rownames(All_subtypes_centroids), function(x) which(x == rownames(expr_z))))
match.name <- names(unlist(sapply(rownames(All_subtypes_centroids), function(x) which(x == rownames(expr_z)))))
sub.centr.idx <- as.numeric(sapply(match.name, function(x) which(x == rownames(All_subtypes_centroids))))

centr.by.mut <- list()
centr.by.mut[[1]] <- All_subtypes_centroids[sub.centr.idx, wt_col]
centr.by.mut[[2]] <- All_subtypes_centroids[sub.centr.idx, mt_col]

print(length(centr.by.mut[[1]]))  # should be number of overlapping genes (e.g., 786)
print(length(centr.by.mut[[2]]))

## Build the matched expression matrix in the SAME gene order used above
##    Note: match.idx is positions in expr_z rownames that correspond to All_subtypes_centroids order
expr.matched <- expr_z[match.idx, , drop = FALSE]

## Gene-wise z-score normalization - If your ge matrix is already z-score normalized, skip this part. 
##   z = (x - mean) / sd, per gene across samples
#z.expr.matched <- (expr.matched - apply(expr.matched, 1, mean, na.rm = TRUE)) / apply(expr.matched, 1, sd, na.rm = TRUE)

## Handle sd==0 to avoid Inf/NaN (optional but recommended)
sdv <- apply(expr.matched, 1, sd, na.rm = TRUE)
sdv[is.na(sdv) | sdv == 0] <- 0.01
z.expr.matched <- (expr.matched - apply(expr.matched, 1, mean, na.rm = TRUE)) / sdv

## Distance-based prediction (keep your for-loop style)
pred.ind <- rep(NA_integer_, ncol(z.expr.matched))
pred.score <- rep(NA_real_, ncol(z.expr.matched))

for (i in seq_len(ncol(z.expr.matched))) {
  d1 <- sum((centr.by.mut[[1]] - z.expr.matched[, i])^2, na.rm = TRUE)
  d2 <- sum((centr.by.mut[[2]] - z.expr.matched[, i])^2, na.rm = TRUE)
  
  if (is.na(d1) | is.na(d2)) next
  
  if (d2 < d1) {
    pred.ind[i] <- 2
    pred.score[i] <- d2
  } else {
    pred.ind[i] <- 1
    pred.score[i] <- d1
  }
}

pred.ind.bak <- pred.ind

## Build output table (same as your TCGA table, but using ge_t colnames)
output <- data.frame(
  SampleID = colnames(z.expr.matched),
  EGFR_mSig_class = ifelse(pred.ind == 2, "EGFR_mSig_mut_like", "EGFR_mSig_WT_like"),
  EGFR_mSig_distance = pred.score,
  n_genes_used = length(match.idx),
  stringsAsFactors = FALSE
)

table(output$EGFR_mSig_class, useNA = "ifany")
summary(output$EGFR_mSig_distance)

write.table(
  output,
  file = "title.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

