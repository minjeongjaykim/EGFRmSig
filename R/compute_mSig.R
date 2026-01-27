#' Compute the predicted EGFR mSig score using predefined LUAD centroids
#'
#' @param centroids All_subtypes_centroids (1020 genes x 10 means), z-score space
#' @param expr gene expression matrix (gene symbols x sample IDs), values = RSEM TPM
#' @return pred.ind (1=WT-like, 2=mt-like), pred.score (distance)
#' @export

compute_mSig <- function(centroids=NULL, expr, meanCol="all", fname="Pred_mSig", export=TRUE) {

  if (is.null(centroids)) {
    centroids <- All_subtypes_centroids
  }

  ## ----- Preprocessing
  # Log transform
  expr_log2 <- log2(expr+1)

  # Gene-wise z-score (handle sd=0 to avoid Inf/NaN)
  sdv <- apply(expr_log2, 1, sd, na.rm=TRUE)
  sdv_floor <- quantile(sdv[sdv>0], 0.01, na.rm=TRUE)

  expr_z <- t(scale(t(expr_log2)))
  expr_z[is.na(expr_z)] <- 0      # sd=0
  expr_z[sdv<=sdv_floor, ] <- 0   # near-zero sd

  # Choose which centroid columns to use (keep same style: column names)
  # If you want the exact old behavior (using meanCol="subtype.i")
  if (!(meanCol %in% c("subtype.i", "all"))) {
    stop(meanCol," is not an option for meanCol.")
  }

  wt_col <- paste0("mean_", meanCol, "_wt")
  mt_col <- paste0("mean_", meanCol, "_mt")
  stopifnot(
    wt_col %in% colnames(All_subtypes_centroids),
    mt_col %in% colnames(All_subtypes_centroids)
  )

  # Match genes between centroids and expression (keep your "sub.centr.idx" style)
  match.idx <- unlist(sapply(rownames(All_subtypes_centroids), function(x) which(x==rownames(expr_z))))
  match.name <- names(unlist(sapply(rownames(All_subtypes_centroids), function(x) which(x==rownames(expr_z)))))
  sub.centr.idx <- as.numeric(sapply(match.name, function(x) which(x==rownames(All_subtypes_centroids))))

  centr.by.mut <- list()
  centr.by.mut[[1]] <- All_subtypes_centroids[sub.centr.idx, wt_col]
  centr.by.mut[[2]] <- All_subtypes_centroids[sub.centr.idx, mt_col]
  stopifnot(
    identical( names(centr.by.mut[[1]]), names(centr.by.mut[[2]]) )
  )

  message("Overlapping genes (n_genes_used): ", length(match.idx))

  # Build the matched expression matrix in the SAME gene order used above
  z.expr.matched <- expr_z[match.idx, , drop=FALSE]
  stopifnot(
    identical( rownames(z.expr.matched), rownames(All_subtypes_centroids)[sub.centr.idx] )
  )

  ## ----- Distance-based prediction (keep your for-loop style)
  pred.ind <- rep(NA_integer_, ncol(z.expr.matched))
  pred.score <- rep(NA_real_, ncol(z.expr.matched))

  for (i in seq_len(ncol(z.expr.matched))) {
    d1 <- sum((centr.by.mut[[1]]-z.expr.matched[, i])^2, na.rm=TRUE)
    d2 <- sum((centr.by.mut[[2]]-z.expr.matched[, i])^2, na.rm=TRUE)

    if (is.na(d1) | is.na(d2)) next

    if (d2 < d1) {
      pred.ind[i] <- 2
      pred.score[i] <- d2
    } else {
      pred.ind[i] <- 1
      pred.score[i] <- d1
    }
  }

  # pred.ind.bak <- pred.ind

  # Build output table
  output <- data.frame(
    SampleID = colnames(z.expr.matched),
    # EGFR_mSig_class = ifelse(pred.ind==2, "EGFR_mSig_mut_like", "EGFR_mSig_WT_like"),
    EGFR_mSig_class = dplyr::case_when(
      pred.ind==1 ~ "EGFR_mSig_WT_like",
      pred.ind==2 ~ "EGFR_mSig_mut_like",
      TRUE        ~ NA_character_
    ),
    EGFR_mSig_distance = pred.score,
    # n_genes_used = length(match.idx),
    stringsAsFactors = FALSE
  )

  # Export output to txt file
  if (export == TRUE) {
    write.table(
      output,
      file = paste0(fname, "_", Sys.Date(), ".txt"), # fname_yyyy-mm-dd.txt
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    message("Output saved to: ", normalizePath(paste0(fname, "_", Sys.Date(), ".txt")))
  }

  else if (export == FALSE) {
    return(output)
  }

  else {
    stop(export," is not an option for export.")
  }
}
