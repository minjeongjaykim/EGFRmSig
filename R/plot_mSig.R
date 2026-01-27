#' Plot the distribution of EGFR mSig distance by class
#'
#' @param res output object returned by \code{compute_mSig()}
#' @return violin and boxplot
#' @import ggplot2
#' @export

plot_mSig <- function(res) {

  # Fix order: WT, mut
  res$EGFR_mSig_class <- factor(
    res$EGFR_mSig_class,
    levels = c("EGFR_mSig_WT_like", "EGFR_mSig_mut_like")
  )

  # n by class
  label_df <- res %>%
    count(EGFR_mSig_class) %>%
    mutate(label = paste0(EGFR_mSig_class, "\n(n=", n, ")"))

  # Plot
  ggplot(res, aes(x=EGFR_mSig_class, y=EGFR_mSig_distance)) +
    geom_violin(trim=FALSE, fill="gray80") +
    geom_boxplot(width=0.15, outlier.shape=NA) +
    scale_x_discrete(
      labels = setNames(label_df$label, label_df$EGFR_mSig_class)
    ) +
    theme_bw(base_size=18) +
    labs(
      x = "",
      y = "EGFR mSig distance",
      title = ""
    )
}
