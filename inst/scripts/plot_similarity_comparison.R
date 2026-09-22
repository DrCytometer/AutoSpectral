# plot_similarity_comparison.R

#' @title Compare Two Cosine-Similarity Distributions
#'
#' @description
#' Overlays the density distributions of two sets of per-cell cosine
#' similarity values (for example, `per.cell.cosine.similarity()` versus
#' `compute.percell.variant.cosine.similarity()` for the same control), and
#' annotates the plot with the p-value of a two-sample Kolmogorov-Smirnov
#' test (`stats::ks.test()`) comparing the two distributions.
#'
#' @importFrom ggplot2 ggplot aes geom_density scale_fill_manual theme_minimal
#' @importFrom ggplot2 labs xlim theme element_text annotate
#' @importFrom stats ks.test setNames
#'
#' @param values.1 Numeric vector, the first set of cosine similarity values.
#' @param values.2 Numeric vector, the second set of cosine similarity
#' values.
#' @param label.1 Character, legend label for `values.1`. Default `"OLS"`.
#' @param label.2 Character, legend label for `values.2`. Default
#' `"PerCell"`.
#' @param color.1 Fill color for `values.1`'s density. Default
#' `"lightblue"`.
#' @param color.2 Fill color for `values.2`'s density. Default `"gold"`.
#' @param x.lim Numeric vector of length 2, x-axis limits. Default
#' `c(0.95, 1)`.
#' @param x.lab Label for the x-axis. Default `"Cosine Similarity"`.
#' @param title Optional plot title.
#' @param text.size Numeric, font size for the KS-test p-value annotation.
#' Default `3.5`.
#'
#' @return A ggplot object.
#'
#' @export

plot.similarity.comparison <- function(
    values.1,
    values.2,
    label.1 = "OLS",
    label.2 = "PerCell",
    color.1 = "lightblue",
    color.2 = "gold",
    x.lim = c( 0.95, 1 ),
    x.lab = "Cosine Similarity",
    title = NULL,
    text.size = 3.5
) {

  plot.data <- data.frame(
    value = c( values.1, values.2 ),
    Method = c(
      rep( label.1, length( values.1 ) ),
      rep( label.2, length( values.2 ) )
    )
  )

  ks.p <- stats::ks.test( values.1, values.2 )$p.value
  ks.label <- if ( ks.p < 2.2e-16 )
    "KS p < 2.2e-16" else sprintf( "KS p = %.2g", ks.p )

  density.plot <- ggplot( plot.data, aes( x = value, fill = Method ) ) +
    geom_density( color = "black", alpha = 0.5 ) +
    scale_fill_manual( values = stats::setNames(
      c( color.1, color.2 ), c( label.1, label.2 )
    ) ) +
    xlim( x.lim ) +
    labs( x = x.lab, y = NULL, title = title ) +
    theme_minimal() +
    theme( axis.text.x = element_text( angle = 45, hjust = 1 ) ) +
    annotate(
      "text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
      label = ks.label, size = text.size
    )

  density.plot
}
