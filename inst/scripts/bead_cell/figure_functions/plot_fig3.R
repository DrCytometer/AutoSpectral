# plot_fig3.R
#
# Plotting helpers for Figure 3:
#   plot.difference.heatmap()  detectors x rows heatmap of spectral differences
#                              with detector labels coloured by laser and a
#                              side column of per-row annotation (e.g. cosine
#                              similarity with the median AF spectrum)
#   plot.af.collinearity()     fraction of each spectral delta lying in the AF
#                              subspace, with bootstrap CI, against the same
#                              quantity for every dye in the reference library

.asp.internal <- function( name ) utils::getFromNamespace( name, "AutoSpectral" )

#' Heatmap of spectral differences
#'
#' @param mat Numeric matrix, rows = spectra (row names are the labels, in
#'   display order top to bottom), columns = detectors (display order).
#' @param group Optional character vector (one per row): a separator line is
#'   drawn where it changes.
#' @param annot Optional named numeric vector (names = row names) shown as a
#'   coloured side column with the value printed.
#' @param annot.text Optional character vector (names = row names) printed in
#'   the side column instead of the formatted value (e.g. to add P values).
#' @param annot.title Title of the side column.
#' @param annot.limits Colour limits of the side column.
#' @param scale.factor Multiply the differences by this (100 = percent of peak).
#' @param limits Symmetric fill limit; NULL = largest |value|.
#' @param fill.label Legend title.
#' @param asp AutoSpectral parameter list (laser colours of detector labels).
#' @param label.size Detector label size (pt).
#' @param rel.width.annot Relative width of the side column.
#' @return A cowplot object.
plot.difference.heatmap <- function(
    mat,
    group = NULL,
    annot = NULL,
    annot.text = NULL,
    annot.title = "cos(Δ, AF)",
    annot.limits = c( -1, 1 ),
    scale.factor = 100,
    limits = NULL,
    fill.label = "Δ (% of peak)",
    title = NULL,
    asp = NULL,
    base.size = 9,
    label.size = 6,
    rel.width.annot = 0.09,
    low = "#2166AC", mid = "white", high = "#B2182B"
) {

  rows <- rownames( mat )
  dets <- colnames( mat )
  vals <- mat * scale.factor
  if ( is.null( limits ) ) limits <- max( abs( vals ), na.rm = TRUE )
  if ( !is.finite( limits ) || limits == 0 ) limits <- 1

  long <- data.frame(
    Row      = factor( rep( rows, times = length( dets ) ), levels = rev( rows ) ),
    Detector = factor( rep( dets, each = length( rows ) ), levels = dets ),
    value    = as.vector( vals )
  )

  laser.labels <- .asp.internal( ".laser.axis.labels" )(
    x = factor( dets, levels = dets ), labels = dets, asp = asp,
    color.by.laser = TRUE, angle = 90, hjust = 1, vjust = 0.5, size = label.size,
    verbose = FALSE
  )

  sep <- NULL
  if ( !is.null( group ) ) {
    change <- which( group[ -1 ] != group[ -length( group ) ] )
    if ( length( change ) > 0 ) sep <- length( rows ) - change + 0.5
  }

  p <- ggplot2::ggplot( long, ggplot2::aes( Detector, Row, fill = value ) ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(
      low = low, mid = mid, high = high, midpoint = 0,
      limits = c( -limits, limits ), oob = scales::squish, name = fill.label
    ) +
    ggplot2::scale_x_discrete( expand = c( 0, 0 ) ) +
    ggplot2::scale_y_discrete( expand = c( 0, 0 ) ) +
    ggplot2::coord_cartesian( clip = "off" ) +
    ggplot2::labs( x = NULL, y = NULL, title = title ) +
    ggplot2::theme_classic( base_size = base.size ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text( angle = 90, hjust = 1, vjust = 0.5, size = label.size ),
      axis.text.y = ggplot2::element_text( size = base.size - 1 ),
      axis.ticks.x = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.key.height = grid::unit( 0.25, "cm" ),
      legend.key.width = grid::unit( 1.2, "cm" ),
      plot.title = ggplot2::element_text( face = "bold", size = base.size + 1 )
    )
  if ( !is.null( sep ) )
    p <- p + ggplot2::geom_hline( yintercept = sep, linewidth = 0.4, colour = "grey20" )
  p <- p + laser.labels$components

  if ( is.null( annot ) ) return( p )

  a <- data.frame(
    Row = factor( rows, levels = rev( rows ) ),
    x = annot.title,
    value = unname( annot[ rows ] ),
    label = if ( is.null( annot.text ) )
      ifelse( is.na( annot[ rows ] ), "", formatC( annot[ rows ], digits = 2, format = "f" ) ) else
      unname( annot.text[ rows ] )
  )
  pa <- ggplot2::ggplot( a, ggplot2::aes( x, Row, fill = value ) ) +
    ggplot2::geom_tile( colour = "white" ) +
    ggplot2::geom_text( ggplot2::aes( label = label ), size = ( base.size - 2.5 ) / ggplot2::.pt ) +
    ggplot2::scale_fill_gradient2(
      low = "#5E3C99", mid = "white", high = "#E66101", midpoint = 0,
      limits = annot.limits, oob = scales::squish, na.value = "grey90", guide = "none"
    ) +
    ggplot2::scale_x_discrete( expand = c( 0, 0 ), position = "top" ) +
    ggplot2::scale_y_discrete( expand = c( 0, 0 ) ) +
    ggplot2::labs( x = NULL, y = NULL ) +
    ggplot2::theme_classic( base_size = base.size ) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      axis.text.x.top = ggplot2::element_text( size = base.size - 2 ),
      legend.position = "bottom"
    )
  if ( !is.null( sep ) )
    pa <- pa + ggplot2::geom_hline( yintercept = sep, linewidth = 0.4, colour = "grey20" )

  cowplot::plot_grid( p, pa, nrow = 1, align = "h", axis = "tb",
                      rel_widths = c( 1, rel.width.annot ) )
}

#' Fraction of each delta in the AF subspace, against the reference library
#'
#' @param stats data.frame with columns fluor, case, r2.af, r2.af.lo, r2.af.hi,
#'   cos.af, p.value.
#' @param null data.frame with columns fluor, r2 (library dyes).
#' @param k,var.explained For the subtitle.
plot.af.collinearity <- function( stats, null, k, var.explained, base.size = 9,
                                  title = "Wrong-background error lies in the AF subspace" ) {

  stats$label <- sprintf( "cos %.2f", stats$cos.af )
  stats$Sign <- ifelse( stats$cos.af >= 0, "Excess AF (cos > 0)", "AF deficit (cos < 0)" )
  null$case <- "Library dyes"
  lev <- c( "Library dyes", unique( stats$case ) )
  stats$case <- factor( stats$case, levels = lev )
  null$case <- factor( null$case, levels = lev )
  stats$fluor <- factor( stats$fluor, levels = unique( stats$fluor ) )
  null$fluor <- factor( null$fluor, levels = levels( stats$fluor ) )

  ggplot2::ggplot() +
    ggplot2::geom_boxplot( data = null, ggplot2::aes( case, r2 ), outlier.shape = NA,
                           fill = "grey90", colour = "grey55", width = 0.6 ) +
    ggplot2::geom_jitter( data = null, ggplot2::aes( case, r2 ), width = 0.15, height = 0,
                          size = 0.4, colour = "grey55", alpha = 0.6 ) +
    ggplot2::geom_errorbar( data = stats, ggplot2::aes( case, ymin = r2.af.lo, ymax = r2.af.hi ),
                            width = 0.25 ) +
    ggplot2::geom_point( data = stats, ggplot2::aes( case, r2.af, fill = Sign ),
                         shape = 21, size = 2.4 ) +
    ggplot2::geom_text( data = stats, ggplot2::aes( case, r2.af.lo, label = label ),
                        vjust = 1.6, size = ( base.size - 3 ) / ggplot2::.pt ) +
    ggplot2::scale_fill_manual( values = c( "Excess AF (cos > 0)" = "#E66101",
                                            "AF deficit (cos < 0)" = "#5E3C99" ), name = NULL,
                                drop = FALSE ) +
    ggplot2::scale_y_continuous( limits = c( -0.08, 1 ), breaks = seq( 0, 1, 0.25 ) ) +
    ggplot2::facet_wrap( ~ fluor, nrow = 1 ) +
    ggplot2::labs(
      x = NULL,
      y = expression( "Fraction of " * "||" * Delta * "||"^2 * " in AF subspace" ),
      title = title,
      subtitle = sprintf( "AF subspace: first %d uncentred PCs of unstained splenocytes (%s of AF variance)",
                          k, paste0( formatC( 100 * sum( var.explained ), digits = 1, format = "f" ), "%" ) )
    ) +
    ggplot2::theme_bw( base_size = base.size ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text( angle = 45, hjust = 1 ),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      plot.title = ggplot2::element_text( face = "bold", size = base.size + 1 ),
      plot.subtitle = ggplot2::element_text( size = base.size - 1 )
    )
}
