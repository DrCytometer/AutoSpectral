# spectral_heatmap.r

#' @title Spectral Heatmap
#'
#' @description
#' This function plots the spectral matrix as a heatmap using the specified
#' color palette, and saves it as a JPEG file.
#'
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_viridis_c theme_classic
#' @importFrom ggplot2 coord_fixed element_text labs ggsave theme
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate %>%
#'
#' @param spectra Matrix or dataframe containing spectral data
#' format: fluorophores x detectors.
#' @param title Optional prefix for the plot filename.
#' @param legend.label Character string that will appear on the heatmap legend.
#' @param plot.dir Optional output directory. Default is NULL, in which case
#' the spectra figure folder will be used.
#' @param color.palette Optional character string defining the viridis color
#' palette to be used for the fluorophore traces. Default is `viridis`. Options
#' are the viridis color options: `magma`, `inferno`, `plasma`, `viridis`,
#' `cividis`, `rocket`, `mako` and `turbo`.
#' @param save Logical, if `TRUE`, saves a JPEG file to the `plot.dir`.
#' Otherwise, the plot will simply be created in the Viewer.
#' @param show.legend Logical. If `TRUE`, figure legend will be included.
#' @param plot.width Width for the output plot. Default is `NULL`, in which case
#' the width will be scaled automatically based on the number of detectors in
#' `spectra`.
#' @param plot.height Height for the output plot. Default is `NULL`, in which
#' case the height will be scaled automatically based on the number of rows in
#' `spectra` plus a safety margin.
#'
#' @return Saves the heatmap plot as a JPEG file in the specified directory.
#'
#' @export

spectral.heatmap <- function(
    spectra,
    title = NULL,
    plot.dir = NULL,
    legend.label = "Intensity",
    color.palette = "viridis",
    save = TRUE,
    show.legend = TRUE,
    plot.width = NULL,
    plot.height = NULL
  ) {

  if ( !is.null( title ) )
    heatmap.filename <- paste( title, "spectral_heatmap.jpg" )
  else
    heatmap.filename <- "spectral_heatmap.jpg"

  if ( is.null( plot.dir ) )
    plot.dir <- getwd()

  heatmap.df <- data.frame( spectra, check.names = FALSE )

  if ( is.null( plot.width ) )
    plot.width <- max( ( ( ncol( heatmap.df ) - 1 ) / 64 * 12 ), 3 )
  if ( is.null( plot.height ) )
    plot.height <- 5 + round( nrow( heatmap.df ) / 8, 0 )

  row.levels <- rownames( heatmap.df )
  col.levels <- colnames( heatmap.df )
  heatmap.df$Fluorophore <- row.levels

  heatmap.long <- heatmap.df %>%
    tidyr::pivot_longer(
      cols = -Fluorophore,
      names_to = "Detector",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      Fluorophore = factor( Fluorophore, levels = rev( row.levels ) ),
      Detector = factor( Detector, levels = col.levels )
    )

  heatmap.plot <- ggplot( heatmap.long, aes( Detector, Fluorophore, fill = value ) ) +
    geom_tile() +
    theme_classic() +
    coord_fixed( ratio = 1 ) +
    theme( axis.text.x = element_text( angle = 45, hjust = 1 ) ) +
    labs( x = NULL, y = NULL, fill = legend.label )+
    scale_fill_viridis_c( option = color.palette )

  if ( !show.legend )
    heatmap.plot <- heatmap.plot + theme( legend.position = "none" )

  if ( save )
    ggsave(
      filename = file.path( plot.dir, heatmap.filename ),
      plot = heatmap.plot,
      width = plot.width,
      height = plot.height
    )
  else
    return( heatmap.plot )
}
