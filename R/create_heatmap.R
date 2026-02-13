# create_heatmap.r


#' @title Create Heatmap Plot
#'
#' @description
#' This function plots a matrix as a heatmap and saves it as a JPEG file.
#'
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_viridis_c theme_minimal
#' @importFrom ggplot2 coord_fixed element_text labs ggsave
#' @importFrom ragg agg_jpeg
#'
#' @param matrix Matrix or dataframe containing spectral data.
#' @param number.labels Logical indicating whether to add number labels to
#' the heatmap. Default is `FALSE`.
#' @param title Optional prefix for the plot filename. Default is `NULL`,
#' in which case the file will just be called `heatmap.jpg`
#' @param legend.label Character string that will appear on the heatmap legend.
#' Default is `heatmap`
#' @param triangular Logical. Plot the lower triangle of the matrix only,
#' diagonal included. Default is `FALSE`.
#' @param plot.dir Optional output directory. Default is `NULL`, in which case
#' the working directory will be used.
#' @param fixed.scale Logical, determines whether to use an externally supplied
#' fixed scale (min and max) for the heatmap color scale. Useful for putting
#' multiple plots on the same scale. Default is `FALSE`
#' @param scale.min Optional numeric. Minimum for the fixed color scale.
#' Default is `NULL`, for no fixed scaling.
#' @param scale.max Optional numeric. Maximum for the fixed color scale.
#' Default is `NULL`, for no fixed scaling.
#' @param color.palette Optional character string defining the viridis color
#' palette to be used for the fluorophore traces. Default is `viridis`. Options
#' are the viridis color options: `magma`, `inferno`, `plasma`, `viridis`,
#' `cividis`, `rocket`, `mako` and `turbo`.
#' @param show.legend Logical. If `TRUE`, figure legend will be included.
#' @param figure.width Numeric. Width of the heatmap figure. Default is `8`.
#' @param figure.height Numeric. Height of the heatmap figure. Default is `6`.
#' @param save Logical, if `TRUE`, saves a JPEG file to the `output.dir`.
#' Otherwise, the plot will simply be created in the Viewer.
#'
#' @return Saves the heatmap plot as a JPEG file and the SSM data as a CSV file
#' in the specified directory.
#'
#' @export

create.heatmap <- function(
    matrix,
    number.labels = FALSE,
    title = NULL,
    legend.label = "heatmap",
    triangular = FALSE,
    plot.dir = NULL,
    fixed.scale = FALSE,
    scale.min = NULL,
    scale.max = NULL,
    color.palette = "viridis",
    show.legend = TRUE,
    figure.width = 8,
    figure.height = 6,
    save = TRUE
  ) {

  if ( !is.null( title ) )
    heatmap.filename <- paste( title, "heatmap.jpg" )
  else
    heatmap.filename <- "heatmap.jpg"

  if ( is.null( plot.dir ) )
    plot.dir <- getwd()

  # rearrange data
  heatmap.df <- as.data.frame( matrix, check.names = FALSE )
  row.levels <- rownames( heatmap.df )
  col.levels <- colnames( heatmap.df )
  heatmap.df$Fluor1 <- row.levels

  # pivot longer
  heatmap.long <- data.frame(
    Fluor1 = rep( heatmap.df$Fluor1, times = ncol( heatmap.df ) - 1 ),
    Fluor2 = rep( col.levels, each = nrow( heatmap.df ) ),
    value  = as.vector( as.matrix( heatmap.df[, col.levels] ) ),
    stringsAsFactors = FALSE
  )

  # convert to factors with original levels
  heatmap.long$Fluor1 <- factor( heatmap.long$Fluor1, levels = row.levels )
  heatmap.long$Fluor2 <- factor( heatmap.long$Fluor2, levels = col.levels )

  # triangular or full heatmap, depending on user-input
  if ( triangular ) {
    # keep only upper triangle (Fluor1 <= Fluor2)
    keep <- as.integer( heatmap.long$Fluor1 ) <= as.integer( heatmap.long$Fluor2 )
    heatmap.long <- heatmap.long[ keep, ]

    # reverse Fluor2 levels for plotting
    heatmap.long$Fluor2 <- factor(
      heatmap.long$Fluor2,
      levels = rev( col.levels )
    )
  } else {
    # full heatmap: reverse Fluor1 for plotting
    heatmap.long$Fluor1 <- factor(
      heatmap.long$Fluor1,
      levels = rev( row.levels )
    )
  }

  # plot and save
  heatmap.plot <- ggplot(
    heatmap.long,
    aes(
      Fluor1,
      Fluor2,
      fill = value
      )
    ) +
    geom_tile() +
    theme_classic() +
    coord_fixed( ratio = 1 ) +
    theme( axis.text.x = element_text( angle = 45, hjust = 1 ) ) +
    labs( x = NULL, y = NULL, fill = legend.label )

  if ( fixed.scale ) {
    heatmap.plot <- heatmap.plot +
      scale_fill_viridis_c(
        option = color.palette,
        limits = c( scale.min, scale.max ) )
  } else {
    heatmap.plot <- heatmap.plot +
      scale_fill_viridis_c( option = color.palette )
  }

  if ( number.labels ) {
    heatmap.long$label <- round( heatmap.long$value, 2 )

    heatmap.plot <- heatmap.plot +
      geom_text(
        aes( x = Fluor1, y = Fluor2, label = label ),
        data = heatmap.long,
        color = "white",
        size = 3
      )
  }

  if ( !show.legend )
    heatmap.plot <- heatmap.plot + theme( legend.position = "none" )

  if ( save )
    ggsave(
      filename = file.path( plot.dir, heatmap.filename ),
      plot = heatmap.plot,
      device = ragg::agg_jpeg,
      width = figure.width,
      height = figure.height,
      limitsize = FALSE
    )
  else
    return( heatmap.plot )
}
