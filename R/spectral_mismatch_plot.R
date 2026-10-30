# spectral_mismatch_plot.R

#' @title Spectral Mismatch Plot
#'
#' @description
#' Produces a two-panel diagnostic figure comparing two fluorophore spectra. The
#' upper panel overlays the normalised spectral traces for a reference and a test
#' spectrum across all detectors. The lower panel shows the per-detector spectral
#' difference (reference minus test) as a bar chart, making it easy to identify
#' channels where the spectra diverge.
#'
#' By default the combined figure is saved as a JPEG to `plot.dir` and returned
#' invisibly as a `ggplot` / cowplot object so it can be embedded in downstream
#' reports or further modified. Optional arguments normalise the spectra,
#' hide the detectors of some lasers, print the cosine similarity, thin the
#' detector labels, add a laser colour key and skip saving, so the figure can
#' be used directly as a manuscript panel.
#'
#' @param ref.spectrum A named numeric vector (or single-row matrix) of
#' normalised spectral intensities for the reference spectrum. Names must
#' correspond to detector channel names.
#' @param test.spectrum A named numeric vector (or single-row matrix) of
#' normalised spectral intensities for the test spectrum. Must have the same
#' names as `ref.spectrum`.
#' @param ref.label Character string used as the legend label for the reference
#' spectrum. If `NULL` (default), the `rownames()` of `ref.spectrum` is used
#' (`"Reference"` for a vector).
#' @param test.label Character string used as the legend label for the test
#' spectrum. If `NULL` (default), the `rownames()` of `test.spectrum` is used
#' (`"Test"` for a vector).
#' @param title Character string used as both the plot title and the JPEG
#' filename stem. Default is `"Spectral comparison"`.
#' @param bar.color Fill color for the spectral-difference bar chart.
#' Default is `"lightblue"`.
#' @param fluorophore.colors Character vector of length 2 giving the line colors
#' for the reference and test spectra respectively. Overrides `color.palette`
#' when non-`NULL`. Default is `c("black", "red")`.
#' @param plot.dir Directory to save the JPEG file. If `NULL` (default), the
#' current working directory is used. Created automatically if absent.
#' @param line.width Numeric. Line width for the spectral traces. Default `1`.
#' @param point.size Numeric. Point size for the spectral traces. Default `1`.
#' @param color.palette Optional viridis palette name (e.g., `"viridis"`,
#' `"plasma"`) used when `fluorophore.colors` is `NULL`. Ignored when
#' `fluorophore.colors` is set. Default is `NULL`.
#' @param show.legend Logical. Whether to display the color legend on the
#' upper panel. Default is `TRUE`.
#' @param plot.width Numeric. Width of the saved figure in inches. If `NULL`
#' (default), an automatic width is derived from the number of detectors
#' (approximately 12 inches per 64 channels, minimum 3 inches).
#' @param plot.height Numeric. Height of the saved figure in inches. If `NULL`
#' (default), an automatic height is derived from the number of spectra rows.
#' @param color.by.laser Logical, default `TRUE`. Colour the detector axis
#'   labels by excitation laser (from `cytometer_database.csv`). Falls back to
#'   black labels if the detectors do not all map to the instrument's detectors.
#' @param asp Optional AutoSpectral parameter list, used to identify the
#'   cytometer for the laser of each detector (`color.by.laser`,
#'   `exclude.lasers`, `show.laser.key`). If `NULL` (default), the cytometer is
#'   identified from the detector names.
#' @param normalization How both spectra are scaled before plotting and before
#'   the difference is taken, always over all detectors (before
#'   `exclude.lasers`): `"none"` (default; as supplied), `"max"` (peak = 1) or
#'   `"L2"` (unit Euclidean norm).
#' @param exclude.lasers Optional character vector of lasers (as named in
#'   `cytometer_database.csv`, e.g. `"DeepUV"`) whose detectors are left out of
#'   the display. Normalisation, and by default the cosine similarity, still
#'   use all detectors. Default `NULL`.
#' @param show.cosine Logical, default `FALSE`. Print the cosine similarity of
#'   the two spectra in the upper panel.
#' @param cosine.on.plotted Logical, default `FALSE`. Calculate the cosine
#'   similarity on the displayed detectors only (otherwise on all detectors).
#' @param cosine.digits Integer, default `3`. Decimal places of the printed
#'   cosine similarity.
#' @param label.every Show every n-th detector label. Default `1` (every
#'   label); `NULL` chooses n so that at most about 70 labels are shown.
#' @param label.angle Numeric, default `45`. Angle of the detector labels
#'   (`90` for vertical labels).
#' @param trace.labels Logical, default `TRUE`. Show detector labels under the
#'   upper panel as well as under the difference panel. `FALSE` labels the
#'   difference panel only.
#' @param y.label Optional title of the upper panel's y axis. Default `NULL`:
#'   `"Intensity (L2-normalized)"` with `normalization = "L2"`, otherwise
#'   `"Normalized Intensity"`.
#' @param diff.label Optional title of the difference panel's y axis. Default
#'   `NULL`: `"Spectral Difference (ref - test)"`.
#' @param legend.position Position of the trace legend. Default `"bottom"`.
#' @param legend.title Title of the trace legend. Default `"Fluorophore"`;
#'   `NULL` for none.
#' @param bar.line.width Numeric, default `0.5`. Outline width of the
#'   difference bars (thinner outlines suit instruments with many detectors).
#' @param show.laser.key Logical, default `FALSE`. Add a one-row key of the
#'   laser colours under the figure.
#' @param base.size Numeric, default `11`. Base font size (points).
#' @param axis.text.size Numeric. Detector label size in points. Default
#'   `NULL`: `0.8 * base.size`.
#' @param save Logical, default `TRUE`. Save the figure as
#'   `<plot.dir>/<title>.jpg`. `FALSE` only returns it.
#'
#' @return The combined cowplot object is returned invisibly, with the cosine
#' similarity as `attr(x, "cosine.similarity")` and the lasers of the
#' displayed detectors as `attr(x, "lasers")` (`NULL` if the lasers could not
#' be assigned). With `save = TRUE` (default) the figure is also saved to disk
#' as `<plot.dir>/<title>.jpg`.
#'
#' @importFrom ggplot2 ggplot aes geom_path geom_point geom_bar labs
#' @importFrom ggplot2 theme_minimal theme element_text scale_color_manual
#' @importFrom ggplot2 scale_color_viridis_d
#' @importFrom cowplot plot_grid
#' @importFrom ragg agg_jpeg
#'
#' @export

spectral.mismatch.plot <- function(
    ref.spectrum,
    test.spectrum,
    ref.label          = NULL,
    test.label         = NULL,
    title              = "Spectral comparison",
    bar.color          = "lightblue",
    fluorophore.colors = c( "black", "red" ),
    plot.dir           = NULL,
    line.width         = 1,
    point.size         = 1,
    color.palette      = NULL,
    show.legend        = TRUE,
    plot.width         = NULL,
    plot.height        = NULL,
    color.by.laser     = TRUE,
    asp                = NULL,
    normalization      = c( "none", "max", "L2" ),
    exclude.lasers     = NULL,
    show.cosine        = FALSE,
    cosine.on.plotted  = FALSE,
    cosine.digits      = 3,
    label.every        = 1,
    label.angle        = 45,
    trace.labels       = TRUE,
    y.label            = NULL,
    diff.label         = NULL,
    legend.position    = "bottom",
    legend.title       = "Fluorophore",
    bar.line.width     = 0.5,
    show.laser.key     = FALSE,
    base.size          = 11,
    axis.text.size     = NULL,
    save               = TRUE
  ) {

  normalization <- match.arg( normalization )

  # --- input validation ---
  if ( !is.numeric( as.matrix( ref.spectrum ) ) || !is.numeric( as.matrix( test.spectrum ) ) ) {
    stop( "`ref.spectrum` and `test.spectrum` must be numeric.", call. = FALSE )
  }

  ref.spectrum  <- .mismatch.as.row( ref.spectrum, "Reference" )
  test.spectrum <- .mismatch.as.row( test.spectrum, "Test" )

  if ( !identical( colnames( ref.spectrum ), colnames( test.spectrum ) ) ) {
    stop(
      "`ref.spectrum` and `test.spectrum` must have identical column names.",
      call. = FALSE
    )
  }

  ref.name  <- if ( !is.null( ref.label ) )  ref.label  else rownames( ref.spectrum )
  test.name <- if ( !is.null( test.label ) ) test.label else rownames( test.spectrum )

  # --- normalisation, over all detectors ---
  normalize <- function( v ) switch(
    normalization,
    none = v,
    max  = v / max( v ),
    L2   = v / sqrt( sum( v^2 ) )
  )
  ref.vals  <- normalize( stats::setNames( as.numeric( ref.spectrum[ 1, ] ), colnames( ref.spectrum ) ) )
  test.vals <- normalize( stats::setNames( as.numeric( test.spectrum[ 1, ] ), colnames( test.spectrum ) ) )
  detectors <- names( ref.vals )

  # --- lasers: display exclusion and key ---
  lasers <- NULL
  if ( !is.null( exclude.lasers ) || show.laser.key || isTRUE( color.by.laser ) )
    lasers <- .detector.lasers( detectors, asp = asp, verbose = FALSE )
  show <- rep( TRUE, length( detectors ) )
  if ( !is.null( exclude.lasers ) ) {
    if ( is.null( lasers ) ) {
      warning( "Lasers could not be assigned to the detectors; `exclude.lasers` ignored.",
               call. = FALSE )
    } else {
      show <- !( lasers %in% exclude.lasers )
    }
  }

  # --- cosine similarity ---
  cos.idx <- if ( cosine.on.plotted ) show else rep( TRUE, length( detectors ) )
  cos.sim <- cosine.similarity( rbind( ref.vals[ cos.idx ], test.vals[ cos.idx ] ) )[ 1, 2 ]

  pivot.cols <- detectors[ show ]
  shown.lasers <- if ( is.null( lasers ) ) NULL else unique( unname( lasers[ show ] ) )

  # --- detector labels ---
  if ( is.null( label.every ) ) label.every <- max( 1L, ceiling( length( pivot.cols ) / 70 ) )
  label.breaks <- pivot.cols[ seq( 1, length( pivot.cols ), by = label.every ) ]
  if ( is.null( axis.text.size ) ) axis.text.size <- 0.8 * base.size
  label.vjust <- if ( label.angle == 90 ) 0.5 else 1

  # --- output directory ---
  if ( save ) {
    if ( is.null( plot.dir ) ) plot.dir <- getwd()
    if ( !dir.exists( plot.dir ) ) dir.create( plot.dir, recursive = TRUE )
  }

  # --- automatic figure dimensions ---
  if ( is.null( plot.width ) )
    plot.width <- max( ( ( length( pivot.cols ) ) / 64 * 12 ), 3 )

  if ( is.null( plot.height ) )
    plot.height <- 5   # 5 + round( 2 spectra / 8 ) * 2

  # --- long format (base R, no tidyr dependency) ---
  spectra.long <- data.frame(
    Fluorophore      = rep( c( ref.name, test.name ), times = length( pivot.cols ) ),
    Detector         = rep( pivot.cols, each = 2 ),
    Intensity        = as.vector( rbind( ref.vals[ pivot.cols ], test.vals[ pivot.cols ] ) ),
    stringsAsFactors = FALSE
  )

  spectra.long$Detector <- factor(
    spectra.long$Detector,
    levels  = pivot.cols,
    ordered = TRUE
  )
  spectra.long$Fluorophore <- factor( spectra.long$Fluorophore,
                                      levels = unique( c( ref.name, test.name ) ) )

  # --- laser-coloured detector labels ---
  laser.labels <- .laser.axis.labels(
    x = factor( label.breaks, levels = pivot.cols, ordered = TRUE ), labels = label.breaks,
    asp = asp, color.by.laser = color.by.laser, angle = label.angle, hjust = 1,
    vjust = label.vjust, size = axis.text.size
  )
  label.clip <- if ( laser.labels$mapped ) "off" else "on"

  if ( is.null( y.label ) )
    y.label <- if ( normalization == "L2" ) "Intensity (L2-normalized)" else "Normalized Intensity"
  if ( is.null( diff.label ) ) diff.label <- "Spectral Difference (ref − test)"

  x.text <- ggplot2::element_text( angle = label.angle, hjust = 1, vjust = label.vjust,
                                   size = axis.text.size )

  # --- upper panel: spectral traces ---
  spectra.plot <- ggplot2::ggplot(
    spectra.long,
    ggplot2::aes(
      x     = Detector,
      y     = Intensity,
      group = Fluorophore,
      color = Fluorophore
    )
  ) +
    ggplot2::geom_path( linewidth = line.width ) +
    ggplot2::geom_point( size = point.size ) +
    ggplot2::scale_x_discrete( breaks = label.breaks ) +
    ggplot2::labs(
      title = title,
      x     = if ( trace.labels ) "Detector" else NULL,
      y     = y.label,
      color = legend.title
    ) +
    ggplot2::theme_minimal( base_size = base.size ) +
    ggplot2::theme(
      axis.text.x    = x.text,
      legend.position = legend.position
    ) +
    ggplot2::coord_cartesian( clip = label.clip )

  spectra.plot <- if ( trace.labels ) {
    spectra.plot + laser.labels$components
  } else {
    spectra.plot + ggplot2::theme( axis.text.x = ggplot2::element_blank() )
  }

  if ( show.cosine )
    spectra.plot <- spectra.plot + ggplot2::annotate(
      "text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.6, size = base.size * 0.32,
      label = sprintf( "cosine similarity = %s",
                       formatC( cos.sim, digits = cosine.digits, format = "f" ) )
    )

  # apply color palette
  if ( !is.null( fluorophore.colors ) ) {
    spectra.plot <- spectra.plot +
      ggplot2::scale_color_manual( values = fluorophore.colors )
  } else if ( !is.null( color.palette ) ) {
    spectra.plot <- spectra.plot +
      ggplot2::scale_color_viridis_d( option = color.palette )
  }

  if ( !show.legend ) {
    spectra.plot <- spectra.plot +
      ggplot2::theme( legend.position = "none" )
  }

  # --- lower panel: spectral difference ---
  # difference is always ref - test (after normalisation)
  delta.long <- data.frame(
    Detector         = pivot.cols,
    Intensity        = unname( ref.vals[ pivot.cols ] - test.vals[ pivot.cols ] ),
    stringsAsFactors = FALSE
  )
  delta.long$Detector <- factor(
    delta.long$Detector,
    levels  = pivot.cols,
    ordered = TRUE
  )

  delta.plot <- ggplot2::ggplot(
    delta.long,
    ggplot2::aes( x = Detector, y = Intensity )
  ) +
    ggplot2::geom_bar(
      stat  = "identity",
      fill  = bar.color,
      color = "black",
      alpha = 0.7,
      linewidth = bar.line.width
    ) +
    ggplot2::scale_x_discrete( breaks = label.breaks ) +
    ggplot2::labs( x = "Detector", y = diff.label ) +
    ggplot2::theme_minimal( base_size = base.size ) +
    ggplot2::theme(
      axis.text.x = x.text
    ) +
    ggplot2::coord_cartesian( clip = label.clip ) +
    laser.labels$components

  # --- combine panels ---
  combined.plot <- cowplot::plot_grid(
    spectra.plot,
    delta.plot,
    ncol        = 1,
    align       = "v",
    axis        = "lr",
    rel_heights = c( 2, 1 )
  )

  if ( show.laser.key ) {
    if ( is.null( shown.lasers ) ) {
      warning( "Lasers could not be assigned to the detectors; no laser key drawn.",
               call. = FALSE )
    } else {
      combined.plot <- cowplot::plot_grid(
        combined.plot, .laser.key.legend( shown.lasers, base.size = base.size ),
        ncol = 1, rel_heights = c( 1, 0.07 )
      )
    }
  }

  # --- save ---
  if ( save )
    ggplot2::ggsave(
      file.path( plot.dir, sprintf( "%s.jpg", title ) ),
      combined.plot,
      device    = ragg::agg_jpeg,
      width     = plot.width,
      height    = plot.height,
      limitsize = FALSE
    )

  attr( combined.plot, "cosine.similarity" ) <- cos.sim
  attr( combined.plot, "lasers" ) <- shown.lasers

  invisible( combined.plot )
}

# One-row matrix with detector column names. A named vector becomes a row
# named `default.name`; for a matrix or data frame the first row is used.
.mismatch.as.row <- function( x, default.name ) {
  if ( is.null( dim( x ) ) ) {
    m <- matrix( as.numeric( x ), nrow = 1, dimnames = list( default.name, names( x ) ) )
  } else {
    m <- as.matrix( x )[ 1, , drop = FALSE ]
    if ( is.null( rownames( m ) ) ) rownames( m ) <- default.name
  }
  m
}

#' One-row legend of laser colours
#'
#' @param lasers Laser names, in display order.
#' @param laser.colors Named colour vector (default `.laser.colors`).
#' @param base.size Base font size.
#' @return A legend grob for `cowplot::plot_grid()`.
#' @noRd
.laser.key.legend <- function( lasers, laser.colors = .laser.colors, base.size = 11 ) {
  lasers <- unique( lasers )
  key.cols <- laser.colors[ lasers ]
  key.cols[ is.na( key.cols ) ] <- "grey50"
  key.df <- data.frame( Laser = factor( lasers, levels = lasers ), x = 1, y = 1 )
  key.plot <- ggplot2::ggplot( key.df, ggplot2::aes( x, y, colour = Laser ) ) +
    ggplot2::geom_point( size = 2.5 ) +
    ggplot2::scale_colour_manual( values = stats::setNames( unname( key.cols ), lasers ),
                                  name = "Laser" ) +
    ggplot2::guides( colour = ggplot2::guide_legend( nrow = 1 ) ) +
    ggplot2::theme_void( base_size = base.size ) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.margin = ggplot2::margin( 0, 0, 0, 0 )
    )
  cowplot::get_legend( key.plot )
}
