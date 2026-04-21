# plot_unmixed_mxn.R

# =============================================================================
# Internal helpers (not exported)
# =============================================================================

# -----------------------------------------------------------------------------
# .read.unmixed.input
#   Accepts either a path to an FCS file or an already-loaded matrix/data.frame.
#   Optionally subsets columns to `channels`.  Returns a plain matrix.
# -----------------------------------------------------------------------------
.read.unmixed.input <- function( unmixed.data, channels = NULL ) {

  if ( is.character( unmixed.data ) ) {
    if ( length( unmixed.data ) != 1L )
      stop( "`unmixed.data` must be a single file path or a matrix/data.frame.",
            call. = FALSE )
    if ( !file.exists( unmixed.data ) )
      stop( paste0( "FCS file not found: '", unmixed.data, "'" ), call. = FALSE )
    mat <- AutoSpectral::readFCS( unmixed.data )
  } else if ( is.matrix( unmixed.data ) || is.data.frame( unmixed.data ) ) {
    mat <- as.matrix( unmixed.data )
  } else {
    stop(
      "`unmixed.data` must be a file path, matrix, or data.frame.",
      call. = FALSE
    )
  }

  if ( !is.null( channels ) ) {
    missing.ch <- setdiff( channels, colnames( mat ) )
    if ( length( missing.ch ) > 0 )
      stop(
        paste0(
          "The following requested channels are absent from the data: ",
          paste( missing.ch, collapse = ", " )
        ),
        call. = FALSE
      )
    mat <- mat[ , channels, drop = FALSE ]
  }

  mat
}


# -----------------------------------------------------------------------------
# .make.biexp.transforms
#   Builds a pair of flowjo_biexp transform functions for x and y axes,
#   replicating the excess-width-basis workaround from create.biplot().
# -----------------------------------------------------------------------------
.make.biexp.transforms <- function(
    asp,
    x.min, x.max,
    y.min, y.max,
    x.width.basis,
    y.width.basis
  ) {

  # x transform
  if ( x.width.basis < -1000 ) {
    x.pos.log <- pmax( log10( x.max ) - 1 - log10( abs( x.width.basis + 1000 ) ), 2 )
  } else {
    x.pos.log <- log10( x.max ) - 1
  }

  # y transform
  if ( y.width.basis < -1000 ) {
    y.pos.log <- pmax( log10( y.max ) - 1 - log10( abs( y.width.basis + 1000 ) ), 2 )
  } else {
    y.pos.log <- log10( y.max ) - 1
  }

  tx <- flowWorkspace::flowjo_biexp(
    channelRange = asp$default.transformation.param$length,
    maxValue     = x.max,
    pos          = x.pos.log,
    neg          = asp$default.transformation.param$neg,
    widthBasis   = x.width.basis,
    inverse      = FALSE
  )

  ty <- flowWorkspace::flowjo_biexp(
    channelRange = asp$default.transformation.param$length,
    maxValue     = y.max,
    pos          = y.pos.log,
    neg          = asp$default.transformation.param$neg,
    widthBasis   = y.width.basis,
    inverse      = FALSE
  )

  list( x = tx, y = ty )
}


# -----------------------------------------------------------------------------
# .axis.setup
#   Returns named list: breaks (data space), limits (data space), labels.
# -----------------------------------------------------------------------------
.axis.setup <- function( asp, dim.min, dim.max ) {
  brks   <- asp$ribbon.breaks[ asp$ribbon.breaks < dim.max ]
  labels <- sapply( brks, function( v ) {
    if ( v == 0 ) "0" else parse( text = paste0( "10^", log10( abs( v ) ) ) )
  } )
  list( breaks = brks, limits = c( dim.min, dim.max ), labels = labels )
}


# -----------------------------------------------------------------------------
# .make.panel
#   Produces a single ggplot panel (no save).  x.dim is the target fluorophore;
#   y.dim is the channel on the y-axis.  Dispatches to geom_hex or
#   geom_scattermore + stat_density_2d according to `use.hex`.
# -----------------------------------------------------------------------------
.make.panel <- function(
    plot.data,          # already-downsampled matrix / data.frame
    x.dim,
    y.dim,
    asp,
    tx,                 # pre-built biexp transform for x
    ty,                 # pre-built biexp transform for y
    x.axis,             # output of .axis.setup() for x
    y.axis,             # output of .axis.setup() for y
    use.hex    = TRUE,
    hex.bins   = 64,
    color.palette = "viridis",
    strip.axes = FALSE  # if TRUE, suppress axis text/title (for interior panels)
  ) {

  viridis.options <- c(
    "magma", "inferno", "plasma", "viridis",
    "cividis", "rocket", "mako", "turbo"
  )

  df <- data.frame(
    x.trans = tx( plot.data[ , x.dim ] ),
    y.trans = ty( plot.data[ , y.dim ] )
  )

  panel <- ggplot2::ggplot( df, ggplot2::aes( x = x.trans, y = y.trans ) )

  if ( use.hex ) {
    panel <- panel +
      ggplot2::geom_hex( bins = hex.bins, na.rm = TRUE )

    if ( color.palette %in% viridis.options ) {
      panel <- panel +
        ggplot2::scale_fill_viridis_c( option = color.palette )
    } else {
      panel <- panel +
        ggplot2::scale_fill_gradientn( colors = asp$density.palette.base.color )
    }
  } else {
    panel <- panel +
      scattermore::geom_scattermore(
        pointsize = asp$figure.gate.point.size,
        color     = "black",
        alpha     = 1,
        na.rm     = TRUE
      ) +
      ggplot2::stat_density_2d(
        ggplot2::aes( fill = ggplot2::after_stat( level ) ),
        geom  = "polygon",
        na.rm = TRUE
      )

    if ( color.palette %in% viridis.options ) {
      panel <- panel +
        ggplot2::scale_fill_viridis_c( option = color.palette )
    } else {
      panel <- panel +
        ggplot2::scale_fill_gradientn( colors = asp$density.palette.base.color )
    }
  }

  panel <- panel +
    ggplot2::scale_x_continuous(
      name   = x.dim,
      breaks = tx( x.axis$breaks ),
      limits = tx( x.axis$limits ),
      labels = x.axis$labels
    ) +
    ggplot2::scale_y_continuous(
      name   = y.dim,
      breaks = ty( y.axis$breaks ),
      limits = ty( y.axis$limits ),
      labels = y.axis$labels
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.margin  = ggplot2::margin(
        asp$figure.margin, asp$figure.margin,
        asp$figure.margin, asp$figure.margin
      ),
      legend.position   = "none",
      axis.ticks        = ggplot2::element_line(
        linewidth = asp$figure.panel.line.size
      ),
      axis.text         = ggplot2::element_text(
        size = asp$figure.axis.text.size
      ),
      axis.title        = ggplot2::element_text(
        size = asp$figure.axis.title.size
      ),
      panel.border      = ggplot2::element_rect(
        fill = NA, linewidth = asp$figure.panel.line.size
      ),
      panel.grid.major  = ggplot2::element_blank(),
      panel.grid.minor  = ggplot2::element_blank()
    )

  if ( strip.axes ) {
    panel <- panel +
      ggplot2::theme(
        axis.text  = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank()
      )
  }

  panel
}


# -----------------------------------------------------------------------------
# .make.diagonal.hist
#   1-D density histogram for a single channel, used on the n×n diagonal.
# -----------------------------------------------------------------------------
.make.diagonal.hist <- function(
    plot.data,
    dim,
    asp,
    tx,
    x.axis,
    color.palette = "viridis",
    fill.color    = "steelblue"
  ) {

  df <- data.frame( x.trans = tx( plot.data[ , dim ] ) )

  ggplot2::ggplot( df, ggplot2::aes( x = x.trans ) ) +
    ggplot2::geom_histogram(
      fill     = fill.color,
      color    = "white",
      bins     = 60,
      na.rm    = TRUE
    ) +
    ggplot2::scale_x_continuous(
      name   = dim,
      breaks = tx( x.axis$breaks ),
      limits = tx( x.axis$limits ),
      labels = x.axis$labels
    ) +
    ggplot2::labs( y = NULL ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.margin       = ggplot2::margin(
        asp$figure.margin, asp$figure.margin,
        asp$figure.margin, asp$figure.margin
      ),
      legend.position   = "none",
      axis.ticks.y      = ggplot2::element_blank(),
      axis.text.y       = ggplot2::element_blank(),
      axis.ticks.x      = ggplot2::element_line(
        linewidth = asp$figure.panel.line.size
      ),
      axis.text.x       = ggplot2::element_text(
        size = asp$figure.axis.text.size
      ),
      axis.title        = ggplot2::element_text(
        size = asp$figure.axis.title.size
      ),
      panel.border      = ggplot2::element_rect(
        fill = NA, linewidth = asp$figure.panel.line.size
      ),
      panel.grid.major  = ggplot2::element_blank(),
      panel.grid.minor  = ggplot2::element_blank()
    )
}


# =============================================================================
# Exported function: unmixed.mxn.plot
# =============================================================================

#' @title m \eqn{\times} n Unmixed Biplot Grid
#'
#' @description
#' Produces a grid of scatter biplots comparing a set of m target fluorophores
#' (x-axis) against every other channel in the dataset (y-axis), and saves the
#' result as a multi-panel PDF.  Each target fluorophore defines one row of
#' panels; the columns correspond to all other selected channels.
#'
#' Data are downsampled once to `max.points` events before any panel is drawn,
#' so all panels share the same cell population.  Axis minima are auto-scaled
#' per channel using the 1st percentile of the data, floored at the
#' user-supplied `x.min` / `y.min`.
#'
#' Panels can use either `geom_hex` (default, fastest) or the standard
#' `geom_scattermore` + `stat_density_2d` combination from [create.biplot()].
#'
#' @param unmixed.data Either a character string giving the path to an unmixed
#' FCS file, or a numeric matrix / data.frame whose columns are fluorophore
#' channels and whose rows are cells.
#' @param fluorophore Character vector of length \eqn{m}. Names of the target
#' fluorophores to place on the x-axis (one row of panels per entry). Every
#' element must be present in `colnames(unmixed.data)` (or the FCS columns).
#' @param asp The AutoSpectral parameter list, created by
#' [get.autospectral.param()].
#' @param channels Optional character vector. Subset of channels to use as
#' y-axis columns. When `NULL` (default) all columns except those in
#' `fluorophore` are used. Must be present in the data if supplied.
#' @param max.points Integer. Total number of events to retain after a single
#' random downsample applied before drawing any panel. Default `5e4`.
#' @param title Character string used as the PDF filename stem and as the
#' overall plot title. Default `"mxn_biplot"`.
#' @param biplot.size Numeric. Edge length in inches of each individual panel.
#' Default `3`.
#' @param n.col Integer. Number of panel columns per page. When `NULL`
#' (default) this is set to the number of y-axis channels, giving one row of
#' panels per target fluorophore on a single page.
#' @param x.min Numeric. Floor for the auto-scaled x-axis minimum (data
#' units). The actual minimum used is `min(auto, x.min)`. Default `-1000`.
#' @param y.min Numeric. Floor for the auto-scaled y-axis minimum (data
#' units). Default `-1000`.
#' @param x.width.basis Numeric. Width basis for the biexponential x-axis
#' transform (passed to [flowWorkspace::flowjo_biexp()]). Default `-1000`.
#' @param y.width.basis Numeric. Width basis for the biexponential y-axis
#' transform. Default `-1000`.
#' @param use.hex Logical. When `TRUE` (default), panels are rendered with
#' [ggplot2::geom_hex()] for speed. When `FALSE`, reverts to
#' [scattermore::geom_scattermore()] + [ggplot2::stat_density_2d()].
#' @param hex.bins Integer. Number of bins for `geom_hex` in each dimension.
#' Ignored when `use.hex = FALSE`. Default `64`.
#' @param color.palette Character string. Viridis palette name (`"viridis"`,
#' `"magma"`, `"inferno"`, `"plasma"`, `"cividis"`, `"rocket"`, `"mako"`,
#' `"turbo"`) or `"rainbow"` to use the package default gradient. Default
#' `"viridis"`.
#' @param output.dir Character string. Directory for the output PDF. Created
#' automatically if absent. Default `"."`.
#'
#' @return The combined [cowplot::plot_grid()] object is returned invisibly.
#' The PDF is always written to `output.dir`.
#'
#' @importFrom ggplot2 ggplot aes geom_hex geom_histogram scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous scale_fill_viridis_c scale_fill_gradientn
#' @importFrom ggplot2 stat_density_2d after_stat theme_bw theme element_line
#' @importFrom ggplot2 element_text element_rect element_blank margin labs ggsave
#' @importFrom scattermore geom_scattermore
#' @importFrom flowWorkspace flowjo_biexp
#' @importFrom cowplot plot_grid
#'
#' @seealso [unmixed.nxn.plot()], [create.biplot()]
#'
#' @export

unmixed.mxn.plot <- function(
    unmixed.data,
    fluorophore,
    asp,
    channels      = NULL,
    max.points    = 5e4,
    title         = "mxn_biplot",
    biplot.size   = 3,
    n.col         = NULL,
    x.min         = -1000,
    y.min         = -1000,
    x.width.basis = -1000,
    y.width.basis = -1000,
    use.hex       = TRUE,
    hex.bins      = 64,
    color.palette = "viridis",
    output.dir    = "."
  ) {

  # --- read and validate data ---
  mat <- .read.unmixed.input( unmixed.data, channels = NULL )

  missing.fluor <- setdiff( fluorophore, colnames( mat ) )
  if ( length( missing.fluor ) > 0 )
    stop(
      paste0(
        "The following `fluorophore` entries are absent from the data: ",
        paste( missing.fluor, collapse = ", " )
      ),
      call. = FALSE
    )

  # build y-axis channel set (everything not in `fluorophore`, then intersect
  # with user-supplied `channels` if provided)
  if ( is.null( channels ) ) {
    y.channels <- setdiff( colnames( mat ), fluorophore )
  } else {
    missing.ch <- setdiff( channels, colnames( mat ) )
    if ( length( missing.ch ) > 0 )
      stop(
        paste0(
          "The following `channels` are absent from the data: ",
          paste( missing.ch, collapse = ", " )
        ),
        call. = FALSE
      )
    y.channels <- setdiff( channels, fluorophore )
  }

  if ( length( y.channels ) == 0 )
    stop( "No y-axis channels remain after removing target `fluorophore`.",
          call. = FALSE )

  all.channels.needed <- unique( c( fluorophore, y.channels ) )
  mat <- mat[ , all.channels.needed, drop = FALSE ]

  # --- shared downsample ---
  if ( nrow( mat ) > max.points ) {
    set.seed( 42 )
    mat <- mat[ sample( nrow( mat ), max.points ), , drop = FALSE ]
  }
  message( "\033[34mPlotting ", nrow( mat ), " events across ",
           length( fluorophore ), " \u00d7 ", length( y.channels ),
           " panels\033[0m" )

  # --- auto-scale axis minima (floor at user-supplied value) ---
  x.mins <- vapply( fluorophore, function( ch ) {
    min( stats::quantile( mat[ , ch ], 0.01 ) * 2, x.min )
  }, numeric( 1 ) )

  y.mins <- vapply( y.channels, function( ch ) {
    min( stats::quantile( mat[ , ch ], 0.01 ) * 2, y.min )
  }, numeric( 1 ) )
  names( y.mins ) <- y.channels

  data.max <- asp$expr.data.max

  # --- layout ---
  n.x <- length( fluorophore )
  n.y <- length( y.channels )
  if ( is.null( n.col ) ) n.col <- n.y
  n.col <- min( n.col, n.y )

  # --- build all panels ---
  panel.list <- vector( "list", n.x * n.y )
  idx <- 1L

  for ( fl in fluorophore ) {

    tx <- .make.biexp.transforms(
      asp,
      x.min = x.mins[ fl ], x.max = data.max,
      y.min = y.min,         y.max = data.max,   # y transform recomputed per channel below
      x.width.basis = x.width.basis,
      y.width.basis = y.width.basis
    )$x   # only need the x transform here; y is per-channel

    x.axis <- .axis.setup( asp, x.mins[ fl ], data.max )

    for ( ch in y.channels ) {

      transforms <- .make.biexp.transforms(
        asp,
        x.min = x.mins[ fl ], x.max = data.max,
        y.min = y.mins[ ch ], y.max = data.max,
        x.width.basis = x.width.basis,
        y.width.basis = y.width.basis
      )
      y.axis <- .axis.setup( asp, y.mins[ ch ], data.max )

      panel.list[[ idx ]] <- .make.panel(
        plot.data     = mat,
        x.dim         = fl,
        y.dim         = ch,
        asp           = asp,
        tx            = transforms$x,
        ty            = transforms$y,
        x.axis        = x.axis,
        y.axis        = y.axis,
        use.hex       = use.hex,
        hex.bins      = hex.bins,
        color.palette = color.palette
      )
      idx <- idx + 1L
    }
  }

  # --- arrange and save ---
  combined <- cowplot::plot_grid(
    plotlist = panel.list,
    ncol     = n.col
  )

  if ( !dir.exists( output.dir ) ) dir.create( output.dir, recursive = TRUE )

  page.width  <- n.col * biplot.size
  page.height <- ceiling( ( n.x * n.y ) / n.col ) * biplot.size

  ggplot2::ggsave(
    filename  = file.path( output.dir, paste0( title, ".pdf" ) ),
    plot      = combined,
    device    = grDevices::cairo_pdf,
    width     = page.width,
    height    = page.height,
    limitsize = FALSE
  )

  message( "\033[32mSaved: ", file.path( output.dir, paste0( title, ".pdf" ) ),
           "\033[0m" )

  invisible( combined )
}
