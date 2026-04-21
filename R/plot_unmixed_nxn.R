# plot_unmixed_nxn.R

#' @title n \eqn{\times} n Unmixed Biplot Triangle
#'
#' @description
#' Produces a lower-triangular grid of scatter biplots for every pair of
#' channels in an unmixed dataset, with 1-D histograms on the diagonal, and
#' saves the result as a single large-format PDF.
#'
#' For \eqn{n} channels the grid is \eqn{n \times n} panels.  The upper
#' triangle and all cells above the diagonal are left blank, giving
#' \eqn{n(n-1)/2} scatter panels and \eqn{n} histogram panels.  Panel
#' \eqn{(i, j)} (row \eqn{i}, column \eqn{j}) shows column \eqn{j} on the
#' x-axis and row \eqn{i} on the y-axis, consistent with the convention used
#' by [unmixed.mxn.plot()]: the named fluorophore always appears on x.
#'
#' Data are downsampled once to `max.points` events before any panel is drawn
#' so all panels share an identical cell population.  Axis minima are
#' auto-scaled per channel from the 1st percentile of the data, floored at
#' the user-supplied `x.min` / `y.min`.
#'
#' @param unmixed.data Either a character string giving the path to an unmixed
#' FCS file, or a numeric matrix / data.frame whose columns are fluorophore
#' channels and whose rows are cells.
#' @param asp The AutoSpectral parameter list, created by
#' [get.autospectral.param()].
#' @param channels Optional character vector. Subset of columns to include.
#' When `NULL` (default) all columns are used. Large panels sets (\eqn{n > 20})
#' will produce very large PDFs; consider subsetting.
#' @param max.points Integer. Total number of events retained after a single
#' random downsample applied before any panel is drawn. Default `5e4`.
#' @param title Character string used as the PDF filename stem. Default
#' `"nxn_biplot"`.
#' @param biplot.size Numeric. Edge length in inches of each individual panel
#' (both scatter and histogram panels). The saved PDF will be
#' `n * biplot.size` \eqn{\times} `n * biplot.size` inches. Default `3`.
#' @param x.min Numeric. Floor for the auto-scaled x-axis minimum (data
#' units). The actual minimum used is `min(auto, x.min)`. Default `-1000`.
#' @param y.min Numeric. Floor for the auto-scaled y-axis minimum (data
#' units). Default `-1000`.
#' @param x.width.basis Numeric. Width basis for the biexponential x-axis
#' transform. Default `-1000`.
#' @param y.width.basis Numeric. Width basis for the biexponential y-axis
#' transform. Default `-1000`.
#' @param use.hex Logical. When `TRUE` (default), scatter panels are rendered
#' with [ggplot2::geom_hex()] for speed. When `FALSE`, reverts to
#' [scattermore::geom_scattermore()] + [ggplot2::stat_density_2d()].
#' @param hex.bins Integer. Number of bins for `geom_hex` in each dimension.
#' Ignored when `use.hex = FALSE`. Default `64`.
#' @param color.palette Character string. Viridis palette name (`"viridis"`,
#' `"magma"`, `"inferno"`, `"plasma"`, `"cividis"`, `"rocket"`, `"mako"`,
#' `"turbo"`) or `"rainbow"` to use the package default gradient. Applies to
#' scatter panels. Default `"viridis"`.
#' @param hist.fill Character string. Fill color for the diagonal histogram
#' bars. Default `"steelblue"`.
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
#' @seealso [unmixed.mxn.plot()], [create.biplot()]
#'
#' @export

unmixed.nxn.plot <- function(
    unmixed.data,
    asp,
    channels      = NULL,
    max.points    = 5e4,
    title         = "nxn_biplot",
    biplot.size   = 3,
    x.min         = -1000,
    y.min         = -1000,
    x.width.basis = -1000,
    y.width.basis = -1000,
    use.hex       = TRUE,
    hex.bins      = 64,
    color.palette = "viridis",
    hist.fill     = "steelblue",
    output.dir    = "."
  ) {

  # --- read and validate data ---
  mat <- .read.unmixed.input( unmixed.data, channels = channels )

  n <- ncol( mat )
  if ( n < 2 )
    stop( "At least 2 channels are required for an n\u00d7n plot.", call. = FALSE )

  ch.names <- colnames( mat )

  message( "\033[34mBuilding ", n, "\u00d7", n, " triangle: ",
           ( n * ( n - 1 ) ) / 2, " scatter panels + ", n,
           " histograms\033[0m" )

  # --- shared downsample ---
  if ( nrow( mat ) > max.points ) {
    set.seed( 42 )
    mat <- mat[ sample( nrow( mat ), max.points ), , drop = FALSE ]
  }

  # --- auto-scale axis minima per channel (floor at user-supplied value) ---
  ch.mins <- vapply( ch.names, function( ch ) {
    min( stats::quantile( mat[ , ch ], 0.01 ) * 2, x.min )
  }, numeric( 1 ) )
  # y.min floor is applied separately; since the grid is symmetric we use the
  # same per-channel minimum for both axes — whichever role the channel plays.
  ch.mins.y <- vapply( ch.names, function( ch ) {
    min( stats::quantile( mat[ , ch ], 0.01 ) * 2, y.min )
  }, numeric( 1 ) )

  data.max <- asp$expr.data.max

  # pre-compute per-channel axis setup objects (x role and y role may differ
  # only in the floor applied, so we keep both)
  x.axes <- lapply( ch.names, function( ch )
    .axis.setup( asp, ch.mins[ ch ], data.max )
  )
  names( x.axes ) <- ch.names

  y.axes <- lapply( ch.names, function( ch )
    .axis.setup( asp, ch.mins.y[ ch ], data.max )
  )
  names( y.axes ) <- ch.names

  # pre-compute per-channel biexp transforms (x and y roles)
  # We need a transform object per channel per role; because the only thing
  # that differs is x.min / y.min we store them by channel name.
  x.transforms <- lapply( ch.names, function( ch ) {
    .make.biexp.transforms(
      asp,
      x.min = ch.mins[ ch ],   x.max = data.max,
      y.min = 0,               y.max = data.max,   # y unused here
      x.width.basis = x.width.basis,
      y.width.basis = y.width.basis
    )$x
  } )
  names( x.transforms ) <- ch.names

  y.transforms <- lapply( ch.names, function( ch ) {
    .make.biexp.transforms(
      asp,
      x.min = 0,               x.max = data.max,   # x unused here
      y.min = ch.mins.y[ ch ], y.max = data.max,
      x.width.basis = x.width.basis,
      y.width.basis = y.width.basis
    )$y
  } )
  names( y.transforms ) <- ch.names

  # --- build the n×n cell list (row-major) ---
  # Row i, column j:
  #   i == j  → diagonal histogram for channel i
  #   i >  j  → lower triangle: x = channel j, y = channel i
  #   i <  j  → upper triangle: NULL (blank)
  total.cells <- n * n
  panel.list  <- vector( "list", total.cells )

  for ( i in seq_len( n ) ) {
    for ( j in seq_len( n ) ) {
      cell.idx <- ( i - 1L ) * n + j

      if ( i == j ) {
        # diagonal: 1-D histogram
        panel.list[[ cell.idx ]] <- .make.diagonal.hist(
          plot.data     = mat,
          dim           = ch.names[ i ],
          asp           = asp,
          tx            = x.transforms[[ ch.names[ i ] ]],
          x.axis        = x.axes[[ ch.names[ i ] ]],
          fill.color    = hist.fill
        )

      } else if ( i > j ) {
        # lower triangle: x = ch j, y = ch i
        panel.list[[ cell.idx ]] <- .make.panel(
          plot.data     = mat,
          x.dim         = ch.names[ j ],
          y.dim         = ch.names[ i ],
          asp           = asp,
          tx            = x.transforms[[ ch.names[ j ] ]],
          ty            = y.transforms[[ ch.names[ i ] ]],
          x.axis        = x.axes[[ ch.names[ j ] ]],
          y.axis        = y.axes[[ ch.names[ i ] ]],
          use.hex       = use.hex,
          hex.bins      = hex.bins,
          color.palette = color.palette
        )

      }
      # upper triangle: leave NULL (cowplot treats NULL as blank space)
    }
  }

  # --- assemble grid ---
  combined <- cowplot::plot_grid(
    plotlist = panel.list,
    ncol     = n,
    nrow     = n
  )

  # --- save ---
  if ( !dir.exists( output.dir ) ) dir.create( output.dir, recursive = TRUE )

  canvas.size <- n * biplot.size

  ggplot2::ggsave(
    filename  = file.path( output.dir, paste0( title, ".pdf" ) ),
    plot      = combined,
    device    = grDevices::cairo_pdf,
    width     = canvas.size,
    height    = canvas.size,
    limitsize = FALSE
  )

  message( "\033[32mSaved: ", file.path( output.dir, paste0( title, ".pdf" ) ),
           "\033[0m" )

  invisible( combined )
}
