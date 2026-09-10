# compare_unmixed_data.R

#' @title Compare Two or More Sets of Already-Unmixed Data
#'
#' @description
#' Lays out biplots for two or more pre-unmixed data sets (matrices or data
#' frames of the same panel) as a grid, one row per data set, so that the
#' same channel pairs can be compared side-by-side across conditions (for
#' example, cells vs. beads, or before vs. after a correction step).
#'
#' Channels are paired up sequentially rather than plotted N x N: for
#' `channels = c("APC", "FITC", "PE", "PerCP", "BV421")` the pairs are
#' APC x FITC, PE x PerCP, and BV421 x APC (the first channel is repeated
#' to complete the last pair when an odd number of channels is requested).
#' Each row uses [create.biplot()] for every pair, and the panels are
#' stitched together with `cowplot::plot_grid()`.
#'
#' If `gate.boundary` is supplied and the `scatter.param` columns are present
#' in every data set, [apply.gate()] is used to gate every data set to that
#' boundary before plotting. If any data set is missing a scatter column,
#' the gate is not applied to any of them (a warning names which data
#' set(s) were missing it), so that every row in the comparison reflects the
#' same population.
#'
#' If `data.list` is a named list, the names are used as row labels on the
#' right of the grid. If unnamed, and `data.list` was constructed inline as
#' `list(...)` at the call site, the deparsed expression for each element
#' (i.e. the variable name passed in) is used instead.
#'
#' @param data.list A list of two or more matrices or data frames of
#' unmixed flow cytometry data, one per data set to compare. Every element
#' must have named columns. Column names do not need to match exactly
#' across data sets -- see `channels` below.
#' @param channels Character vector of channel names to plot, in the order
#' they should be paired (see Description). Must resolve to at least 2
#' channels present in every element of `data.list` after the intersect
#' described below.
#' @param asp The AutoSpectral parameter list, prepared using
#' `get.autospectral.param()`. Passed through to [create.biplot()] for every
#' panel, and supplies the default for `scatter.param`.
#' @param gate.boundary Optional gate boundary, as returned by
#' `define.gate.landmarks()`, `define.gate.density()`, or `do.gate()` -- a
#' list containing at least numeric `x` and `y` components describing the
#' polygon vertices. When supplied, every data set is gated with
#' [apply.gate()] before plotting, provided all data sets carry
#' `scatter.param`. Default `NULL` (no gating).
#' @param scatter.param Character vector of length 2 giving the names of the
#' two scatter columns to gate on, passed to [apply.gate()]. Only used when
#' `gate.boundary` is supplied. Default `asp$default.scatter.parameter`.
#' @param variants The variant list returned by `get.spectral.variants()`,
#' passed through to every [create.biplot()] call for reference curves.
#' Default `NULL` (no reference curves).
#' @param spread.kappa Numeric, passed through to [create.biplot()]. Default `2`.
#' @param x.min,y.min Numeric, axis minima (data units) passed through to
#' [create.biplot()]. Default `-5000` for both.
#' @param x.max,y.max Numeric, axis maxima (data units) passed through to
#' [create.biplot()]. Default `asp$expr.data.max` for both.
#' @param x.width.basis,y.width.basis Numeric, biexponential width bases
#' passed through to [create.biplot()]. Default `-1000` for both.
#' @param max.points Numeric, per-panel point cap passed through to
#' [create.biplot()]. Default `5e5`.
#' @param color.palette Character, passed through to [create.biplot()].
#' Default `"rainbow"`.
#' @param panel.width,row.height Numeric, width and height (inches) of a
#' single biplot panel, used to size the saved figure. Defaults `3` and `3`.
#' @param label.width Numeric, width (inches) reserved for the row-label
#' column on the right of the figure. Default `1`.
#' @param label.size Numeric, font size (points) for the row labels.
#' Default `12`.
#' @param label.angle Numeric, rotation angle (degrees) for the row labels.
#' Default `-90` (reads bottom-to-top, matching a facet_grid right strip).
#' @param save Logical, if `TRUE` (default), saves a JPEG file to
#' `output.dir`. Otherwise the combined plot is only printed.
#' @param title Character, used as the JPEG filename stem. Default
#' `"unmix_comparison"`.
#' @param output.dir Optional output directory. Default `NULL`, in which
#' case the current working directory is used.
#'
#' @return Invisibly, the combined `cowplot` object.
#'
#' @seealso
#' * [create.biplot()]
#' * [apply.gate()]
#' * [compare.unmix()]
#'
#' @importFrom cowplot plot_grid ggdraw draw_label
#' @importFrom ggplot2 ggsave
#' @importFrom ragg agg_jpeg
#'
#' @export

compare.unmixed.data <- function(
    data.list,
    channels,
    asp,
    gate.boundary   = NULL,
    scatter.param   = asp$default.scatter.parameter,
    variants        = NULL,
    spread.kappa    = 2,
    x.min           = -5000,
    x.max           = asp$expr.data.max,
    y.min           = -5000,
    y.max           = asp$expr.data.max,
    x.width.basis   = -1000,
    y.width.basis   = -1000,
    max.points      = 5e5,
    color.palette   = "rainbow",
    panel.width     = 3,
    row.height      = 3,
    label.width     = 1,
    label.size      = 12,
    label.angle     = -90,
    save            = TRUE,
    title           = "unmix_comparison",
    output.dir      = NULL
) {

  # capture the call expression for data.list before any modification, so we
  # can fall back to it for row labels if data.list is unnamed
  data.list.expr <- substitute( data.list )

  # --- basic input checks ---
  if ( !is.list( data.list ) || length( data.list ) < 2 ) {
    stop( "data.list must be a list of two or more matrices or data frames.", call. = FALSE )
  }

  for ( i in seq_along( data.list ) ) {
    if ( is.null( dim( data.list[[ i ]] ) ) || is.null( colnames( data.list[[ i ]] ) ) ) {
      stop(
        paste0( "Element ", i, " of data.list must be a matrix or data frame with named columns." ),
        call. = FALSE
      )
    }
  }

  # --- resolve row labels ---
  supplied.names <- names( data.list )
  if ( is.null( supplied.names ) ) supplied.names <- rep( "", length( data.list ) )
  missing.name.idx <- which( is.na( supplied.names ) | supplied.names == "" )

  if ( length( missing.name.idx ) > 0 ) {

    fallback.names <- rep( NA_character_, length( data.list ) )

    if ( is.call( data.list.expr ) && identical( as.character( data.list.expr[[ 1 ]] ), "list" ) ) {
      arg.exprs <- as.list( data.list.expr )[ -1 ]
      if ( length( arg.exprs ) == length( data.list ) ) {
        fallback.names <- vapply( arg.exprs, function( e ) deparse( e ), character( 1 ) )
      }
    } else if ( is.symbol( data.list.expr ) ) {
      fallback.names <- paste0( deparse( data.list.expr ), seq_along( data.list ) )
    }

    for ( i in missing.name.idx ) {
      supplied.names[ i ] <- if ( !is.na( fallback.names[ i ] ) ) fallback.names[ i ] else paste0( "dataset", i )
    }
  }

  names( data.list ) <- supplied.names

  # --- resolve channels present in every data set ---
  common.cols <- Reduce( intersect, lapply( data.list, colnames ) )
  missing.channels <- setdiff( channels, common.cols )

  if ( length( missing.channels ) > 0 ) {
    warning(
      paste0(
        "The following requested channel(s) are not present in every data set and have been dropped: ",
        paste( missing.channels, collapse = ", " )
      ),
      call. = FALSE
    )
  }

  channels <- intersect( channels, common.cols )

  if ( length( channels ) < 2 ) {
    stop( "At least 2 channels present in every data set are required to plot.", call. = FALSE )
  }

  # repeat the first channel to complete the last pair on an odd count
  if ( length( channels ) %% 2 == 1 ) {
    channels <- c( channels, channels[ 1 ] )
  }

  pair.idx <- matrix( seq_along( channels ), ncol = 2, byrow = TRUE )
  n.pairs  <- nrow( pair.idx )

  # --- optional gating, all-or-nothing across data sets ---
  if ( !is.null( gate.boundary ) ) {

    if ( is.null( scatter.param ) || length( scatter.param ) != 2 ) {
      stop(
        "scatter.param must be a length-2 character vector (e.g. asp$default.scatter.parameter) to apply gate.boundary.",
        call. = FALSE
      )
    }

    scatter.present <- vapply(
      data.list, function( d ) all( scatter.param %in% colnames( d ) ), logical( 1 )
    )

    if ( all( scatter.present ) ) {
      data.list <- lapply(
        data.list,
        function( d ) apply.gate( d, gate.boundary, scatter.param = scatter.param, asp = asp, min.fraction = 0 )
      )
    } else {
      warning(
        paste0(
          "gate.boundary was supplied but scatter.param (", paste( scatter.param, collapse = ", " ),
          ") is missing from: ", paste( names( data.list )[ !scatter.present ], collapse = ", " ),
          ". Gate not applied to any data set."
        ),
        call. = FALSE
      )
    }
  }

  # --- build one biplot per data set x channel-pair, row-major ---
  plot.list <- list()

  for ( i in seq_along( data.list ) ) {
    for ( p in seq_len( n.pairs ) ) {

      x.dim <- channels[ pair.idx[ p, 1 ] ]
      y.dim <- channels[ pair.idx[ p, 2 ] ]

      plot.list[[ length( plot.list ) + 1 ]] <- create.biplot(
        data.list[[ i ]],
        x.dim         = x.dim,
        y.dim         = y.dim,
        asp           = asp,
        variants      = variants,
        spread.kappa  = spread.kappa,
        x.min         = x.min,
        x.max         = x.max,
        y.min         = y.min,
        y.max         = y.max,
        x.width.basis = x.width.basis,
        y.width.basis = y.width.basis,
        max.points    = max.points,
        color.palette = color.palette,
        save          = FALSE
      )
    }
  }

  main.grid <- cowplot::plot_grid(
    plotlist = plot.list,
    nrow     = length( data.list ),
    ncol     = n.pairs,
    align    = "hv",
    axis     = "tblr"
  )

  # --- row-label column on the right ---
  label.plots <- lapply( names( data.list ), function( n ) {
    cowplot::ggdraw() +
      cowplot::draw_label( n, size = label.size, angle = label.angle )
  } )

  label.col <- cowplot::plot_grid( plotlist = label.plots, nrow = length( data.list ), ncol = 1 )

  grid.width <- n.pairs * panel.width
  combined.plot <- cowplot::plot_grid(
    main.grid, label.col,
    ncol       = 2,
    rel_widths = c( grid.width, label.width )
  )

  # --- save or print ---
  if ( save ) {

    if ( is.null( output.dir ) ) output.dir <- getwd()
    if ( !dir.exists( output.dir ) ) dir.create( output.dir, recursive = TRUE )

    ggsave(
      file.path( output.dir, sprintf( "%s.jpg", title ) ),
      plot      = combined.plot,
      device    = ragg::agg_jpeg,
      width     = grid.width + label.width,
      height    = length( data.list ) * row.height,
      limitsize = FALSE
    )
  }

  print( combined.plot )

  invisible( combined.plot )
}
