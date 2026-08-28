# get_scatter_occupancy.r

#' @title Get Scatter Occupancy
#'
#' @description
#' Determines the region of 2D scatter space (typically FSC vs SSC) that a
#' sample's events actually occupy, using a density-relative threshold rather
#' than a fixed instrument range. A 2D kernel density estimate is built over a
#' robust, quantile-trimmed span of the data, so that a handful of extreme
#' events cannot inflate the estimation grid and coarsen its resolution over
#' the region that actually matters. Events falling in density bins below
#' `density.threshold` of the peak density are treated as sparse noise or
#' outlier events and excluded; the remaining "occupied" events define a
#' data-driven working range for gating and plotting.
#'
#' This identifies dense clusters generally, and does not by itself
#' distinguish a real debris population from the population of interest -- a
#' distinct, well-populated debris cloud is typically just as dense as real
#' cells or beads and will not be excluded here. Separating debris from the
#' target population remains the job of the downstream target-maximum search
#' (e.g. the bottom-left exclusion logic in `do.gate()`); this function only
#' keeps sparse, spread-out outlier events from distorting the working scale
#' used to set that search up.
#'
#' @importFrom MASS kde2d
#'
#' @param scatter.data A matrix or data frame with (at least) two columns:
#' Forward Scatter and Side Scatter, in that order.
#' @param density.threshold Numeric 0-1, default `0.02`. Bins with a kernel
#' density below this fraction of the peak density are excluded as sparse.
#' @param grid.n Numeric, default `128`. Binning grid for the kernel density
#' estimation.
#' @param max.events Numeric, default `1e5`. Maximum number of events used to
#' fit the kernel density surface. Above this, a random subsample is used to
#' fit the surface; every event is still evaluated against it.
#' @param trim.quantile Numeric 0-0.5, default `0.001`. Quantile trimmed from
#' each end of the data range before setting the density grid limits, so a
#' small number of extreme events cannot stretch the grid.
#' @param bird.seed Integer, seed for reproducible subsampling. Default `NULL`
#' skips seeding.
#'
#' @return A list with:
#' - `keep`: Logical vector, length `nrow(scatter.data)`, `TRUE` for events in
#' an occupied (dense) region.
#' - `x.range`, `y.range`: Two-element numeric vectors giving the occupied
#' range on each axis, taken from the retained events.
#' - `density`: The `x`/`y`/`z` kernel density surface used to determine
#' occupancy (same structure as `MASS::kde2d()` output), useful for
#' diagnostic plotting.
#'
#' @seealso
#' * [do.gate()]
#' * [define.gate.density()]
#' * [define.gate.landmarks()]
#'
#' @export

get.scatter.occupancy <- function(
    scatter.data,
    density.threshold = 0.02,
    grid.n = 128,
    max.events = 1e5,
    trim.quantile = 0.001,
    bird.seed = NULL
) {

  sc <- as.matrix( scatter.data[ , 1:2, drop = FALSE ] )
  n  <- nrow( sc )

  if ( n < 2 ) {
    return( list(
      keep = rep( TRUE, n ),
      x.range = if ( n == 1 ) rep( sc[ 1, 1 ], 2 ) else c( 0, 0 ),
      y.range = if ( n == 1 ) rep( sc[ 1, 2 ], 2 ) else c( 0, 0 ),
      density = NULL
    ) )
  }

  # robust grid limits: trim a small quantile from each end so a handful of
  # extreme events cannot stretch the density grid and coarsen its resolution
  # over the region that actually matters
  x.lim <- stats::quantile( sc[ , 1 ], c( trim.quantile, 1 - trim.quantile ), names = FALSE )
  y.lim <- stats::quantile( sc[ , 2 ], c( trim.quantile, 1 - trim.quantile ), names = FALSE )

  fit.idx <- seq_len( n )
  if ( length( fit.idx ) > max.events ) {
    if ( !is.null( bird.seed ) ) set.seed( bird.seed )
    fit.idx <- sample( fit.idx, max.events )
  }

  bw <- .safe.bandwidth( sc[ fit.idx, , drop = FALSE ] )

  if ( requireNamespace( "AutoSpectralRcpp", quietly = TRUE ) &&
       "fast_kde2d_cpp" %in% ls( getNamespace( "AutoSpectralRcpp" ) ) &&
       length( fit.idx ) > 1e4 ) {
    kde <- AutoSpectralRcpp::fast_kde2d_cpp(
      sc[ fit.idx, 1 ], sc[ fit.idx, 2 ],
      n = grid.n,
      h = bw * 0.1,
      x_limits = x.lim,
      y_limits = y.lim
    )
  } else {
    kde <- MASS::kde2d(
      sc[ fit.idx, 1 ], sc[ fit.idx, 2 ],
      h = bw, n = grid.n,
      lims = c( x.lim, y.lim )
    )
  }

  ix <- findInterval( sc[ , 1 ], kde$x, all.inside = TRUE )
  iy <- findInterval( sc[ , 2 ], kde$y, all.inside = TRUE )
  keep <- kde$z[ cbind( ix, iy ) ] >= density.threshold * max( kde$z )

  if ( !any( keep ) ) {
    # pathological case (e.g. threshold too strict for very sparse data);
    # fall back to keeping everything rather than returning an empty gate
    keep <- rep( TRUE, n )
  }

  x.range <- range( sc[ keep, 1 ] )
  y.range <- range( sc[ keep, 2 ] )
  if ( diff( x.range ) == 0 )
    x.range <- x.range + c( -1, 1 ) * max( abs( x.range[ 1 ] ), 1 ) * 0.01
  if ( diff( y.range ) == 0 )
    y.range <- y.range + c( -1, 1 ) * max( abs( y.range[ 1 ] ), 1 ) * 0.01

  list(
    keep = keep,
    x.range = x.range,
    y.range = y.range,
    density = kde
  )
}
