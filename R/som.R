# som.R

#' @title Cluster Data with a Self-Organizing Map
#'
#' @description
#' Replaces `FlowSOM::SOM()` / `EmbedSOM::SOM()` as the SOM clustering engine
#' behind `get.af.spectra()` and `get.fluor.variants()`. Uses the OpenMP-
#' accelerated batch SOM (`som_train_batch_cpp()`) from AutoSpectralRcpp when
#' available, falling back to `FlowSOM::SOM()` (pure R, no compiled
#' dependency) otherwise -- so a plain `AutoSpectral` install with no
#' `AutoSpectralRcpp` still works, just without the OpenMP speedup.
#'
#' @importFrom FlowSOM SOM
#'
#' @param data Numeric matrix, training events x features. Must have colnames.
#' @param som.dim Integer, side length of the square SOM grid.
#' @param rlen Integer, number of epochs (full passes over the data).
#'   Default `10`.
#' @param radius Length-2 numeric vector, start/end neighbourhood radius.
#'   Default `NULL`, which derives both ends from the grid's own neighbour
#'   distances (67th percentile down to a small non-zero floor). Only used
#'   on the AutoSpectralRcpp path -- `FlowSOM::SOM()` manages its own radius
#'   schedule internally.
#' @param dist Integer 1:4, distance function (1 manhattan, 2 euclidean,
#'   3 chebyshev, 4 cosine). Default `4`.
#' @param seed Integer, RNG seed for the initial codebook sample. Callers
#'   should pass `asp$bird.seed`.
#' @param threads Integer, OpenMP threads for the accelerated path. Default
#'   `0` (all available cores). Ignored on the `FlowSOM::SOM()` fallback
#'   path, which is single-threaded.
#' @param unit.norm Logical, default `FALSE`. Scale every training row to unit
#'   L2 length before training. Appropriate with `dist = 4` (cosine), where the
#'   distortion-minimising code is the mean of unit event vectors rather than
#'   the mean of raw ones; without it, bright events dominate each node's code.
#'   Leaves cosine distances, and therefore node assignments, unchanged. Has no
#'   principled justification under `dist = 2` and is left off by default so
#'   existing callers are unaffected. Note that the returned codes are then on
#'   the unit-normalised scale.
#'
#' @return A list with `codes` (matrix, SOM nodes x features), `grid`, and
#'   `nNodes`, matching the subset of `FlowSOM::SOM()`'s return value
#'   actually consumed elsewhere in AutoSpectral (`map$codes`).
#'
#' @keywords internal
#' @export

get.som.codes <- function(
    data,
    som.dim,
    rlen      = 10L,
    radius    = NULL,
    dist      = 4L,
    seed      = 1337L,
    threads   = 0L,
    unit.norm = FALSE
) {

  if ( is.null( colnames( data ) ) )
    stop( "`data` must have colnames.", call. = FALSE )

  # Under cosine distance the code minimising within-node distortion is the
  # direction of the sum of unit event vectors, not the arithmetic mean of the
  # raw vectors. The batch trainer accumulates raw values, so scaling each row
  # to unit length here makes the update consistent with the distance. Row
  # scaling leaves cosine distances unchanged, so assignments are unaffected;
  # only the relative influence of bright and dim events on each code changes.
  if ( unit.norm ) {
    row.norm <- sqrt( rowSums( data^2 ) )
    keep     <- row.norm > 0
    if ( !all( keep ) ) {
      data     <- data[ keep, , drop = FALSE ]
      row.norm <- row.norm[ keep ]
    }
    data <- data / row.norm
  }

  grid   <- expand.grid( seq_len( som.dim ), seq_len( som.dim ) )
  ncodes <- nrow( grid )

  if ( nrow( data ) < ncodes )
    stop(
      paste0(
        "Not enough events (", nrow( data ), ") to initialise a ",
        som.dim, "x", som.dim, " SOM (", ncodes, " nodes)."
      ),
      call. = FALSE
    )

  nhbrdist <- as.matrix( stats::dist( grid, method = "maximum" ) )

  if ( is.null( radius ) )
    radius <- as.numeric( stats::quantile( nhbrdist, 0.67 ) * c( 1, 0.1 ) )
  # NOTE: batch SOM's neighbourhood kernel is exp(-d^2/radius^2), so unlike
  # the retired online trainer the end radius must stay strictly positive --
  # collapsing to 0 starves every node outside the exact best-matching unit.

  has.rcpp.som <- requireNamespace( "AutoSpectralRcpp", quietly = TRUE ) &&
    ( "som_train_batch_cpp" %in% ls( getNamespace( "AutoSpectralRcpp" ) ) )

  set.seed( seed )

  if ( has.rcpp.som ) {

    init.codes <- data[ sample( seq_len( nrow( data ) ), ncodes, replace = FALSE ), , drop = FALSE ]
    radii      <- seq( radius[ 1 ], radius[ 2 ], length.out = rlen )

    codes <- AutoSpectralRcpp::som_train_batch_cpp(
      data       = data,
      init_codes = init.codes,
      nhbrdist   = nhbrdist,
      radii      = radii,
      dist       = dist,
      n_threads  = threads
    )

    colnames( codes ) <- colnames( data )

  } else {

    if ( !requireNamespace( "FlowSOM", quietly = TRUE ) )
      stop(
        "get.som.codes() requires either AutoSpectralRcpp (fast, ",
        "recommended) or the FlowSOM package, and neither is installed.",
        call. = FALSE
      )

    map   <- FlowSOM::SOM(
      data,
      xdim = som.dim,
      ydim = som.dim,
      rlen = rlen,
      distf = dist,
      silent = TRUE
    )
    codes <- map$codes

  }

  list( codes = codes, grid = grid, nNodes = ncodes )
}
