# select_noise_events.R

#' @title Select Events for the Noise Model
#'
#' @description
#' The scatter core of a sample and its singlets, independent of any gate,
#' for `estimate.noise.model()`:
#'
#' 1. The density mode of log10 FSC-A x log10 SSC-A: an 80 x 80 histogram
#'    over the 0.5-99.5th percentiles, 5 x 5 box-smoothed.
#' 2. Concentration steps from there: 15 steps of h = 50%, each centred on
#'    the median of the inner half, with the covariance scaled by the
#'    consistency factor 3.26 and 0.01^2 added to its diagonal, so that a
#'    saturated axis does not collapse it.
#' 3. The gate: Mahalanobis^2 below 5.991 (95%).
#' 4. Singlets, when FSC-H is given: |log10(A / H) - median| < 4 MAD.
#'
#' Saturated events are kept (bright beads saturate SSC, or FSC and SSC, on
#' some cytometers); events with non-positive or non-finite FSC-A or SSC-A
#' are not. Doublets and off-core particles add variance that grows with
#' signal, which the noise fit would otherwise read as shot noise.
#'
#' @param fsc.a,ssc.a Numeric vectors, forward and side scatter areas.
#' @param fsc.h Optional numeric vector, forward scatter height. Default
#'   `NULL` (no singlet band).
#' @param cutoff Numeric, Mahalanobis^2 cut. Default `5.991`.
#'
#' @return A list: `mask` (logical, one per event) and `counts`
#'   (`n.root`, `n.scatter`, `n.singlet`).
#'
#' @noRd

.select.noise.events <- function( fsc.a, ssc.a, fsc.h = NULL, cutoff = 5.991 ) {
  core <- .scatter.core.gate( fsc.a, ssc.a, cutoff )
  mask <- if ( is.null( fsc.h ) ) core else .singlet.band( fsc.a, fsc.h, core )
  list( mask = mask, counts = c( n.root = length( mask ), n.scatter = sum( core ), n.singlet = sum( mask ) ) )
}

.scatter.core.gate <- function( fsc.a, ssc.a, cutoff = 5.991, n.steps = 15L,
                                consistency = 3.26, hist.bins = 80L, smooth = 5L ) {
  fsc.a <- as.numeric( fsc.a )
  ssc.a <- as.numeric( ssc.a )
  ok <- is.finite( fsc.a ) & is.finite( ssc.a ) & fsc.a > 0 & ssc.a > 0
  if ( sum( ok ) < 3L ) return( ok )
  xy <- cbind( log10( pmax( fsc.a, 1 ) ), log10( pmax( ssc.a, 1 ) ) )
  xy[ !is.finite( xy ) ] <- 0

  # density mode
  pts <- xy[ ok, , drop = FALSE ]
  lo  <- apply( pts, 2, stats::quantile, probs = 0.005, names = FALSE )
  hi  <- apply( pts, 2, stats::quantile, probs = 0.995, names = FALSE )
  upper <- hi + 1e-9
  width <- ( upper - lo ) / hist.bins
  inside <- pts[ , 1 ] >= lo[ 1 ] & pts[ , 1 ] <= upper[ 1 ] &
    pts[ , 2 ] >= lo[ 2 ] & pts[ , 2 ] <= upper[ 2 ]
  ix <- pmin( floor( ( pts[ inside, 1 ] - lo[ 1 ] ) / width[ 1 ] ), hist.bins - 1L ) + 1L
  iy <- pmin( floor( ( pts[ inside, 2 ] - lo[ 2 ] ) / width[ 2 ] ), hist.bins - 1L ) + 1L
  hist <- matrix( tabulate( ( iy - 1L ) * hist.bins + ix, nbins = hist.bins^2 ), hist.bins, hist.bins )
  pad <- smooth %/% 2L
  padded <- matrix( 0, hist.bins + 2L * pad, hist.bins + 2L * pad )
  padded[ pad + seq_len( hist.bins ), pad + seq_len( hist.bins ) ] <- hist
  dens <- matrix( 0, hist.bins, hist.bins )
  for ( i in seq_len( smooth ) ) for ( j in seq_len( smooth ) )
    dens <- dens + padded[ i - 1L + seq_len( hist.bins ), j - 1L + seq_len( hist.bins ) ]
  top <- which( dens == max( dens ), arr.ind = TRUE )
  top <- top[ order( top[ , 1 ], top[ , 2 ] ), , drop = FALSE ][ 1, ]
  centre <- c( lo[ 1 ] + ( top[ 1 ] - 0.5 ) * width[ 1 ], lo[ 2 ] + ( top[ 2 ] - 0.5 ) * width[ 2 ] )
  cov.m  <- diag( ( ( hi - lo ) / 20 )^2 )

  reg <- diag( c( 0.01^2, 0.01^2 ) )
  for ( step in seq_len( n.steps ) ) {
    md2 <- stats::mahalanobis( xy, centre, cov.m )
    md2[ !ok ] <- Inf
    inner <- md2 <= stats::quantile( md2[ ok ], 0.5, names = FALSE ) & md2 < 4 * cutoff
    if ( sum( inner ) < 3L ) break
    centre <- apply( xy[ inner, , drop = FALSE ], 2, stats::median )
    cov.m  <- stats::cov( xy[ inner, , drop = FALSE ] ) * consistency + reg
  }
  ok & stats::mahalanobis( xy, centre, cov.m ) < cutoff
}

.singlet.band <- function( fsc.a, fsc.h, core, n.mads = 4 ) {
  if ( !any( core ) ) return( core )
  ratio <- log10( pmax( as.numeric( fsc.a ), 1 ) ) - log10( pmax( as.numeric( fsc.h ), 1 ) )
  med <- stats::median( ratio[ core ] )
  mad <- 1.4826 * stats::median( abs( ratio[ core ] - med ) )
  core & abs( ratio - med ) < n.mads * max( mad, 1e-3 )
}
