# nnls_normal.R

#' @title Non-Negative Least Squares from Normal Equations
#'
#' @description
#' Solves \eqn{\min_{x \ge 0} \| A x - b \|^2} given only the normal
#' equations \eqn{G = A^\top A} and \eqn{h = A^\top b}, by the Lawson-Hanson
#' active-set method. Working from `G` and `h` keeps the cost in the number
#' of columns, so a problem with many rows can be accumulated block by block
#' without ever holding `A`.
#'
#' Columns are scaled to unit norm (`sqrt(diag(G))`) before solving, which
#' leaves the solution unchanged and evens out the conditioning. A column
#' with zero norm stays at zero. Each passive-set solve uses a Cholesky
#' factorisation, falling back to a pivoted QR solve when the passive block
#' is not numerically positive definite.
#'
#' @param G Numeric matrix, \eqn{A^\top A} (p x p, symmetric).
#' @param h Numeric vector, \eqn{A^\top b} (length p).
#' @param max.iter Integer, maximum passive-set changes. Default `3 * p`.
#' @param tol Numeric, a column enters the passive set only when its scaled
#'   gradient exceeds this. Default `10 * eps * p * max(1, max(|h_s|))`,
#'   with `h_s` the scaled `h`.
#'
#' @return Numeric vector, the non-negative solution (length p).
#'
#' @noRd

.nnls.normal <- function( G, h, max.iter = NULL, tol = NULL ) {

  p <- length( h )
  if ( p == 0L ) return( numeric( 0 ) )

  scale <- sqrt( pmax( diag( G ), 0 ) )
  scale[ !is.finite( scale ) | scale == 0 ] <- 1
  Gs <- G / tcrossprod( scale )
  hs <- h / scale

  if ( is.null( max.iter ) ) max.iter <- 3L * p
  if ( is.null( tol ) ) tol <- 10 * .Machine$double.eps * p * max( 1, max( abs( hs ) ) )

  spd.solve <- function( idx ) {
    a <- Gs[ idx, idx, drop = FALSE ]
    r <- tryCatch( chol( a ), error = function( e ) NULL )
    if ( !is.null( r ) )
      return( backsolve( r, forwardsolve( t( r ), hs[ idx ] ) ) )
    qr.coef( qr( a, LAPACK = TRUE ), hs[ idx ] )
  }

  x       <- numeric( p )
  passive <- logical( p )
  blocked <- logical( p )
  w       <- hs
  iter    <- 0L

  while ( iter < max.iter ) {

    candidate <- !passive & !blocked & w > tol
    if ( !any( candidate ) ) break

    t.idx <- which( candidate )[ which.max( w[ candidate ] ) ]
    passive[ t.idx ] <- TRUE

    z <- numeric( p )
    z[ passive ] <- spd.solve( which( passive ) )

    # the entering column cannot be positive (rounding): leave it out until x moves
    if ( z[ t.idx ] <= 0 ) {
      passive[ t.idx ] <- FALSE
      blocked[ t.idx ] <- TRUE
      next
    }

    while ( any( z[ passive ] <= 0 ) && iter < max.iter ) {
      iter  <- iter + 1L
      neg   <- which( passive & z <= 0 )
      ratio <- x[ neg ] / ( x[ neg ] - z[ neg ] )
      alpha <- min( ratio )
      x     <- x + alpha * ( z - x )
      x[ neg[ ratio <= alpha ] ] <- 0
      passive[ passive & x <= 0 ] <- FALSE
      x[ !passive ] <- 0
      z <- numeric( p )
      if ( any( passive ) ) z[ passive ] <- spd.solve( which( passive ) )
    }

    iter    <- iter + 1L
    x       <- z
    blocked <- logical( p )
    w       <- hs - as.vector( Gs %*% x )
  }

  x / scale
}
