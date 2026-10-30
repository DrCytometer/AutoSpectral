# fix_my_unmix_anchored.R

#' Half-sample mode of a numeric vector (Bickel's estimator).
#'
#' Repeatedly keeps the densest half of the sorted values until three or
#' fewer remain. Needs no bandwidth, so the long positive tail of a channel
#' does not blur the mode of its negative population.
#'
#' @param x Numeric vector.
#'
#' @return Numeric, the mode; `NA` for no finite values.
#'
#' @noRd
.fix.half.sample.mode <- function( x ) {

  y <- sort( x[ is.finite( x ) ] )

  if ( length( y ) == 0L ) return( NA_real_ )

  while ( length( y ) > 3L ) {
    n <- length( y )
    h <- ceiling( n / 2 )
    widths <- y[ h:n ] - y[ 1:( n - h + 1L ) ]
    i <- which.min( widths )
    y <- y[ i:( i + h - 1L ) ]
  }

  if ( length( y ) == 3L ) {
    lo <- y[ 2 ] - y[ 1 ]
    hi <- y[ 3 ] - y[ 2 ]
    return( if ( lo < hi ) mean( y[ 1:2 ] ) else if ( hi < lo ) mean( y[ 2:3 ] ) else y[ 2 ] )
  }

  mean( y )
}


#' Upper bound of a sample's own negative population in one channel.
#'
#' The autogating tail threshold: the events below their own positivity
#' boundary `t.event` (so a majority-positive marker still gives the negative
#' mode) are mirrored about their mode, and the bound is the mirrored
#' population's `1 - tail` quantile, `2 mode - quantile(lower half, 2 tail)`.
#' The spread of the negative population in the stained sample -- other dyes'
#' spillover spread, the noise covariance of collinear rows -- is in this
#' bound; the unstained threshold does not see it.
#'
#' @param x Numeric vector, one channel's compensated abundance.
#' @param t.event Numeric vector, the same length, each event's positivity
#'   boundary in that channel.
#' @param tail Numeric in (0, 0.5), the share of the negative population
#'   allowed above the bound.
#' @param max.events Integer, cap on the negative events used.
#' @param min.events Integer, fewest negative events for a bound.
#'
#' @return Numeric, the bound; `-Inf` when too few events are negative.
#'
#' @noRd
.fix.negative.bound <- function( x, t.event, tail, max.events = 50000L,
                                 min.events = 200L ) {

  neg <- x[ x <= t.event ]
  neg <- neg[ is.finite( neg ) ]

  if ( length( neg ) < min.events ) return( -Inf )

  if ( length( neg ) > max.events ) neg <- sample( neg, max.events )

  mode <- .fix.half.sample.mode( neg )
  lower <- neg[ neg <= mode ]

  if ( length( lower ) < 5L ) return( -Inf )

  2 * mode - stats::quantile( lower, probs = min( max( 2 * tail, 0 ), 1 ),
                              names = FALSE )
}


#' Smooth step from 0 at `low` to 1 at `high`.
#' @noRd
.fix.smoothstep <- function( value, low, high ) {
  t <- pmin( pmax( ( value - low ) / max( high - low, .Machine$double.eps ), 0 ), 1 )
  t * t * ( 3 - 2 * t )
}


#' Winning fluorophore per event (0 for none), each abundance scored against
#' its own event's positivity boundary over the channel's dynamic range.
#' @noRd
.fix.dominance <- function( unmixed, threshold.matrix, dyn.range ) {

  score <- sweep( pmax( unmixed - threshold.matrix, 0 ), 2, dyn.range, "/" )

  dominant <- max.col( score, ties.method = "first" )
  dominant[ score[ cbind( seq_along( dominant ), dominant ) ] <= 0 ] <- 0L
  names( dominant ) <- NULL

  dominant
}


#' Clearly positive events of every source fluorophore and the weight they
#' earn.
#'
#' An event is clearly positive for a source when it is above both its
#' spread-scaled unstained boundary and the upper bound of the stained
#' sample's own negative population in that channel
#' (`.fix.negative.bound()`), and (with `source.dominant`) the source is its
#' dominant fluorophore. The tilted tail of a collinear negative population is
#' then never mistaken for positive events. The unstained sample, compensated
#' the same way and judged by the same rule, gives the rate of spurious clear
#' positives (AF outliers, noise tails); `excess` is the count in the stained
#' sample minus that rate times `n.events`. The weight rises smoothly from 0
#' at `source.events[1]` excess events to 1 at `source.events[2]`.
#'
#' @param unmixed.comp,unstained.comp Numeric matrices (events x
#'   fluorophores), the stained and unstained compensated abundance.
#' @param thresholds Named numeric vector, flat positivity thresholds.
#' @param spillover.spread,spread.kappa As in `get.spread.thresholds()`.
#' @param negative.tail,source.events,source.dominant See `fix.my.unmix()`.
#' @param n.events Integer, the stained event count the unstained rate is
#'   scaled to; the whole sample when `unmixed.comp` is a stratified
#'   subsample, which keeps small populations whole.
#'
#' @return A list: `threshold.matrix` and `t.clear` (events x fluorophores),
#'   per fluorophore `bound`, `n.clear`, `baseline`, `excess` and `weight`,
#'   and `dominant` (`NULL` without `source.dominant`).
#'
#' @noRd
.fix.source.evidence <- function( unmixed.comp, unstained.comp, thresholds,
                                  spillover.spread, spread.kappa,
                                  negative.tail, source.events,
                                  source.dominant = TRUE,
                                  n.events = nrow( unmixed.comp ) ) {

  fluorophores <- colnames( unmixed.comp )
  n.fluor <- length( fluorophores )

  tm <- get.spread.thresholds(
    unmixed = unmixed.comp, thresholds = thresholds,
    spillover.spread = spillover.spread, spread.kappa = spread.kappa,
    verbose = FALSE )
  tm0 <- get.spread.thresholds(
    unmixed = unstained.comp, thresholds = thresholds,
    spillover.spread = spillover.spread, spread.kappa = spread.kappa,
    verbose = FALSE )

  dyn.range <- pmax(
    apply( unmixed.comp, 2, stats::quantile, probs = 0.999, names = FALSE ) -
      thresholds[ fluorophores ], .Machine$double.eps )

  bound <- vapply( seq_len( n.fluor ), function( f )
    .fix.negative.bound( unmixed.comp[ , f ], tm[ , f ], negative.tail ),
    numeric( 1 ) )
  names( bound ) <- fluorophores

  bound.row <- matrix( bound, nrow( tm ), n.fluor, byrow = TRUE )
  t.clear <- pmax( tm, bound.row )
  dimnames( t.clear ) <- dimnames( unmixed.comp )

  clear  <- unmixed.comp > t.clear
  clear0 <- unstained.comp > pmax( tm0, matrix( bound, nrow( tm0 ), n.fluor, byrow = TRUE ) )

  dominant <- NULL

  if ( source.dominant ) {
    dominant <- .fix.dominance( unmixed.comp, tm, dyn.range )
    dominant0 <- .fix.dominance( unstained.comp, tm0, dyn.range )
    clear  <- clear  & outer( dominant,  seq_len( n.fluor ), "==" )
    clear0 <- clear0 & outer( dominant0, seq_len( n.fluor ), "==" )
  }

  n.clear  <- stats::setNames( colSums( clear ), fluorophores )
  baseline <- stats::setNames( colMeans( clear0 ) * n.events, fluorophores )
  excess   <- n.clear - baseline

  list(
    threshold.matrix = tm,
    t.clear  = t.clear,
    bound    = bound,
    n.clear  = n.clear,
    baseline = baseline,
    excess   = excess,
    weight   = .fix.smoothstep( excess, source.events[ 1 ], source.events[ 2 ] ),
    dominant = dominant )
}


#' Huber IRLS slope of a line through the origin, with its standard error.
#'
#' Warm-started at `start`; the standard error is the robust residual scale
#' over `sqrt(sum(w x^2))`. Returns slope `0` and standard error `NA` when it
#' does not converge.
#'
#' @noRd
.fix.huber.origin <- function( x, y, start = 0, k = 1.345, max.iter = 100L,
                               tol = 1e-4 ) {

  eps <- .Machine$double.eps
  b <- start

  robust.weights <- function( r ) {
    s <- max( stats::median( abs( r - stats::median( r ) ) ) / 0.6745, eps )
    list( scale = s, w = pmin( 1, k / pmax( abs( r / s ), eps ) ) )
  }

  for ( it in seq_len( as.integer( max.iter ) ) ) {

    rw <- robust.weights( y - b * x )
    sxx <- sum( rw$w * x * x )
    if ( !( sxx > 0 ) ) return( c( slope = 0, se = NA_real_ ) )

    b.new <- sum( rw$w * x * y ) / sxx
    settled <- abs( b.new - b ) < tol * max( abs( b.new ), 1 )
    b <- b.new

    if ( settled ) {
      rw <- robust.weights( y - b * x )
      return( c( slope = b, se = rw$scale / sqrt( max( sum( rw$w * x * x ), eps ) ) ) )
    }
  }

  c( slope = 0, se = NA_real_ )
}


#' Anchored slopes in R, the fallback for `fix_anchored_batch_rcpp()` with
#' identical results. Arrays as in `.fix.anchored.slopes()`, already
#' subsampled.
#' @noRd
.fix.anchored.r <- function( xb, Yb, Tb, Lb, xk, Yk, Tk, start, floor.n,
                             max.coefficient, max.mask.passes, mask.tolerance,
                             k = 1.345, max.iter = 100L, tol = 1e-4 ) {

  m <- ncol( Yb )
  slope <- rep( NA_real_, m )
  se    <- rep( NA_real_, m )
  n.out <- integer( m )

  for ( j in seq_len( m ) ) {

    curr <- start[ j ]
    if ( !is.finite( curr ) || abs( curr ) > max.coefficient ) curr <- 0

    for ( pass in 0:as.integer( max.mask.passes ) ) {

      neg.k   <- ( Yk[ , j ] - curr * xk ) < Tk[ , j ]
      resid.b <- Yb[ , j ] - curr * xb
      neg.b   <- resid.b < Tb[ , j ]
      if ( !is.null( Lb ) && pass > 0L ) neg.b <- neg.b & resid.b > Lb[ , j ]

      if ( sum( neg.k ) < floor.n || sum( neg.b ) < floor.n ) break

      x0 <- stats::median( xk[ neg.k ] )
      y0 <- stats::median( Yk[ neg.k, j ] )

      fit <- .fix.huber.origin( xb[ neg.b ] - x0, Yb[ neg.b, j ] - y0, curr,
                                k = k, max.iter = max.iter, tol = tol )

      n.out[ j ] <- sum( neg.b )
      slope[ j ] <- fit[[ "slope" ]]
      se[ j ]    <- fit[[ "se" ]]

      valid   <- abs( fit[[ "slope" ]] ) <= max.coefficient
      settled <- abs( fit[[ "slope" ]] - curr ) <= mask.tolerance * abs( fit[[ "slope" ]] )
      curr <- fit[[ "slope" ]]

      if ( !valid || settled ) break
    }
  }

  list( slope = slope, n = n.out, se = se )
}


#' Residual spillover slope from one source into every target, anchored at
#' the source-negative bulk and fitted on the source's clearly positive
#' events.
#'
#' The bulk (`x <= t.source`) fixes where the line passes -- its robust centre
#' among events negative for the target -- but not its slope. Its own spread
#' in source and target comes from other dyes' spillover spread, the noise
#' covariance of collinear rows and AF residuals, none of which is spillover
#' from the source, and in a fit over every target-negative event it would
#' dominate a small positive population. The slope is a Huber line through
#' that centre over the clearly positive events (`x > t.source`) that are
#' negative for the target, with negativity re-decided with the current slope
#' removed up to `max.mask.passes` times. With `Threshold.lower`, every pass
#' after the first also drops positive events below the target's lower
#' boundary: a one-sided cut removes only the upper tail of the source's own
#' cell-to-cell spectral variation and so pulls the slope down, while a band
#' around the current line trims both tails. The first pass stays one-sided,
#' so a large negative coefficient is found before the band centres on it.
#'
#' The bulk and the positive events are each subsampled to
#' `max.truncated.events`. Runs in `AutoSpectralRcpp::fix_anchored_batch_rcpp()`
#' when it is installed, else in R, with the same results.
#'
#' @param x Numeric vector (events), the source abundance.
#' @param Y,Threshold.target Numeric matrices (events x targets), the target
#'   abundance and its upper boundary.
#' @param t.source Numeric vector (events), the source's clear-positive
#'   boundary.
#' @param start Numeric vector (targets), warm starts.
#' @param Threshold.lower Numeric matrix (events x targets) or `NULL`, the
#'   target's lower boundary (band).
#'
#' @return A list of per-target vectors: `slope` and `se` (`NA` where fewer
#'   than `floor.n` positive or bulk events are target-negative) and `n`.
#'
#' @noRd
.fix.anchored.slopes <- function( x, Y, Threshold.target, t.source, start,
                                  Threshold.lower      = NULL,
                                  max.truncated.events = 20000L,
                                  floor.n              = 25L,
                                  max.coefficient      = 0.5,
                                  max.mask.passes      = 3L,
                                  mask.tolerance       = 0.05,
                                  huber.k              = 1.345,
                                  huber.max.iter       = 100L,
                                  huber.tol            = 1e-4,
                                  n.threads            = 1L ) {

  m <- ncol( Y )
  floor.n <- max( as.integer( floor.n ), 3L )
  empty <- list( slope = rep( NA_real_, m ), n = integer( m ), se = rep( NA_real_, m ) )

  bright.all <- x > t.source
  bright <- which( bright.all )
  bulk   <- which( !bright.all )

  if ( m == 0L || length( bright ) < floor.n || length( bulk ) < floor.n )
    return( empty )

  if ( length( bulk ) > max.truncated.events )
    bulk <- sort( sample( bulk, max.truncated.events ) )
  if ( length( bright ) > max.truncated.events )
    bright <- sort( sample( bright, max.truncated.events ) )

  xb <- x[ bright ]
  Yb <- Y[ bright, , drop = FALSE ]
  Tb <- Threshold.target[ bright, , drop = FALSE ]
  Lb <- if ( is.null( Threshold.lower ) ) NULL else Threshold.lower[ bright, , drop = FALSE ]
  xk <- x[ bulk ]
  Yk <- Y[ bulk, , drop = FALSE ]
  Tk <- Threshold.target[ bulk, , drop = FALSE ]

  if ( requireNamespace( "AutoSpectralRcpp", quietly = TRUE ) &&
       "fix_anchored_batch_rcpp" %in% ls( getNamespace( "AutoSpectralRcpp" ) ) ) {

    AutoSpectralRcpp::fix_anchored_batch_rcpp(
      xb = xb, Yb = Yb, Tb = Tb, Lb = Lb, xk = xk, Yk = Yk, Tk = Tk,
      start_slope     = as.numeric( start ),
      max_coefficient = max.coefficient,
      max_mask_passes = as.integer( max.mask.passes ),
      mask_tolerance  = mask.tolerance,
      min_count       = floor.n,
      k               = huber.k,
      max_iter        = as.integer( huber.max.iter ),
      tol             = huber.tol,
      n_threads       = as.integer( n.threads ) )

  } else {

    .fix.anchored.r( xb, Yb, Tb, Lb, xk, Yk, Tk, as.numeric( start ), floor.n,
                     max.coefficient, max.mask.passes, mask.tolerance,
                     k = huber.k, max.iter = huber.max.iter, tol = huber.tol )
  }
}


#' Residual spillover matrix by `correct.spectra.glasso()`'s phase one, used
#' as a check on the pair coefficients.
#'
#' Every target's compensated abundance over its own negative events is
#' regressed jointly on every other fluorophore's abundance with an L1
#' penalty (`.glasso.neighborhood.row()`), so a source that only correlates
#' with the target through a third dye has that correlation absorbed by the
#' third dye's coefficient. Negativity is re-decided with the fitted row
#' removed up to `max.mask.passes` times. With `refit`, the selected sources
#' are refitted by least squares, without the penalty's shrinkage.
#' Coefficients below `tolerance` or above `max.coefficient` are zeroed, and a
#' selected source must span `min.span` times its threshold over the
#' target-negative events and rise by `min.rise` noise SDs across it.
#'
#' @return Numeric matrix (source x target), zero diagonal.
#'
#' @noRd
.fix.glasso.matrix <- function( unmixed.comp, threshold.matrix, thresholds,
                                neg.var,
                                min.negative.events  = 200L,
                                max.mask.passes      = 3L,
                                mask.tolerance       = 0.05,
                                max.truncated.events = 20000L,
                                max.coefficient      = 0.5,
                                min.span             = 5,
                                min.rise             = 1,
                                n.lambda             = 40L,
                                lambda.min.ratio     = 1e-3,
                                margin.frac          = 0.05,
                                tolerance            = 1e-4,
                                refit                = TRUE ) {

  fluorophores <- colnames( unmixed.comp )
  out <- matrix( 0, length( fluorophores ), length( fluorophores ),
                 dimnames = list( fluorophores, fluorophores ) )

  for ( target in fluorophores ) {

    sources <- setdiff( fluorophores, target )
    neg.idx <- which( unmixed.comp[ , target ] <= threshold.matrix[ , target ] )

    if ( length( neg.idx ) < min.negative.events ) next

    beta <- NULL

    for ( pass in seq_len( max( as.integer( max.mask.passes ), 1L ) ) ) {

      idx <- if ( length( neg.idx ) > max.truncated.events )
        sort( sample( neg.idx, max.truncated.events ) ) else neg.idx

      row.fit <- .glasso.neighborhood.row(
        abundance = unmixed.comp, target = target, sources = sources,
        neg.idx = idx, n.lambda = n.lambda, lambda.min.ratio = lambda.min.ratio,
        margin.frac = margin.frac, max.events = length( idx ) )

      b <- row.fit$beta[ sources ]
      support <- which( b != 0 )

      if ( refit && length( support ) > 0L ) {
        design <- cbind( 1, unmixed.comp[ idx, sources[ support ], drop = FALSE ] )
        b[ support ] <- stats::lm.fit( design, unmixed.comp[ idx, target ] )$coefficients[ -1L ]
      }

      b[ !is.finite( b ) ] <- 0
      b[ abs( b ) < tolerance ] <- 0
      b[ abs( b ) > max.coefficient ] <- 0
      beta <- b

      fitted <- as.numeric( unmixed.comp[ , sources, drop = FALSE ] %*% b )
      neg.next <- which( ( unmixed.comp[ , target ] - fitted ) <= threshold.matrix[ , target ] )

      if ( length( neg.next ) < min.negative.events ) break

      moved <- length( union( neg.idx, neg.next ) ) - length( intersect( neg.idx, neg.next ) )
      settled <- moved / max( length( neg.idx ), 1L ) <= mask.tolerance
      neg.idx <- neg.next

      if ( settled ) break
    }

    noise <- sqrt( max( neg.var[ target ], 0 ) )

    for ( src in sources[ beta != 0 ] ) {
      span <- if ( length( neg.idx ) > 1L ) diff( range( unmixed.comp[ neg.idx, src ] ) ) else 0
      if ( !( span > min.span * abs( thresholds[ src ] ) &&
              abs( beta[ src ] ) * span >= min.rise * noise ) )
        beta[ src ] <- 0
    }

    out[ sources, target ] <- beta
  }

  out
}


#' How much a spectral change looks like autofluorescence: the largest squared
#' cosine between `change` (detectors) and any row of `af.spectra`; 0 for no
#' change.
#' @noRd
.fix.af.likeness <- function( change, af.spectra ) {

  norm.c <- sqrt( sum( change^2 ) )
  if ( !( norm.c > 0 ) || is.null( af.spectra ) ) return( 0 )

  af.spectra <- as.matrix( af.spectra )
  cosine <- as.numeric( af.spectra %*% change ) /
    pmax( sqrt( rowSums( af.spectra^2 ) ) * norm.c, .Machine$double.eps )

  max( cosine^2 )
}
