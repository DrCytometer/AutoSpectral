# fix_my_unmix_vectorized.R

#' Batched multi-target version of `.fix.envelope.slope()`, sharing one
#' source fluorophore's abundance across every target it might spill into.
#'
#' Everything in the per-pair estimator that depends only on the source --
#' the dominance filter, the source-positivity split, and the abundance bins
#' and their breakpoints -- is identical for every target sharing that
#' source, so it is computed once here instead of once per pair. Only the
#' per-target pieces (the target's own threshold, its own negative-
#' population variance, its own spillover-spread term, and the robust slope
#' itself) vary by column.
#'
#' The truncated estimator -- `select.negative()`, the `max.mask.passes`
#' refinement loop, and `.fix.huber.slope()` -- is delegated to
#' `fix_envelope_truncated_batch_rcpp()` for every target sharing this
#' source in one call. This piece does not have a profitable pure-R batched
#' form: batching across targets does not reduce the number of `median()`
#' calls or the arithmetic they sit inside (still one Huber IRLS per target
#' per pass, exactly as many as a per-pair loop would run), and R's
#' `sweep()` / `outer()` reallocate a working-set-sized object on every one
#' of up to `targets x mask.passes x irls.iterations` steps, which a per-pair
#' loop never needed to do at all. The only way to make batching pay off here
#' is to remove that allocation and the R-level `median()` dispatch
#' entirely, which needs compiled code. Everything else in this function --
#' the binned envelope estimator -- was never the bottleneck (it runs on
#' `n.levels.pair` bins, not the full event count) and stays in R.
#'
#' Unlike `.fix.envelope.slope()`, the truncated estimator here does not
#' reproduce the bulk subsample when a target's negative population exceeds
#' `max.truncated.events`; it always fits the full negative-selected
#' population. Per that function's own documented rationale, the bulk sits
#' at the origin and carries no leverage on the slope, so this is expected
#' to be equivalent, not an approximation -- but it means a target whose
#' subsample would have triggered will not be bit-identical to
#' `.fix.envelope.slope()`, only statistically equivalent to it.
#'
#' @param x.source Numeric vector, length `n`, the source fluorophore's
#'   current compensated abundance.
#' @param X.target Numeric matrix, `n x m`, one column per target
#'   fluorophore, column-named by fluorophore.
#' @param threshold.source Numeric vector, length `n`, the source's own
#'   per-event positivity boundary.
#' @param Threshold.target Numeric matrix, `n x m`, the per-event,
#'   per-target positivity boundary.
#' @param spread.var Numeric vector, length `m`, this source's contribution
#'   to each target's spillover-spread variance.
#' @param neg.var Numeric vector, length `m`, each target's negative-
#'   population variance.
#' @param source.mask Logical vector, length `n`, or `NULL`. Shared across
#'   every target, since it depends only on the source.
#' @param quantiles,n.levels,min.events,min.bin.negative,spread.addback,
#'   anchor.weight,max.coefficient,max.mask.passes,mask.tolerance
#'   As in `.fix.envelope.slope()`.
#' @param start.slope Numeric vector, length `m`, per-target warm starts.
#'   `NULL` starts every target at zero.
#' @param n.threads Integer, OpenMP threads for
#'   `fix_envelope_truncated_batch_rcpp()`. Default `1`. Only raise this
#'   when the outer call to `fix.my.unmix()` is not itself running inside
#'   another parallel context (e.g. one sample of several under
#'   `mclapply()`); see that function's `n_threads` documentation.
#'
#' @return A data frame, one row per target in the order of
#'   `colnames(X.target)`, with columns `target`, `slope`,
#'   `slope.truncated`, `se`, `slope.alt`, `disagreement`, `coverage`, `n`,
#'   `span`, `noise` -- the same quantities `.fix.envelope.slope()` returns
#'   in its list -- or `NULL` if the shared source-side population is too
#'   small to fit anything.
#'
#' @noRd
.fix.envelope.slope.batch <- function(
    x.source, X.target, threshold.source, Threshold.target,
    spread.var, neg.var,
    source.mask          = NULL,
    quantiles            = c( 0.05, 0.5 ),
    n.levels             = 10L,
    min.events           = 200L,
    min.bin.negative     = 25L,
    spread.addback       = TRUE,
    anchor.weight        = 1,
    max.coefficient      = 0.2,
    max.mask.passes      = 3L,
    mask.tolerance       = 0.05,
    start.slope          = NULL,
    n.threads            = 1L ) {

  n <- length( x.source )
  m <- ncol( X.target )

  if ( is.null( start.slope ) ) start.slope <- rep( 0, m )

  if ( length( threshold.source ) == 1L )
    threshold.source <- rep( threshold.source, n )

  if ( !is.null( source.mask ) ) {

    keep.event <- source.mask | x.source <= threshold.source

    x.source         <- x.source[ keep.event ]
    X.target         <- X.target[ keep.event, , drop = FALSE ]
    threshold.source <- threshold.source[ keep.event ]
    Threshold.target <- Threshold.target[ keep.event, , drop = FALSE ]

    n <- length( x.source )
  }

  if ( n < min.events ) return( NULL )

  if (
    requireNamespace( "AutoSpectralRcpp", quietly = TRUE ) &&
    "fix_envelope_truncated_batch_rcpp" %in% ls( getNamespace( "AutoSpectralRcpp" ) )
  ) {
    truncated <- AutoSpectralRcpp::fix_envelope_truncated_batch_rcpp(
      x                = x.source,
      Y                = X.target,
      Threshold_target = Threshold.target,
      start_slope      = start.slope,
      max_coefficient  = max.coefficient,
      max_mask_passes  = as.integer( max.mask.passes ),
      mask_tolerance   = mask.tolerance,
      min_events       = as.integer( min.events ),
      n_threads        = as.integer( n.threads ) )
  } else {
    stop( "Install AutoSpectralRcpp" )
  }

  slope.truncated <- truncated$slope

  slope.mask <- ifelse( is.finite( slope.truncated ) &
                          abs( slope.truncated ) <= max.coefficient,
                        slope.truncated, 0 )

  # Recomputed once, cheaply, from the final slope -- not inside the pass
  # loop, so this is a single O(n x m) comparison, not a hot path.
  Neg <- ( X.target - outer( x.source, slope.mask ) ) < Threshold.target
  target.negative <- Neg
  source.positive <- x.source > threshold.source

  coverage <- rep( 1, m )

  if ( sum( source.positive ) >= min.bin.negative ) {

    bright.index <- which( source.positive )
    bright.cut   <- stats::quantile( x.source[ bright.index ], probs = 2 / 3,
                                     names = FALSE )
    top.index    <- bright.index[ x.source[ bright.index ] >= bright.cut ]

    if ( length( top.index ) >= min.bin.negative )
      coverage <- colMeans( target.negative[ top.index, , drop = FALSE ] )
  }

  span.truncated <- vapply( seq_len( m ), function( k ) {
    idx <- which( Neg[ , k ] )
    if ( length( idx ) > 1L ) diff( range( x.source[ idx ] ) ) else 0
  }, numeric( 1 ) )

  result <- data.frame(
    target          = colnames( X.target ),
    slope           = NA_real_,
    slope.truncated = slope.truncated,
    se              = NA_real_,
    slope.alt       = NA_real_,
    disagreement    = NA_real_,
    coverage        = coverage,
    n               = truncated$n,
    span            = span.truncated,
    noise           = sqrt( pmax( neg.var, 0 ) ),
    stringsAsFactors = FALSE )

  have.bins <- sum( source.positive ) >= min.events &&
    sum( !source.positive ) >= min.bin.negative

  if ( have.bins ) {

    brk <- unique( stats::quantile(
      x.source[ source.positive ],
      probs = seq( 0, 1, length.out = n.levels + 1 ), names = FALSE ) )

    if ( length( brk ) >= 4 ) {

      bin <- integer( n )
      bin[ source.positive ] <- as.integer( cut(
        x.source[ source.positive ], breaks = brk, include.lowest = TRUE ) )

      bins <- sort( unique( bin[ bin > 0L ] ) )

      if ( length( bins ) >= 5 ) {

        idx.by.bin <- split( which( bin > 0L ), bin[ bin > 0L ] )
        idx.by.bin <- idx.by.bin[ as.character( bins ) ]

        centre <- vapply( idx.by.bin, function( idx )
          stats::median( x.source[ idx ] ), numeric( 1 ) )
        bin.n  <- lengths( idx.by.bin )

        neg.n <- t( vapply( idx.by.bin, function( idx )
          colSums( target.negative[ idx, , drop = FALSE ] ), numeric( m ) ) )

        usable <- neg.n >= min.bin.negative

        sd.bin <- sqrt( pmax(
          outer( pmax( centre, 0 ), spread.var ) +
            matrix( neg.var, length( bins ), m, byrow = TRUE ), 0 ) )

        z <- if ( spread.addback ) -stats::qnorm( quantiles ) else c( 0, 0 )

        anchor.cap <- function( w ) {
          if ( !any( bins > 0L ) || !any( bins == 0L ) ) return( w )
          cap  <- anchor.weight * apply( w[ bins > 0L, , drop = FALSE ], 2, max )
          idx0 <- which( bins == 0L )
          w[ idx0, ] <- pmin( w[ idx0, , drop = FALSE ],
                              matrix( cap, length( idx0 ), m, byrow = TRUE ) )
          w
        }

        weight.envelope <- anchor.cap( neg.n / pmax( sd.bin^2, .Machine$double.eps ) )
        weight.compare  <- anchor.cap(
          matrix( bin.n, length( bins ), m ) /
            pmax( sd.bin^2, .Machine$double.eps ) )

        envelope.value <- t( vapply( idx.by.bin, function( idx )
          vapply( seq_len( m ), function( k ) {
            v <- X.target[ idx, k ][ target.negative[ idx, k ] ]
            if ( length( v ) < min.bin.negative ) return( NA_real_ )
            stats::quantile( v, probs = quantiles[ 1 ], names = FALSE )
          }, numeric( 1 ) ), numeric( m ) ) )

        median.value <- t( vapply( idx.by.bin, function( idx )
          vapply( seq_len( m ), function( k )
            stats::quantile( X.target[ idx, k ], probs = quantiles[ 2 ],
                             names = FALSE ),
            numeric( 1 ) ), numeric( m ) ) )

        fit.trace.batch <- function( value, weight ) {

          slope.out <- rep( NA_real_, m )
          se.out    <- rep( NA_real_, m )

          for ( k in seq_len( m ) ) {

            keep <- usable[ , k ] & is.finite( value[ , k ] )
            if ( sum( keep ) < 4L ) next

            x.bin <- centre[ keep ]
            y.bin <- value[ keep, k ]
            w.bin <- weight[ keep, k ]

            x.mean <- sum( w.bin * x.bin ) / sum( w.bin )
            y.mean <- sum( w.bin * y.bin ) / sum( w.bin )
            sxx    <- sum( w.bin * ( x.bin - x.mean )^2 )
            sxy    <- sum( w.bin * ( x.bin - x.mean ) * ( y.bin - y.mean ) )

            if ( sxx <= 0 ) next

            slope.k     <- sxy / sxx
            intercept.k <- y.mean - slope.k * x.mean
            residuals.k <- y.bin - ( intercept.k + slope.k * x.bin )

            dof <- sum( keep ) - 2L

            slope.out[ k ] <- slope.k
            se.out[ k ]    <- if ( dof > 0 && sxx > 0 )
              sqrt( sum( w.bin * residuals.k^2 ) / dof / sxx ) else NA_real_
          }

          list( slope = slope.out, se = se.out )
        }

        envelope <- fit.trace.batch( envelope.value + z[ 1 ] * sd.bin,
                                     weight.envelope )
        compare  <- fit.trace.batch( median.value   + z[ 2 ] * sd.bin,
                                     weight.compare )

        disagreement <- abs( envelope$slope - compare$slope ) /
          ( abs( envelope$slope ) + abs( compare$slope ) + .Machine$double.eps )

        result$slope        <- envelope$slope
        result$se           <- envelope$se
        result$slope.alt    <- compare$slope
        result$disagreement <- disagreement
        result$span         <- max( centre ) - min( centre )
      }
    }
  }

  result
}
