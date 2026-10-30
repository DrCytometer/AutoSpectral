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
#' The robust line is either the anchored estimator (`estimator =
#' "anchored"`, `.fix.anchored.slopes()`: the source-negative bulk only fixes
#' the origin and the clearly positive events set the slope, with a band
#' around the current line when `Threshold.lower` is given) or the truncated
#' estimator -- `select.negative()`, the `max.mask.passes` refinement loop,
#' and `.fix.huber.slope()` -- delegated to
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
#' Reproduces `select.negative()`'s bulk cap: every bright (source-positive)
#' event is kept regardless of file size, and the origin-hugging bulk below
#' the source threshold is subsampled down to `max.truncated.events` events
#' total once the negative-selected population exceeds it, delegated to
#' `fix_envelope_truncated_batch_rcpp()`'s own per-target random-number
#' stream. Each target's stream is seeded independently from R's ambient RNG
#' (one integer per target, drawn once before the call), so results are
#' reproducible under `set.seed()` and identical regardless of `n.threads`
#' -- the subsample a target draws depends only on its own seed, never on
#' which thread happened to process it. Not bit-identical to
#' `.fix.envelope.slope()`, which draws its bulk subsample from a single
#' shared RNG stream advanced pair by pair, but statistically equivalent to
#' it.
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
#'   anchor.weight,max.truncated.events,max.coefficient,max.mask.passes,
#'   mask.tolerance As in `.fix.envelope.slope()`.
#' @param start.slope Numeric vector, length `m`, per-target warm starts.
#'   `NULL` starts every target at zero.
#' @param estimator Character, `"anchored"` (default) or `"truncated"`, the
#'   robust line fitted.
#' @param Threshold.lower Numeric matrix, `n x m`, the per-event, per-target
#'   lower boundary, or `NULL`. With the anchored estimator, every pass after
#'   the first also drops positive events below it (band).
#' @param disagreement.reference Character, the slope the binned median is
#'   compared with to give `disagreement`. `"envelope"` (default) is the binned
#'   lower envelope itself, which falls with the spread as the source brightens
#'   and so runs steeper than the line it bounds. `"envelope.corrected"` moves
#'   each bin's envelope back to its centre by the bin's own lower-tail scale,
#'   `(q(envelope.scale.quantile) - q(quantiles[1])) / (z(quantiles[1]) -
#'   z(envelope.scale.quantile))`, which needs no spread model and stays below
#'   most co-expression. `"anchored"` is the fitted robust line, which moves
#'   with any co-expression inside the band, as the median does. `slope` and
#'   `se` are the plain envelope's under every choice.
#' @param envelope.scale.quantile Numeric, the upper quantile of each bin's
#'   lower tail for `"envelope.corrected"`. Default `0.25`.
#' @param n.threads Integer, OpenMP threads for
#'   `fix_envelope_truncated_batch_rcpp()`. Default `1`. Only raise this
#'   when the outer call to `fix.my.unmix()` is not itself running inside
#'   another parallel context (e.g. one sample of several under
#'   `mclapply()`); see that function's `n_threads` documentation.
#'
#' @return A data frame, one row per target in the order of
#'   `colnames(X.target)`, with columns `target`, `slope`,
#'   `slope.truncated`, `slope.anchored`, `se`, `se.fit`, `slope.alt`,
#'   `disagreement`, `coverage`, `n`, `n.positive`, `span`, `noise` -- the
#'   quantities `.fix.envelope.slope()` returns in its list, plus the
#'   anchored slope and its own standard error (`NA` under the truncated
#'   estimator, as `slope.truncated` is under the anchored one) and the
#'   number of positive source events -- or `NULL` if the shared source-side
#'   population is too small to fit anything.
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
    max.truncated.events = 20000L,
    max.coefficient      = 0.2,
    max.mask.passes      = 3L,
    mask.tolerance       = 0.05,
    start.slope          = NULL,
    estimator            = c( "anchored", "truncated" ),
    Threshold.lower      = NULL,
    disagreement.reference  = c( "envelope", "envelope.corrected", "anchored" ),
    envelope.scale.quantile = 0.25,
    n.threads            = 1L ) {

  estimator              <- match.arg( estimator )
  disagreement.reference <- match.arg( disagreement.reference )

  if ( disagreement.reference == "envelope.corrected" &&
       !( quantiles[ 1 ] < envelope.scale.quantile && envelope.scale.quantile < 1 ) )
    stop( "`envelope.scale.quantile` must lie between the envelope quantile and 1.",
          call. = FALSE )

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
    if ( !is.null( Threshold.lower ) )
      Threshold.lower <- Threshold.lower[ keep.event, , drop = FALSE ]

    n <- length( x.source )
  }

  if ( n < min.events ) return( NULL )

  slope.truncated <- rep( NA_real_, m )
  slope.anchored  <- rep( NA_real_, m )
  se.fit          <- rep( NA_real_, m )

  if ( estimator == "anchored" ) {

    anchored <- .fix.anchored.slopes(
      x                    = x.source,
      Y                    = X.target,
      Threshold.target     = Threshold.target,
      t.source             = threshold.source,
      start                = start.slope,
      Threshold.lower      = Threshold.lower,
      max.truncated.events = max.truncated.events,
      floor.n              = min.bin.negative,
      max.coefficient      = max.coefficient,
      max.mask.passes      = max.mask.passes,
      mask.tolerance       = mask.tolerance,
      n.threads            = n.threads )

    slope.anchored <- anchored$slope
    se.fit         <- anchored$se
    n.fit          <- anchored$n
    slope.fit      <- slope.anchored

  } else if (
    requireNamespace( "AutoSpectralRcpp", quietly = TRUE ) &&
    "fix_envelope_truncated_batch_rcpp" %in% ls( getNamespace( "AutoSpectralRcpp" ) )
  ) {

    # One seed per target, drawn from R's ambient RNG stream so the whole
    # pipeline is reproducible under set.seed() upstream of this call.
    seed <- sample.int( .Machine$integer.max, m )

    truncated <- AutoSpectralRcpp::fix_envelope_truncated_batch_rcpp(
      x                     = x.source,
      Y                     = X.target,
      Threshold_target      = Threshold.target,
      threshold_source      = threshold.source,
      start_slope           = start.slope,
      seed                  = as.integer( seed ),
      max_coefficient       = max.coefficient,
      max_mask_passes       = as.integer( max.mask.passes ),
      mask_tolerance        = mask.tolerance,
      min_events            = as.integer( min.events ),
      max_truncated_events  = as.integer( max.truncated.events ),
      n_threads             = as.integer( n.threads ) )

    slope.truncated <- truncated$slope
    n.fit           <- truncated$n
    slope.fit       <- slope.truncated

  } else {
    stop( "Install AutoSpectralRcpp" )
  }

  slope.mask <- ifelse( is.finite( slope.fit ) &
                          abs( slope.fit ) <= max.coefficient,
                        slope.fit, 0 )

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
    slope.anchored  = slope.anchored,
    se              = NA_real_,
    se.fit          = se.fit,
    slope.alt       = NA_real_,
    disagreement    = NA_real_,
    coverage        = coverage,
    n               = n.fit,
    n.positive      = sum( source.positive ),
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

        # vapply() drops the matrix dimension when m == 1 (a single target
        # fluorophore), which would otherwise leave t() returning a 1 x B
        # matrix instead of B x 1 and break every downstream B x m op against
        # sd.bin/weight.* below (the "non-conformable arrays" failure mode).
        # Setting dim() explicitly keeps the m x B shape regardless of m.
        neg.n <- vapply( idx.by.bin, function( idx )
          colSums( target.negative[ idx, , drop = FALSE ] ), numeric( m ) )
        dim( neg.n ) <- c( m, length( bins ) )
        neg.n <- t( neg.n )

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

        negative.quantile <- function( probs ) {
          value <- vapply( idx.by.bin, function( idx )
            vapply( seq_len( m ), function( k ) {
              v <- X.target[ idx, k ][ target.negative[ idx, k ] ]
              if ( length( v ) < min.bin.negative ) return( NA_real_ )
              stats::quantile( v, probs = probs, names = FALSE )
            }, numeric( 1 ) ), numeric( m ) )
          dim( value ) <- c( m, length( bins ) )
          t( value )
        }

        envelope.value <- negative.quantile( quantiles[ 1 ] )

        median.value <- vapply( idx.by.bin, function( idx )
          vapply( seq_len( m ), function( k )
            stats::quantile( X.target[ idx, k ], probs = quantiles[ 2 ],
                             names = FALSE ),
            numeric( 1 ) ), numeric( m ) )
        dim( median.value ) <- c( m, length( bins ) )
        median.value <- t( median.value )

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

        reference.slope <- switch(
          disagreement.reference,
          envelope = envelope$slope,
          anchored = ifelse( is.finite( slope.fit ) &
                               abs( slope.fit ) <= max.coefficient,
                             slope.fit, NA_real_ ),
          envelope.corrected = {
            scale.value <- negative.quantile( envelope.scale.quantile )
            z.tail <- -stats::qnorm( c( quantiles[ 1 ], envelope.scale.quantile ) )
            fit.trace.batch(
              envelope.value + z.tail[ 1 ] * ( scale.value - envelope.value ) /
                ( z.tail[ 1 ] - z.tail[ 2 ] ),
              weight.envelope )$slope
          } )

        disagreement <- abs( reference.slope - compare$slope ) /
          ( abs( reference.slope ) + abs( compare$slope ) + .Machine$double.eps )

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
