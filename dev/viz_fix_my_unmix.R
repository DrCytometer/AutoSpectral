# viz_fix_my_unmix.R
#
# Diagnostic visualizations for fix.my.unmix(). Companion to
# CONTEXT_fix_my_unmix.md.
#
# This is diagnostic-only code, not part of the package, following the same
# local-reimplementation pattern as run_fix_my_unmix_diagnostics_2_1.R and
# _3_1.R: the pair estimator's binning and mask logic are cloned here
# (prefixed fmuviz.*) so the intermediate tables the production function
# does not return - the per-pass mask history and the per-bin stratification
# table - are available for plotting. Exported package functions
# (fit.robust.linear.model, get.spread.thresholds) are called directly
# rather than reimplemented.
#
# Requires AutoSpectral to be loaded (devtools::load_all() or
# library(AutoSpectral)), plus ggplot2.
#
# Two things this script does NOT try to reconstruct after the fact, because
# fix.my.unmix() does not return them: the per-event unmixed abundances and
# the per-pair `noise` used by the rise gate. For those, drive
# fmuviz.pair.stratification() directly on a pair of interest with the same
# unmixed abundances and thresholds you passed into fix.my.unmix(); its
# result carries `noise` alongside everything coefficient.log logs.

# ---------------------------------------------------------------------------
# Pair estimator, instrumented: mask history + stratification bin table
# ---------------------------------------------------------------------------
# Clones .fix.envelope.slope() exactly (including the truncated estimator and
# its iterated negativity mask), but also returns the per-pass mask slopes
# and the per-bin envelope/median table, plus two plots built from them.

fmuviz.pair.stratification <- function(
    x.source, x.target,
    threshold.source, threshold.target,
    spread.var            = 0,
    neg.var                = 125^2,
    source.mask            = NULL,
    quantiles              = c( 0.05, 0.5 ),
    n.levels               = 10L,
    min.events             = 200L,
    min.bin.negative       = 25L,
    spread.addback         = FALSE,
    anchor.weight          = 1,
    max.truncated.events   = 20000L,
    max.coefficient        = 0.2,
    max.mask.passes        = 3L,
    mask.tolerance         = 0.05,
    source.name            = "source",
    target.name            = "target"
) {

  n <- length( x.source )
  if ( n < min.events )
    stop( "Fewer than `min.events` events for this pair.", call. = FALSE )

  if ( length( threshold.source ) == 1L )
    threshold.source <- rep( threshold.source, n )
  if ( length( threshold.target ) == 1L )
    threshold.target <- rep( threshold.target, n )

  if ( !is.null( source.mask ) ) {
    keep.event <- source.mask | x.source <= threshold.source
    x.source         <- x.source[ keep.event ]
    x.target         <- x.target[ keep.event ]
    threshold.source <- threshold.source[ keep.event ]
    threshold.target <- threshold.target[ keep.event ]
    n <- length( x.source )
    if ( n < min.events )
      stop( "Fewer than `min.events` events survive `source.mask`.", call. = FALSE )
  }

  select.negative <- function( slope ) {
    index <- which( x.target - slope * x.source < threshold.target )
    if ( length( index ) <= max.truncated.events ) return( index )
    bright <- index[ x.source[ index ] > threshold.source[ index ] ]
    bulk   <- setdiff( index, bright )
    n.bulk <- max( max.truncated.events - length( bright ), min.events )
    if ( length( bulk ) > n.bulk ) bulk <- sample( bulk, n.bulk )
    c( bright, bulk )
  }

  fit.truncated <- function( index ) {
    if ( length( index ) < min.events ) return( NA_real_ )
    fit.robust.linear.model(
      x.data = x.source[ index ], y.data = x.target[ index ],
      x.name = source.name, y.name = target.name, fix.unmix = TRUE )[ 2 ]
  }

  truncated.index <- select.negative( 0 )
  slope.truncated <- fit.truncated( truncated.index )

  mask.history <- data.frame(
    pass = 0L, slope = slope.truncated, n.negative = length( truncated.index ) )

  for ( pass in seq_len( as.integer( max.mask.passes ) ) ) {

    if ( !is.finite( slope.truncated ) || abs( slope.truncated ) > max.coefficient ) break

    index.next <- select.negative( slope.truncated )
    slope.next <- fit.truncated( index.next )
    if ( !is.finite( slope.next ) ) break

    settled <- abs( slope.next - slope.truncated ) <= mask.tolerance * abs( slope.next )

    slope.truncated <- slope.next
    truncated.index  <- index.next

    mask.history <- rbind( mask.history, data.frame(
      pass = pass, slope = slope.truncated, n.negative = length( truncated.index ) ) )

    if ( settled ) break
  }

  slope.mask <- if ( is.finite( slope.truncated ) && abs( slope.truncated ) <= max.coefficient )
    slope.truncated else 0

  target.negative <- x.target - slope.mask * x.source < threshold.target
  source.positive <- x.source > threshold.source

  coverage <- 1
  if ( sum( source.positive ) >= min.bin.negative ) {
    bright.index <- which( source.positive )
    bright.cut   <- stats::quantile( x.source[ bright.index ], probs = 2 / 3, names = FALSE )
    top.index    <- bright.index[ x.source[ bright.index ] >= bright.cut ]
    if ( length( top.index ) >= min.bin.negative )
      coverage <- mean( target.negative[ top.index ] )
  }

  span.truncated <- if ( length( truncated.index ) > 1L )
    diff( range( x.source[ truncated.index ] ) ) else 0

  # --- mask evolution plot --------------------------------------------------
  # The boundary lines use the median of threshold.target for display: the
  # true boundary is per-event (and spread-scaled when threshold.target came
  # from get.spread.thresholds()); see fmuviz.spread.boundary() for that.

  mask.plot.data <- data.frame(
    x.source = x.source, x.target = x.target, target.negative = target.negative )

  if ( nrow( mask.plot.data ) > 20000L )
    mask.plot.data <- mask.plot.data[ sample.int( nrow( mask.plot.data ), 20000L ), ]

  boundary.lines <- data.frame(
    pass      = mask.history$pass,
    intercept = stats::median( threshold.target ),
    slope     = mask.history$slope )

  plot.mask <- ggplot2::ggplot( mask.plot.data,
      ggplot2::aes( x = x.source, y = x.target, colour = target.negative ) ) +
    ggplot2::geom_point( size = 0.4, alpha = 0.35 ) +
    ggplot2::geom_abline(
      data = boundary.lines,
      ggplot2::aes( intercept = intercept, slope = slope, linetype = factor( pass ) ),
      colour = "black", linewidth = 0.6 ) +
    ggplot2::scale_colour_manual(
      values = c( `TRUE` = "steelblue", `FALSE` = "grey70" ),
      labels = c( `TRUE` = "target-negative (final mask)",
                 `FALSE` = "target-positive (final mask)" ), name = NULL ) +
    ggplot2::labs(
      title    = sprintf( "%s -> %s: negativity mask across %d pass(es)",
                          source.name, target.name, max( mask.history$pass ) ),
      subtitle = sprintf(
        "Final slope %.4f (started at 0); boundary lines use the median threshold for display only",
        slope.mask ),
      x = sprintf( "%s abundance", source.name ), y = sprintf( "%s abundance", target.name ),
      linetype = "Mask pass" ) +
    ggplot2::theme_minimal() +
    ggplot2::theme( legend.position = "bottom" )

  truncated.only <- function() {
    list( slope = NA_real_, slope.truncated = slope.truncated, se = NA_real_,
          slope.alt = NA_real_, disagreement = NA_real_, coverage = coverage,
          n = n, span = span.truncated, noise = sqrt( max( neg.var, 0 ) ),
          mask.history = mask.history, bin.table = NULL,
          plot.mask = plot.mask, plot.strata = NULL )
  }

  if ( sum( source.positive ) < min.events ) return( truncated.only() )
  if ( sum( !source.positive ) < min.bin.negative ) return( truncated.only() )

  brk <- unique( stats::quantile(
    x.source[ source.positive ], probs = seq( 0, 1, length.out = n.levels + 1 ), names = FALSE ) )
  if ( length( brk ) < 4 ) return( truncated.only() )

  bin <- integer( n )
  bin[ source.positive ] <- as.integer( cut(
    x.source[ source.positive ], breaks = brk, include.lowest = TRUE ) )

  bins <- sort( unique( bin ) )
  if ( length( bins ) < 5 ) return( truncated.only() )

  centre <- vapply( bins, function( b ) stats::median( x.source[ bin == b ] ), numeric( 1 ) )
  neg.n  <- vapply( bins, function( b ) sum( target.negative[ bin == b ] ), integer( 1 ) )
  bin.n  <- vapply( bins, function( b ) sum( bin == b ), integer( 1 ) )

  usable <- neg.n >= min.bin.negative
  if ( sum( usable ) < 4L ) return( truncated.only() )

  sd.bin <- sqrt( pmax( neg.var + spread.var * pmax( centre, 0 ), 0 ) )
  z      <- if ( spread.addback ) -stats::qnorm( quantiles ) else c( 0, 0 )

  anchor.cap <- function( w ) {
    if ( !any( bins > 0L ) || !any( bins == 0L ) ) return( w )
    w[ bins == 0L ] <- pmin( w[ bins == 0L ], anchor.weight * max( w[ bins > 0L ] ) )
    w
  }

  weight.envelope <- anchor.cap( neg.n / pmax( sd.bin^2, .Machine$double.eps ) )
  weight.compare  <- anchor.cap( bin.n / pmax( sd.bin^2, .Machine$double.eps ) )

  fit.trace <- function( value, weight ) {
    keep <- usable & is.finite( value )
    if ( sum( keep ) < 4L )
      return( list( slope = NA_real_, se = NA_real_, fitted = rep( NA_real_, length( value ) ) ) )
    x.bin  <- centre[ keep ]
    y.bin  <- value[ keep ]
    w.bin  <- weight[ keep ]
    fit    <- stats::lm.wfit( x = cbind( 1, x.bin ), y = y.bin, w = w.bin )
    x.mean <- sum( w.bin * x.bin ) / sum( w.bin )
    sxx    <- sum( w.bin * ( x.bin - x.mean )^2 )
    dof    <- sum( keep ) - 2L
    se <- if ( dof > 0 && sxx > 0 )
      sqrt( sum( w.bin * fit$residuals^2 ) / dof / sxx ) else NA_real_
    fitted.full <- rep( NA_real_, length( value ) )
    fitted.full[ keep ] <- fit$fitted.values
    list( slope = unname( fit$coefficients[ 2 ] ), se = se, fitted = fitted.full )
  }

  envelope.value <- vapply( bins, function( b ) {
    v <- x.target[ bin == b & target.negative ]
    if ( length( v ) < min.bin.negative ) return( NA_real_ )
    stats::quantile( v, probs = quantiles[ 1 ], names = FALSE )
  }, numeric( 1 ) )

  median.value <- vapply( bins, function( b )
    stats::quantile( x.target[ bin == b ], probs = quantiles[ 2 ], names = FALSE ), numeric( 1 ) )

  envelope <- fit.trace( envelope.value + z[ 1 ] * sd.bin, weight.envelope )
  compare  <- fit.trace( median.value   + z[ 2 ] * sd.bin, weight.compare )

  disagreement <- abs( envelope$slope - compare$slope ) /
    ( abs( envelope$slope ) + abs( compare$slope ) + .Machine$double.eps )

  bin.table <- data.frame(
    bin = bins, centre = centre, neg.n = neg.n, bin.n = bin.n, usable = usable,
    envelope.value = envelope.value, median.value = median.value,
    envelope.fitted = envelope$fitted, compare.fitted = compare$fitted )

  anchor.row <- bin.table[ bin.table$bin == 0L, , drop = FALSE ]

  plot.strata <- ggplot2::ggplot( bin.table, ggplot2::aes( x = centre ) ) +
    ggplot2::geom_point( ggplot2::aes( y = envelope.value, size = neg.n ),
                         colour = "firebrick", na.rm = TRUE ) +
    ggplot2::geom_line( ggplot2::aes( y = envelope.fitted ), colour = "firebrick",
                        linetype = "dashed", na.rm = TRUE ) +
    ggplot2::geom_point( ggplot2::aes( y = median.value, size = bin.n ),
                         colour = "grey40", na.rm = TRUE ) +
    ggplot2::geom_line( ggplot2::aes( y = compare.fitted ), colour = "grey40",
                        linetype = "dashed", na.rm = TRUE ) +
    ggplot2::geom_point( data = anchor.row, ggplot2::aes( y = envelope.value ),
                         colour = "black", shape = 8, size = 3 ) +
    ggplot2::labs(
      title    = sprintf( "%s -> %s: abundance strata", source.name, target.name ),
      subtitle = sprintf(
        "Envelope slope %.4f (red), median slope %.4f (grey), disagreement %.2f; star = anchor bin",
        envelope$slope, compare$slope, disagreement ),
      x = sprintf( "%s abundance (bin median)", source.name ),
      y = sprintf( "%s abundance", target.name ), size = "Events" ) +
    ggplot2::theme_minimal() +
    ggplot2::theme( legend.position = "bottom" )

  list(
    slope = envelope$slope, slope.truncated = slope.truncated, se = envelope$se,
    slope.alt = compare$slope, disagreement = disagreement, coverage = coverage,
    n = n, span = max( centre ) - min( centre ), noise = sqrt( max( neg.var, 0 ) ),
    mask.history = mask.history, bin.table = bin.table,
    plot.mask = plot.mask, plot.strata = plot.strata )
}


# ---------------------------------------------------------------------------
# Spread-scaled positivity boundary
# ---------------------------------------------------------------------------
# Calls the exported get.spread.thresholds() directly and plots the boundary
# it produces for one target channel against a chosen source's abundance,
# next to the flat threshold it replaces. Because the spread term sums the
# contribution of every fluorophore present in each event
# (CONTEXT math: t_c,b = m*t_b + kappa*sqrt(sum_a SS_a,b * x_c,a)), the curve
# traced here reflects the whole event's spread, not `source` alone; it will
# look like a clean sqrt curve only where `source` dominates the spread into
# `target` for the events plotted.

fmuviz.spread.boundary <- function(
    unmixed, thresholds, spillover.spread, spread.kappa = 2,
    source, target, margin = 1
) {

  threshold.matrix <- get.spread.thresholds(
    unmixed = unmixed, thresholds = thresholds, spillover.spread = spillover.spread,
    spread.kappa = spread.kappa, margin = margin, verbose = FALSE )

  x.source    <- unmixed[ , source ]
  y.target    <- unmixed[ , target ]
  y.threshold <- threshold.matrix[ , target ]
  flat        <- thresholds[ target ] * margin

  ord <- order( x.source )
  curve.data <- data.frame( x.source = x.source[ ord ], y.threshold = y.threshold[ ord ] )
  curve.data <- curve.data[ !duplicated( round( curve.data$x.source, 6 ) ), ]

  event.data <- data.frame(
    x.source = x.source, y.target = y.target, clears.spread = y.target > y.threshold )

  if ( nrow( event.data ) > 20000L )
    event.data <- event.data[ sample.int( nrow( event.data ), 20000L ), ]

  ggplot2::ggplot() +
    ggplot2::geom_point( data = event.data,
        ggplot2::aes( x = x.source, y = y.target, colour = clears.spread ),
        size = 0.4, alpha = 0.35 ) +
    ggplot2::geom_hline( yintercept = flat, linetype = "dotted", colour = "grey30" ) +
    ggplot2::geom_line( data = curve.data,
        ggplot2::aes( x = x.source, y = y.threshold ), colour = "black", linewidth = 0.8 ) +
    ggplot2::scale_colour_manual(
      values = c( `TRUE` = "steelblue", `FALSE` = "grey70" ),
      labels = c( `TRUE` = "clears spread-scaled boundary",
                 `FALSE` = "below spread-scaled boundary" ), name = NULL ) +
    ggplot2::labs(
      title    = sprintf( "%s spread into %s: spread-scaled vs flat boundary", source, target ),
      subtitle = sprintf(
        "Flat threshold (dotted) = %.1f; spread-scaled boundary (solid) widens with total spread variance, kappa = %.1f",
        flat, spread.kappa ),
      x = sprintf( "%s abundance", source ), y = sprintf( "%s abundance", target ) ) +
    ggplot2::theme_minimal() +
    ggplot2::theme( legend.position = "bottom" )
}


# ---------------------------------------------------------------------------
# Phase one gate scorecard — from coefficient.log, no recomputation
# ---------------------------------------------------------------------------
# coefficient.log does not carry `noise`, so the rise gate
# (abs(slope) * span >= min.rise * noise) cannot be shown here; drive
# fmuviz.pair.stratification() on a specific pair for that.

fmuviz.pair.gate.scorecard <- function(
    coefficient.log, min.negative.frac = 0.10, max.coefficient = 0.2,
    max.disagreement = 0.5
) {

  d <- coefficient.log
  d$pair      <- paste( d$source, "->", d$channel )
  d$accepted  <- d$trust > 0
  d$slope.use <- ifelse( is.na( d$slope.truncated ), d$slope, d$slope.truncated )

  long <- rbind(
    data.frame( pair = d$pair, accepted = d$accepted, gate = "coverage",
               value = d$coverage, reference = min.negative.frac ),
    data.frame( pair = d$pair, accepted = d$accepted, gate = "|coefficient|",
               value = abs( d$slope.use ), reference = max.coefficient ),
    data.frame( pair = d$pair, accepted = d$accepted, gate = "disagreement",
               value = d$disagreement, reference = max.disagreement ),
    data.frame( pair = d$pair, accepted = d$accepted, gate = "span",
               value = d$span, reference = NA_real_ )
  )

  ggplot2::ggplot( long, ggplot2::aes( x = value, y = pair, colour = accepted ) ) +
    ggplot2::geom_vline( ggplot2::aes( xintercept = reference ), linetype = "dashed",
                         colour = "grey40", na.rm = TRUE ) +
    ggplot2::geom_point( size = 1.6, na.rm = TRUE ) +
    ggplot2::facet_wrap( ~ gate, scales = "free_x" ) +
    ggplot2::scale_colour_manual( values = c( `TRUE` = "steelblue", `FALSE` = "grey70" ) ) +
    ggplot2::labs(
      title    = "Phase one pair gates",
      subtitle = paste(
        "Coverage and disagreement are from the final iteration;",
        "span has no fixed reference (compare to min.span x |threshold|);",
        "the rise gate needs `noise` - see fmuviz.pair.stratification()." ),
      x = NULL, y = NULL, colour = "Coefficient fitted" ) +
    ggplot2::theme_minimal() +
    ggplot2::theme( axis.text.y = ggplot2::element_text( size = 7 ) )
}


# ---------------------------------------------------------------------------
# Phase two gate scorecard — from signature.log, no recomputation
# ---------------------------------------------------------------------------

fmuviz.signature.gate.scorecard <- function(
    signature.log,
    min.explained = 0.8, max.explained = 1.2, max.resid = 0.03,
    max.intercept = 0.03, min.bg.align = -0.9, max.clamp.frac = 0.15,
    max.angle = 10, max.anchor = 0.10, max.vif = 500
) {

  d <- signature.log

  long <- rbind(
    data.frame( fluorophore = d$fluorophore, accepted = d$accepted,
               gate = "explained.total", value = d$explained.total,
               lower = min.explained, upper = max.explained ),
    data.frame( fluorophore = d$fluorophore, accepted = d$accepted,
               gate = "resid.rel", value = d$resid.rel,
               lower = NA_real_, upper = max.resid ),
    data.frame( fluorophore = d$fluorophore, accepted = d$accepted,
               gate = "intercept.rel", value = d$intercept.rel,
               lower = NA_real_, upper = max.intercept ),
    data.frame( fluorophore = d$fluorophore, accepted = d$accepted,
               gate = "bg.align", value = d$bg.align,
               lower = min.bg.align, upper = NA_real_ ),
    data.frame( fluorophore = d$fluorophore, accepted = d$accepted,
               gate = "clamp.frac", value = d$clamp.frac,
               lower = NA_real_, upper = max.clamp.frac ),
    data.frame( fluorophore = d$fluorophore, accepted = d$accepted,
               gate = "deg.change", value = d$deg.change,
               lower = NA_real_, upper = max.angle ),
    data.frame( fluorophore = d$fluorophore, accepted = d$accepted,
               gate = "anchor.rel", value = d$anchor.rel,
               lower = NA_real_, upper = max.anchor ),
    data.frame( fluorophore = d$fluorophore, accepted = d$accepted,
               gate = "vif.target", value = d$vif.target,
               lower = NA_real_, upper = max.vif )
  )

  ggplot2::ggplot( long, ggplot2::aes( x = value, y = fluorophore, colour = accepted ) ) +
    ggplot2::geom_vline( ggplot2::aes( xintercept = lower ), linetype = "dashed",
                         colour = "grey40", na.rm = TRUE ) +
    ggplot2::geom_vline( ggplot2::aes( xintercept = upper ), linetype = "dashed",
                         colour = "grey40", na.rm = TRUE ) +
    ggplot2::geom_point( size = 1.8, na.rm = TRUE ) +
    ggplot2::facet_wrap( ~ gate, scales = "free_x" ) +
    ggplot2::scale_colour_manual( values = c( `TRUE` = "steelblue", `FALSE` = "firebrick" ) ) +
    ggplot2::labs(
      title    = "Phase two signature gates",
      subtitle = "Dashed lines mark the acceptance bounds; colour is the final accepted flag",
      x = NULL, y = NULL, colour = "Accepted" ) +
    ggplot2::theme_minimal()
}


fmuviz.rejection.summary <- function( signature.log ) {

  tab <- as.data.frame( table( reason = signature.log$reason ) )

  ggplot2::ggplot( tab, ggplot2::aes( x = stats::reorder( reason, -Freq ), y = Freq ) ) +
    ggplot2::geom_col( fill = "steelblue" ) +
    ggplot2::labs( title = "Phase two: rejection reasons", x = NULL, y = "Fluorophores" ) +
    ggplot2::theme_minimal() +
    ggplot2::theme( axis.text.x = ggplot2::element_text( angle = 45, hjust = 1 ) )
}


# ---------------------------------------------------------------------------
# Example usage (not executed)
# ---------------------------------------------------------------------------
#
# result <- fix.my.unmix(
#   spectra = spectra, unstained.sample = unstained.file,
#   fully.stained.sample = stained.file, flow.control = flow.control,
#   asp = asp, variants = variants )
#
# fmuviz.pair.gate.scorecard( result$coefficient.log )
# fmuviz.signature.gate.scorecard( result$signature.log )
# fmuviz.rejection.summary( result$signature.log )
#
# To drill into one pair's stratification and mask evolution, supply the
# same unmixed abundances and thresholds fix.my.unmix() would have used
# (e.g. from unmix.ols.fast(raw.data, spectra) and get.spread.thresholds()):
#
# strat <- fmuviz.pair.stratification(
#   x.source          = unmixed.comp[ , "PE" ],
#   x.target          = unmixed.comp[ , "Spark Blue 550" ],
#   threshold.source  = threshold.matrix[ , "PE" ],
#   threshold.target  = threshold.matrix[ , "Spark Blue 550" ],
#   spread.var        = variants$spillover.spread[ "PE", "Spark Blue 550" ],
#   neg.var           = 125^2,
#   source.name       = "PE", target.name = "Spark Blue 550" )
# strat$plot.mask
# strat$plot.strata
#
# fmuviz.spread.boundary(
#   unmixed = unmixed.comp, thresholds = thresholds,
#   spillover.spread = variants$spillover.spread, spread.kappa = 2,
#   source = "PE", target = "Spark Blue 550" )
