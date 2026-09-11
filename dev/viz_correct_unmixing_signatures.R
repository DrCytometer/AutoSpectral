# viz_correct_unmixing_signatures.R
#
# Diagnostic visualizations for correct.unmixing.signatures(). Companion to
# CONTEXT_signature_correction2.md.
#
# This is diagnostic-only code, not part of the package. It reimplements the
# small internal steps of correct.unmixing.signatures() locally (prefixed
# .sigviz.*) so every plot is self-contained and traceable to a named object
# in scope, following the same pattern as the run_signature_correction_*.R
# harnesses. It calls the package's exported functions directly
# (unmix.ols.fast, get.spread.thresholds) rather than reimplementing them.
#
# Requires AutoSpectral to be loaded (devtools::load_all() or
# library(AutoSpectral)), plus ggplot2, MASS and FNN.
#
# Workflow:
#   1. Run correct.unmixing.signatures() as usual and keep its inputs
#      (raw.data, spectra, unmixed.thresholds, scatter, spillover.spread,
#      etc.) and its output (result$spectra, result$fit.log).
#   2. prep <- sigviz.prepare(...) with the same inputs, to reproduce the
#      gated, background-subtracted raw.data and the dominance assignment
#      the correction loop actually used.
#   3. Use prep$raw.data, prep$dominant, prep$above, prep$panel,
#      prep$background.idx to drive sigviz.event.regression() and
#      sigviz.bin.regression() for a fluorophore of interest, against
#      whichever spectra you want to inspect (the starting spectra for an
#      iter-1 view, or result$spectra for the final state).
#   4. sigviz.gate.scorecard() and sigviz.bg.align.summary() work directly
#      from result$fit.log and need no recomputation.

# ---------------------------------------------------------------------------
# Local reimplementations of small internal steps (self-contained by design)
# ---------------------------------------------------------------------------

.sigviz.restricted.unmix <- function( y, spectra, active ) {

  sub <- spectra[ active, , drop = FALSE ]

  sub.plus <- tryCatch(
    solve.default( tcrossprod( sub ), sub ),
    error = function( e ) MASS::ginv( tcrossprod( sub ) ) %*% sub )

  y %*% t( sub.plus )
}


.sigviz.safe.bandwidth <- function( x ) {

  apply( x, 2, function( col ) {

    bw <- tryCatch( MASS::bandwidth.nrd( col ), error = function( e ) NA_real_ )

    if ( is.na( bw ) || !is.finite( bw ) || bw <= 0 )
      bw <- 4 * 1.06 * stats::sd( col ) * length( col ) ^ ( -1 / 5 )

    if ( is.na( bw ) || !is.finite( bw ) || bw <= 0 ) {
      col.range <- diff( range( col ) )
      bw <- if ( col.range > 0 ) col.range * 0.01 else 1
    }

    bw
  } )
}


.sigviz.main.gate <- function( scatter, gate.level = 0.1, grid.n = 128L,
                               max.events = 100000L ) {

  sc <- as.matrix( scatter[ , 1:2, drop = FALSE ] )

  fit.idx <- seq_len( nrow( sc ) )
  if ( length( fit.idx ) > max.events )
    fit.idx <- sample( fit.idx, max.events )

  bw <- .sigviz.safe.bandwidth( sc[ fit.idx, , drop = FALSE ] )
  bw <- pmax( bw, .Machine$double.eps )

  kde <- MASS::kde2d( sc[ fit.idx, 1 ], sc[ fit.idx, 2 ], h = bw, n = grid.n )

  ix   <- findInterval( sc[ , 1 ], kde$x, all.inside = TRUE )
  iy   <- findInterval( sc[ , 2 ], kde$y, all.inside = TRUE )
  dens <- kde$z[ cbind( ix, iy ) ]

  dens >= gate.level * max( kde$z )
}


.sigviz.knn.subtract <- function( raw.data, scatter, unstained, unstained.scatter,
                                  k.neighbors = 20L, max.reference = 50000L ) {

  ref.idx <- seq_len( nrow( unstained ) )
  if ( length( ref.idx ) > max.reference )
    ref.idx <- sample( ref.idx, max.reference )

  knn.idx <- FNN::knnx.index(
    data  = as.matrix( unstained.scatter[ ref.idx, , drop = FALSE ] ),
    query = as.matrix( scatter ),
    k     = k.neighbors )

  bg <- matrix( 0, nrow( raw.data ), ncol( raw.data ), dimnames = dimnames( raw.data ) )

  for ( ki in seq_len( k.neighbors ) )
    bg <- bg + as.matrix( unstained[ ref.idx, , drop = FALSE ] )[
      knn.idx[ , ki ], , drop = FALSE ]

  raw.data - bg / k.neighbors
}


# ---------------------------------------------------------------------------
# Reproduce the gate, dominance assignment and background removal
# ---------------------------------------------------------------------------
# Mirrors correct.unmixing.signatures() up to the start of its correction
# loop, so `dominant` and `above` line up with whichever `raw.data` is
# returned here.

sigviz.prepare <- function(
    raw.data, spectra, unmixed.thresholds,
    af.name           = "AF",
    scatter           = NULL,
    gate.main         = TRUE,
    gate.level        = 0.1,
    spillover.spread  = NULL,
    spread.kappa      = 2,
    bg.mode           = c( "global.mean", "scatter.knn", "none" ),
    unstained         = NULL,
    unstained.scatter = NULL,
    k.neighbors       = 20L,
    min.events        = 200L,
    background.n      = 5000L
) {

  bg.mode  <- match.arg( bg.mode )
  spectra  <- as.matrix( spectra )
  raw.data <- as.matrix( raw.data )
  panel    <- setdiff( rownames( spectra ), af.name )

  gate.keep <- NULL

  if ( gate.main && !is.null( scatter ) ) {

    gate.keep <- .sigviz.main.gate( scatter, gate.level = gate.level )
    raw.data  <- raw.data[ gate.keep, , drop = FALSE ]
    scatter   <- scatter[ gate.keep, , drop = FALSE ]

    if ( !is.null( unstained ) && !is.null( unstained.scatter ) ) {
      unst.keep         <- .sigviz.main.gate( unstained.scatter, gate.level = gate.level )
      unstained         <- unstained[ unst.keep, , drop = FALSE ]
      unstained.scatter <- unstained.scatter[ unst.keep, , drop = FALSE ]
    }
  }

  unmixed <- unmix.ols.fast( raw.data, spectra )

  if ( is.null( spillover.spread ) ) {
    above <- sweep( unmixed[ , panel, drop = FALSE ], 2,
                    unmixed.thresholds[ panel ], ">" )
  } else {
    threshold.matrix <- get.spread.thresholds(
      unmixed          = unmixed,
      thresholds       = unmixed.thresholds,
      spillover.spread = spillover.spread,
      spread.kappa     = spread.kappa,
      verbose          = FALSE )
    above <- unmixed[ , panel, drop = FALSE ] > threshold.matrix[ , panel, drop = FALSE ]
  }

  excess <- pmax( sweep( unmixed[ , panel, drop = FALSE ], 2,
                         unmixed.thresholds[ panel ], "-" ), 0 )
  dyn.range <- apply( unmixed[ , panel, drop = FALSE ], 2, stats::quantile,
                      probs = 0.999 ) - unmixed.thresholds[ panel ]
  dyn.range <- pmax( dyn.range, .Machine$double.eps )

  score    <- sweep( excess, 2, dyn.range, "/" )
  dominant <- max.col( score, ties.method = "first" )
  top      <- score[ cbind( seq_len( nrow( score ) ), dominant ) ]
  dominant[ top <= 0 ] <- 0L

  background.idx <- which( dominant == 0L )

  if ( bg.mode == "scatter.knn" ) {
    raw.data <- .sigviz.knn.subtract( raw.data, scatter, unstained, unstained.scatter,
                                      k.neighbors = k.neighbors )
  } else if ( bg.mode == "global.mean" && length( background.idx ) >= min.events ) {
    raw.data <- sweep( raw.data, 2,
                       colMeans( raw.data[ background.idx, , drop = FALSE ] ), "-" )
  }

  if ( length( background.idx ) > background.n )
    background.idx <- sample( background.idx, background.n )

  list(
    raw.data        = raw.data,
    scatter         = scatter,
    unmixed         = unmixed,
    above           = above,
    dominant        = dominant,
    panel           = panel,
    background.idx  = background.idx,
    gate.keep       = gate.keep
  )
}


# ---------------------------------------------------------------------------
# Event-level intercept-vs-slope regression: the bg.align finding
# ---------------------------------------------------------------------------
# Reproduces the exact regression correct.unmixing.signatures() uses for its
# background-confound gate (CONTEXT_signature_correction2.md, section 3.2,
# item 6): for the dominance population of one fluorophore, restrict-unmix
# against the target plus its co-active nuisance dyes, then regress the
# residual on the target's own abundance alone. The intercept and slope this
# produces are both vectors over detector channels; when they are
# anti-collinear (cosine near -1) the population is dominated by a common-mode
# background residual rather than a genuine spectral error.

sigviz.event.regression <- function(
    raw.data, spectra, dominant, above, panel, fluorophore,
    nuisance.frac = 0.5
) {

  f <- match( fluorophore, panel )
  if ( is.na( f ) )
    stop( "`fluorophore` is not in `panel`.", call. = FALSE )

  idx <- which( dominant == f )
  if ( length( idx ) < 2 )
    stop( "Too few events assigned to this fluorophore.", call. = FALSE )

  co.frac  <- colMeans( above[ idx, , drop = FALSE ] )
  nuisance <- setdiff( panel[ co.frac > nuisance.frac ], fluorophore )
  active   <- c( fluorophore, nuisance )

  y.evt <- raw.data[ idx, , drop = FALSE ]
  x.evt <- .sigviz.restricted.unmix( y.evt, spectra, active )
  r.evt <- y.evt - x.evt %*% spectra[ active, , drop = FALSE ]

  fit.evt   <- stats::lm.fit( x = cbind( 1, x.evt[ , 1 ] ), y = r.evt )
  alpha.evt <- stats::coef( fit.evt )[ 1, ]
  beta.evt  <- stats::coef( fit.evt )[ 2, ]
  alpha.evt[ !is.finite( alpha.evt ) ] <- 0
  beta.evt[ !is.finite( beta.evt ) ]   <- 0

  alpha.norm <- sqrt( sum( alpha.evt^2 ) )
  beta.norm  <- sqrt( sum( beta.evt^2 ) )
  bg.align   <- if ( alpha.norm > 0 && beta.norm > 0 )
    sum( alpha.evt * beta.evt ) / ( alpha.norm * beta.norm ) else NA_real_

  detectors <- colnames( raw.data )

  # Scale beta to the same physical units as alpha via a representative
  # abundance (the population's own median target abundance), so the two
  # traces sit on one axis instead of differing by orders of magnitude.
  scale.x <- stats::median( x.evt[ , 1 ] )

  plot.data <- data.frame(
    Detector  = factor( rep( detectors, 2 ), levels = detectors ),
    Component = rep( c( "Intercept (background)", "Slope x median abundance" ),
                     each = length( detectors ) ),
    Value     = c( alpha.evt, beta.evt * scale.x ) )

  p <- ggplot2::ggplot( plot.data,
      ggplot2::aes( x = Detector, y = Value, group = Component, colour = Component ) ) +
    ggplot2::geom_hline( yintercept = 0, linewidth = 0.3, colour = "grey60" ) +
    ggplot2::geom_path( linewidth = 1 ) +
    ggplot2::geom_point( size = 1 ) +
    ggplot2::labs(
      title    = sprintf( "%s: event-level intercept vs slope, cosine = %.3f",
                          fluorophore, bg.align ),
      subtitle = "Anti-alignment (cosine near -1) marks a background confound, not a spectral error",
      x = "Detector", y = "Regression coefficient" ) +
    ggplot2::theme_minimal() +
    ggplot2::theme( axis.text.x = ggplot2::element_text( angle = 45, hjust = 1 ),
                    legend.position = "bottom" )

  list( alpha = alpha.evt, beta = beta.evt, bg.align = bg.align,
        n.events = length( idx ), active = active, plot = p )
}


# ---------------------------------------------------------------------------
# Bin-level regression that produces the correction slope
# ---------------------------------------------------------------------------
# Reproduces the fit.slope() step inside correct.unmixing.signatures():
# abundance-binned, background-anchored, multivariate residual regression
# against the restricted design. Returns both the full detector-space
# residual by bin (a spectral ribbon, useful for spotting the curvature
# described in CONTEXT_signature_correction2.md section 2.5) and the
# residual projected onto the fitted slope direction (a 1D view of the same
# regression that produced `slope`).

sigviz.bin.regression <- function(
    raw.data, spectra, dominant, above, panel, background.idx, fluorophore,
    nuisance.frac = 0.5, n.levels = 10L, min.events = 200L
) {

  f <- match( fluorophore, panel )
  if ( is.na( f ) )
    stop( "`fluorophore` is not in `panel`.", call. = FALSE )

  idx <- which( dominant == f )
  if ( length( idx ) < min.events )
    stop( "Fewer than `min.events` events assigned to this fluorophore.",
          call. = FALSE )

  co.frac  <- colMeans( above[ idx, , drop = FALSE ] )
  nuisance <- setdiff( panel[ co.frac > nuisance.frac ], fluorophore )
  active   <- c( fluorophore, nuisance )

  y.use <- raw.data[ idx, , drop = FALSE ]
  x.use <- .sigviz.restricted.unmix( y.use, spectra, active )[ , 1 ]

  brk <- unique( stats::quantile( x.use, probs = seq( 0, 1, length.out = n.levels + 1 ) ) )
  if ( length( brk ) < 3 )
    stop( "Too few distinct abundance values to form bins.", call. = FALSE )

  bin   <- as.integer( cut( x.use, breaks = brk, include.lowest = TRUE ) )
  y.bin <- t( vapply( sort( unique( bin ) ), function( b )
    colMeans( y.use[ bin == b, , drop = FALSE ] ), numeric( ncol( y.use ) ) ) )

  has.anchor <- length( background.idx ) >= min.events
  if ( has.anchor )
    y.bin <- rbind( colMeans( raw.data[ background.idx, , drop = FALSE ] ), y.bin )

  x.bin <- .sigviz.restricted.unmix( y.bin, spectra, active )
  r.bin <- y.bin - x.bin %*% spectra[ active, , drop = FALSE ]

  fit   <- stats::lm.fit( x = cbind( 1, x.bin ), y = r.bin )
  slope <- stats::coef( fit )[ 2, ]
  slope[ !is.finite( slope ) ] <- 0

  abundance.level <- x.bin[ , 1 ]
  detectors        <- colnames( raw.data )

  ribbon.data <- data.frame(
    Bin       = rep( seq_len( nrow( r.bin ) ), each = ncol( r.bin ) ),
    Abundance = rep( abundance.level, each = ncol( r.bin ) ),
    Detector  = factor( rep( detectors, nrow( r.bin ) ), levels = detectors ),
    Residual  = as.vector( t( r.bin ) ) )

  plot.ribbon <- ggplot2::ggplot( ribbon.data,
      ggplot2::aes( x = Detector, y = Residual, group = Bin, colour = Abundance ) ) +
    ggplot2::geom_hline( yintercept = 0, linewidth = 0.3, colour = "grey60" ) +
    ggplot2::geom_path( linewidth = 0.6, alpha = 0.8 ) +
    ggplot2::scale_colour_viridis_c() +
    ggplot2::labs(
      title    = sprintf( "%s: restricted residual by abundance bin", fluorophore ),
      subtitle = "A shape that holds constant across bins supports a single-slope correction; a shape that rotates with abundance signals a variant mixture",
      x = "Detector", y = "Residual (background- and nuisance-adjusted)",
      colour = sprintf( "%s abundance", fluorophore ) ) +
    ggplot2::theme_minimal() +
    ggplot2::theme( axis.text.x = ggplot2::element_text( angle = 45, hjust = 1 ),
                    legend.position = "bottom" )

  slope.norm  <- sqrt( sum( slope^2 ) )
  unit.slope  <- if ( slope.norm > 0 ) slope / slope.norm else slope
  projection  <- as.vector( r.bin %*% unit.slope )
  fitted.proj <- as.vector( fit$fitted.values %*% unit.slope )

  is.anchor <- rep( FALSE, nrow( r.bin ) )
  if ( has.anchor ) is.anchor[ 1 ] <- TRUE

  proj.data <- data.frame(
    Abundance = abundance.level, Projection = projection, Fitted = fitted.proj,
    IsAnchor  = is.anchor )

  plot.projection <- ggplot2::ggplot( proj.data,
      ggplot2::aes( x = Abundance, y = Projection ) ) +
    ggplot2::geom_point( ggplot2::aes( shape = IsAnchor ), size = 2.2, colour = "steelblue" ) +
    ggplot2::geom_line( ggplot2::aes( y = Fitted ), colour = "firebrick", linetype = "dashed" ) +
    ggplot2::scale_shape_manual(
      values = c( `TRUE` = 8, `FALSE` = 16 ),
      labels = c( `TRUE` = "background anchor", `FALSE` = "abundance bin" ), name = NULL ) +
    ggplot2::labs(
      title    = sprintf( "%s: residual projected onto the fitted correction direction",
                          fluorophore ),
      subtitle = "Fitted line is the multivariate bin regression projected the same way; curvature away from the line signals abundance-dependent error",
      x = sprintf( "%s abundance (restricted unmix)", fluorophore ),
      y = "Residual . unit(slope)" ) +
    ggplot2::theme_minimal()

  list( x.bin = x.bin, y.bin = y.bin, r.bin = r.bin, slope = slope,
        n.bins = nrow( r.bin ), plot.ribbon = plot.ribbon,
        plot.projection = plot.projection )
}


# ---------------------------------------------------------------------------
# Gate scorecard from fit.log — no recomputation needed
# ---------------------------------------------------------------------------

sigviz.gate.scorecard <- function(
    fit.log, unmixed.thresholds,
    min.span = 5, min.explained = 0.5, min.gain = 0.002,
    max.step = 0.08, max.span.drift = 1.10, max.bg.alignment = -0.9,
    fluorophores = NULL
) {

  d <- fit.log
  if ( !is.null( fluorophores ) )
    d <- d[ d$fluorophore %in% fluorophores, , drop = FALSE ]

  d$span.threshold <- min.span * abs( unmixed.thresholds[ d$fluorophore ] )

  long <- rbind(
    data.frame( iter = d$iter, fluorophore = d$fluorophore, accepted = d$accepted,
               gate = "x.span",     value = d$x.span,     reference = d$span.threshold ),
    data.frame( iter = d$iter, fluorophore = d$fluorophore, accepted = d$accepted,
               gate = "explained",  value = d$explained,  reference = min.explained ),
    data.frame( iter = d$iter, fluorophore = d$fluorophore, accepted = d$accepted,
               gate = "gain",       value = d$gain,        reference = min.gain ),
    data.frame( iter = d$iter, fluorophore = d$fluorophore, accepted = d$accepted,
               gate = "rel.step",   value = d$rel.step,    reference = max.step ),
    data.frame( iter = d$iter, fluorophore = d$fluorophore, accepted = d$accepted,
               gate = "span.drift", value = d$span.drift,  reference = max.span.drift ),
    data.frame( iter = d$iter, fluorophore = d$fluorophore, accepted = d$accepted,
               gate = "bg.align",   value = d$bg.align,    reference = max.bg.alignment )
  )

  ggplot2::ggplot( long,
      ggplot2::aes( x = iter, y = value, colour = accepted, group = fluorophore ) ) +
    ggplot2::geom_hline( ggplot2::aes( yintercept = reference ), linetype = "dashed",
                         colour = "grey40", na.rm = TRUE ) +
    ggplot2::geom_line( colour = "grey80", linewidth = 0.3 ) +
    ggplot2::geom_point( size = 1.8, na.rm = TRUE ) +
    ggplot2::facet_grid( gate ~ fluorophore, scales = "free_y" ) +
    ggplot2::scale_colour_manual( values = c( `TRUE` = "steelblue", `FALSE` = "firebrick" ) ) +
    ggplot2::labs(
      title    = "Signature correction gate scorecard",
      subtitle = "Dashed line is the acceptance bound for that gate, per iteration",
      x = "Iteration", y = NULL, colour = "Accepted" ) +
    ggplot2::theme_minimal() +
    ggplot2::theme( strip.text.x = ggplot2::element_text( angle = 90 ),
                    axis.text.x  = ggplot2::element_text( angle = 0 ) )
}


sigviz.bg.align.summary <- function( fit.log, max.bg.alignment = -0.9 ) {

  ggplot2::ggplot( fit.log,
      ggplot2::aes( x = stats::reorder( fluorophore, bg.align, FUN = stats::median ),
                   y = bg.align, colour = accepted ) ) +
    ggplot2::geom_hline( yintercept = max.bg.alignment, linetype = "dashed", colour = "grey40" ) +
    ggplot2::geom_jitter( width = 0.15, height = 0, size = 1.8, na.rm = TRUE ) +
    ggplot2::coord_flip() +
    ggplot2::scale_colour_manual( values = c( `TRUE` = "steelblue", `FALSE` = "firebrick" ) ) +
    ggplot2::labs(
      title    = "Background-confound cosine across iterations",
      subtitle = sprintf(
        "Points at or below %.2f are rejected as background confounds, not spectral errors",
        max.bg.alignment ),
      x = NULL, y = "cosine( event-level intercept, event-level slope )", colour = "Accepted" ) +
    ggplot2::theme_minimal()
}


# ---------------------------------------------------------------------------
# Before/after row comparison
# ---------------------------------------------------------------------------

sigviz.row.overlay <- function( spectra.old, spectra.new, fluorophore ) {

  if ( !fluorophore %in% rownames( spectra.old ) ||
       !fluorophore %in% rownames( spectra.new ) )
    stop( "`fluorophore` must be a row of both spectra matrices.", call. = FALSE )

  detectors <- colnames( spectra.old )

  plot.data <- data.frame(
    Detector  = factor( rep( detectors, 2 ), levels = detectors ),
    Version   = rep( c( "Before", "After" ), each = length( detectors ) ),
    Intensity = c( spectra.old[ fluorophore, ], spectra.new[ fluorophore, ] ) )

  cs <- sum( spectra.old[ fluorophore, ] * spectra.new[ fluorophore, ] ) /
    ( sqrt( sum( spectra.old[ fluorophore, ]^2 ) ) *
        sqrt( sum( spectra.new[ fluorophore, ]^2 ) ) )
  deg <- 180 / pi * acos( pmin( 1, pmax( -1, cs ) ) )

  ggplot2::ggplot( plot.data,
      ggplot2::aes( x = Detector, y = Intensity, group = Version, colour = Version ) ) +
    ggplot2::geom_path( linewidth = 1 ) +
    ggplot2::geom_point( size = 1 ) +
    ggplot2::labs(
      title = sprintf( "%s: spectrum before and after correction (%.2f deg)", fluorophore, deg ),
      x = "Detector", y = "Normalized Intensity" ) +
    ggplot2::theme_minimal() +
    ggplot2::theme( axis.text.x = ggplot2::element_text( angle = 45, hjust = 1 ),
                    legend.position = "bottom" )
}


# ---------------------------------------------------------------------------
# Example usage (not executed)
# ---------------------------------------------------------------------------
#
# result <- correct.unmixing.signatures(
#   raw.data = raw.data, spectra = spectra, unmixed.thresholds = thresholds,
#   scatter = scatter, spillover.spread = variants$spillover.spread,
#   bg.mode = "scatter.knn", unstained = unstained,
#   unstained.scatter = unstained.scatter, true.spectra = true.spectra )
#
# prep <- sigviz.prepare(
#   raw.data = raw.data, spectra = spectra, unmixed.thresholds = thresholds,
#   scatter = scatter, spillover.spread = variants$spillover.spread,
#   bg.mode = "scatter.knn", unstained = unstained,
#   unstained.scatter = unstained.scatter )
#
# er <- sigviz.event.regression(
#   prep$raw.data, result$spectra, prep$dominant, prep$above, prep$panel,
#   fluorophore = "BV711" )
# er$plot
#
# br <- sigviz.bin.regression(
#   prep$raw.data, result$spectra, prep$dominant, prep$above, prep$panel,
#   prep$background.idx, fluorophore = "BV711" )
# br$plot.ribbon
# br$plot.projection
#
# sigviz.gate.scorecard( result$fit.log, thresholds,
#                        fluorophores = c( "BV711", "PerCP-eFluor 710", "PE" ) )
# sigviz.bg.align.summary( result$fit.log )
# sigviz.row.overlay( spectra, result$spectra, "BV711" )
