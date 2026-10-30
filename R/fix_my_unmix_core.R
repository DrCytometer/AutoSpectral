# fix_my_unmix_core.R

#' Phase one and phase two of `fix.my.unmix()` on events already projected
#' onto the panel.
#'
#' Takes the unstained and fully stained events unmixed with `spectra.fluor`
#' (`.fix.projector()`), so the same core serves a single pass and every round
#' of `.fix.smooth()`. The arguments after `seed` are `fix.my.unmix()`'s, passed
#' through unchanged.
#'
#' @param spectra.fluor Numeric matrix (fluorophores x detectors), the spectra
#'   the events were unmixed with and the run starts from.
#' @param unstained.fit,stained.fit Lists with `abundance` (events x
#'   fluorophores) and `residual` (events x detectors).
#' @param background.basis Numeric matrix of the design's background rows, or
#'   `NULL`.
#' @param af.gate.spectra Numeric matrix the AF likeness of a change is
#'   measured against, or `NULL`.
#' @param freeze.coupled Logical, whether fluorophores inseparable from the
#'   background basis (`max.hotspot`) are frozen.
#' @param spillover.spread Spillover spread matrix, or `NULL`.
#' @param variants `get.spectral.variants()` output, or `NULL`.
#' @param reference.spectra Numeric matrix, the rows the plausibility box is
#'   measured from.
#' @param tandems Named character vector, tandem to parent.
#' @param seed Integer, the base of the leakage gate's seeds.
#'
#' @return A list of the quantities `fix.my.unmix()` returns, with `spectra`
#'   over the fluorophores only, and `accepted`, per fluorophore.
#'
#' @noRd
.fix.my.unmix.core <- function(
    spectra.fluor, unstained.fit, stained.fit, background.basis,
    af.gate.spectra, freeze.coupled, spillover.spread, variants,
    reference.spectra, tandems, seed,
    max.iter, downsample, downsample.background.frac, downsample.min.stratum,
    unstained.threshold, lower.unstained.threshold, unstained.margin,
    spread.kappa, negative.tail, source.events, envelope.quantiles,
    disagreement.reference, envelope.scale.quantile, min.negative.frac,
    max.disagreement, min.hypernegative.frac, max.hypernegative.frac,
    min.hypernegative.events, hypernegative.dominant, leakage.prior,
    span.fraction, min.negative.events, min.bin.negative, min.span, min.rise,
    n.levels, min.bin.events, n.levels.pair, multivariate, ridge,
    max.truncated.events, max.mask.passes, estimator, glasso.check,
    glasso.warm, glasso.n.lambda, glasso.lambda.min.ratio, glasso.margin.frac,
    glasso.tolerance, af.gate, af.gate.max, plausibility.ratio,
    plausibility.floor, plausibility.reject, plausibility.reset,
    phase.two.rows, phase.two.active, phase.two.present, phase.two.fallback,
    phase.two.rescue, rescue.reasons, rescue.min.weight, leakage.channels,
    tandem.candidate, tandem.min.coefficient, source.dominant, spread.addback,
    anchor.weight, max.coefficient, convergence.threshold, convergence.quantile,
    stop.rule, stop.trust, stop.confirm, coefficient.decay, update.spectra,
    step, intercept, min.explained, max.explained, max.resid, max.intercept,
    min.bg.align, max.angle, max.clamp.frac, max.anchor, max.vif,
    max.condition.increase, peak.shift.min.rel, max.hotspot, leakage.margin,
    n.threads, keep.history, verbose ) {

  eps <- .Machine$double.eps

  fluorophores  <- rownames( spectra.fluor )
  fluorophore.n <- length( fluorophores )

  identity.matrix <- diag( fluorophore.n )
  dimnames( identity.matrix ) <- list( fluorophores, fluorophores )
  off.mask <- row( identity.matrix ) != col( identity.matrix )

  design <- if ( is.null( background.basis ) ) spectra.fluor else
    rbind( background.basis, spectra.fluor )

  if ( nrow( design ) > ncol( design ) )
    stop( "The background basis plus the panel exceeds the detector count.",
          call. = FALSE )

  # The plausibility box is measured from the rows the whole correction started
  # from, so the rounds of a smoothed run cannot walk a row out of it.
  reference.spectra <- reference.spectra[ fluorophores, , drop = FALSE ]
  reference.spectra <- reference.spectra /
    pmax( apply( reference.spectra, 1, max ), eps )

  plausible <- plausibility.ratio > 0
  allow <- if ( plausible ) .fix.plausible.change(
    reference.spectra, plausibility.ratio, plausibility.floor ) else NULL

  af.hotspot <- calculate.hotspot.matrix( design )

  af.frozen <- character( 0 )
  if ( freeze.coupled && !is.null( background.basis ) ) {
    coupling  <- apply(
      af.hotspot[ fluorophores, rownames( background.basis ), drop = FALSE ],
      1, max )
    af.frozen <- names( coupling )[ coupling > max.hotspot ]
  }

  if ( verbose && length( af.frozen ) > 0 )
    message( sprintf(
      "\033[31mFrozen, inseparable from the background basis: %s.\033[0m",
      paste( af.frozen, collapse = ", " ) ) )

  if ( af.gate && is.null( af.gate.spectra ) )
    stop( "`af.gate` needs `af.spectra` or a background basis.", call. = FALSE )

  prior <- NULL

  if ( leakage.prior && !is.null( variants ) )
    prior <- get.variant.leakage.prior(
      spectra       = spectra.fluor,
      variants      = variants,
      extra.rows    = background.basis,
      span.fraction = span.fraction,
      verbose       = verbose )

  n.events     <- nrow( stained.fit$abundance )
  stained.full <- stained.fit

  # A uniform downsample would thin a dim or rare fluorophore's own positive
  # population at the same rate as the background bulk, so every fluorophore's
  # stratum gets its own quota, with a floor that keeps a small one whole.
  if ( !isFALSE( downsample ) && !is.null( downsample ) &&
       as.integer( downsample ) < n.events ) {

    thr.rough <- unstained.margin * apply(
      unstained.fit$abundance, 2, stats::quantile,
      probs = unstained.threshold, names = FALSE )
    names( thr.rough ) <- fluorophores

    keep.idx <- .fix.stratified.sample(
      abundance        = stained.fit$abundance,
      thresholds       = thr.rough,
      threshold.matrix = get.spread.thresholds(
        unmixed          = stained.fit$abundance,
        thresholds       = thr.rough,
        spillover.spread = spillover.spread,
        spread.kappa     = spread.kappa,
        verbose          = FALSE ),
      n.total          = as.integer( downsample ),
      background.frac  = downsample.background.frac,
      min.stratum      = downsample.min.stratum )

    stained.fit <- lapply( stained.fit, function( v )
      if ( is.matrix( v ) ) v[ keep.idx, , drop = FALSE ] else v[ keep.idx ] )
  }

  boundaries <- function( unstained.comp ) {
    upper <- unstained.margin * apply(
      unstained.comp, 2, stats::quantile, probs = unstained.threshold,
      names = FALSE )
    lower <- unstained.margin * apply(
      unstained.comp, 2, stats::quantile, probs = 1 - lower.unstained.threshold,
      names = FALSE )
    names( upper ) <- fluorophores
    names( lower ) <- fluorophores
    list( upper = upper, lower = lower )
  }

  # ---------------------------------------------------------------------------
  # Phase one: residual spillover matrix
  # ---------------------------------------------------------------------------

  spillover.curr <- identity.matrix
  trust          <- identity.matrix
  pair.log       <- NULL
  source.log     <- NULL

  # Warm-starts each pair's robust fit from the slope the previous iteration
  # left behind.
  slope.warm <- identity.matrix * 0

  convergence.log <- data.frame(
    iter = integer(), delta = numeric(), delta.quantile = numeric(),
    delta.max = numeric(), n.fitted = integer(), stringsAsFactors = FALSE )

  delta.history     <- rep( NA_real_, 3L )
  below.run         <- 0L
  stop.after        <- if ( stop.rule == "settled" ) as.integer( stop.confirm ) else 1L
  stop.reason       <- "iteration limit"
  spillover.history <- if ( keep.history ) list() else NULL

  for ( iter in seq_len( as.integer( max.iter ) ) ) {

    compensation.curr <- tryCatch( solve( spillover.curr ),
                                   error = function( e ) NULL )

    if ( is.null( compensation.curr ) ) {
      stop.reason <- "singular spillover"
      warning( "Spillover matrix became singular; stopping iteration.",
               call. = FALSE )
      break
    }

    unmixed.comp   <- stained.fit$abundance   %*% compensation.curr
    unstained.comp <- unstained.fit$abundance %*% compensation.curr

    # Recomputed in the current compensated space each iteration, so the
    # boundaries follow the data as the matrix is refined. The lower boundary
    # is the unstained population's own lower tail rather than a mirror of the
    # upper one: AF stretches the positive tail without stretching the negative.
    bounds         <- boundaries( unstained.comp )
    thresholds     <- bounds$upper
    neg.thresholds <- bounds$lower

    neg.var <- apply( unstained.comp, 2, stats::mad )^2
    names( neg.var ) <- fluorophores

    threshold.matrix <- get.spread.thresholds(
      unmixed          = unmixed.comp,
      thresholds       = thresholds,
      spillover.spread = spillover.spread,
      spread.kappa     = spread.kappa,
      verbose          = FALSE )

    neg.threshold.matrix <- get.spread.thresholds(
      unmixed          = unmixed.comp,
      thresholds       = neg.thresholds,
      spillover.spread = spillover.spread,
      spread.kappa     = spread.kappa,
      side             = "lower",
      verbose          = FALSE )

    evidence <- .fix.source.evidence(
      unmixed.comp     = unmixed.comp,
      unstained.comp   = unstained.comp,
      thresholds       = thresholds,
      spillover.spread = spillover.spread,
      spread.kappa     = spread.kappa,
      negative.tail    = negative.tail,
      source.events    = source.events,
      source.dominant  = source.dominant,
      n.events         = n.events )

    t.clear  <- evidence$t.clear
    dominant <- evidence$dominant

    lasso <- function( warm ) .fix.glasso.matrix(
      unmixed.comp         = unmixed.comp,
      threshold.matrix     = threshold.matrix,
      thresholds           = thresholds,
      neg.var              = neg.var,
      warm                 = warm,
      min.negative.events  = min.negative.events,
      max.mask.passes      = max.mask.passes,
      max.truncated.events = max.truncated.events,
      max.coefficient      = max.coefficient,
      min.span             = min.span,
      min.rise             = min.rise,
      n.lambda             = glasso.n.lambda,
      lambda.min.ratio     = glasso.lambda.min.ratio,
      margin.frac          = glasso.margin.frac,
      tolerance            = glasso.tolerance )

    joint <- if ( glasso.check && !glasso.warm ) lasso( NULL ) else NULL

    marker.spillover <- identity.matrix
    trust            <- identity.matrix
    anchored         <- identity.matrix * 0

    pair.log   <- list()
    source.log <- list()
    candidates <- list()

    # One batched call per source: everything that depends only on the source
    # (the dominance mask, the positivity split, the bins) is computed once and
    # shared across every channel it might spill into.
    for ( source in fluorophores ) {

      channel.set <- setdiff( fluorophores, source )

      spread.var.vec <- rep( 0, length( channel.set ) )
      if ( !is.null( spillover.spread ) &&
           source %in% rownames( spillover.spread ) ) {
        avail <- intersect( channel.set, colnames( spillover.spread ) )
        spread.var.vec[ match( avail, channel.set ) ] <-
          spillover.spread[ source, avail ]
        spread.var.vec[ !is.finite( spread.var.vec ) ] <- 0
        spread.var.vec <- pmax( spread.var.vec, 0 )
      }

      w.source <- unname( evidence$weight[ source ] )

      source.log[[ length( source.log ) + 1L ]] <- data.frame(
        source            = source,
        negative.bound    = unname( evidence$bound[ source ] ),
        clear.events      = unname( evidence$n.clear[ source ] ),
        baseline.events   = unname( evidence$baseline[ source ] ),
        excess.events     = unname( evidence$excess[ source ] ),
        source.weight     = w.source,
        af.likeness       = NA_real_,
        af.gated          = FALSE,
        plausible.project = 0,
        row.names = NULL, stringsAsFactors = FALSE )

      est.batch <- if ( w.source <= 0 ) NULL else .fix.envelope.slope.batch(
        x.source                = unmixed.comp[ , source ],
        X.target                = unmixed.comp[ , channel.set, drop = FALSE ],
        threshold.source        = t.clear[ , source ],
        Threshold.target        = threshold.matrix[ , channel.set, drop = FALSE ],
        spread.var              = spread.var.vec,
        neg.var                 = neg.var[ channel.set ],
        source.mask             = if ( is.null( dominant ) ) NULL else
          dominant == match( source, fluorophores ),
        quantiles               = envelope.quantiles,
        n.levels                = n.levels.pair,
        min.events              = min.negative.events,
        min.bin.negative        = min.bin.negative,
        spread.addback          = spread.addback,
        anchor.weight           = anchor.weight,
        max.truncated.events    = max.truncated.events,
        max.coefficient         = max.coefficient,
        max.mask.passes         = max.mask.passes,
        start.slope             = slope.warm[ source, channel.set ],
        estimator               = if ( estimator == "truncated" ) "truncated" else "anchored",
        Threshold.lower         = if ( estimator == "truncated" ) NULL else
          neg.threshold.matrix[ , channel.set, drop = FALSE ],
        disagreement.reference  = disagreement.reference,
        envelope.scale.quantile = envelope.scale.quantile,
        n.threads               = n.threads )

      for ( ci in seq_along( channel.set ) ) {

        channel <- channel.set[ ci ]
        est <- if ( is.null( est.batch ) ) NULL else
          as.list( est.batch[ ci, , drop = FALSE ] )

        slope.fit <- if ( is.null( est ) ) NA_real_ else
          if ( estimator == "truncated" ) est$slope.truncated else est$slope.anchored

        if ( is.finite( slope.fit ) ) {
          slope.warm[ source, channel ] <- slope.fit
          if ( abs( slope.fit ) <= max.coefficient )
            anchored[ source, channel ] <- slope.fit
        }

        w <- 0

        slope.use <- if ( is.null( est ) ) NA_real_ else
          if ( estimator == "envelope" ) est$slope else slope.fit

        glasso.agreement <- NA_real_

        if ( !is.null( est ) && is.finite( slope.use ) &&
             est$coverage >= min.negative.frac &&
             abs( slope.use ) <= max.coefficient &&
             abs( slope.use ) * est$span >= min.rise * est$noise &&
             est$span > min.span * abs( thresholds[ source ] ) ) {

          prior.var <- if ( is.null( prior ) ) slope.use^2 else
            prior$variance[ source, channel ]

          if ( !is.finite( prior.var ) || prior.var <= 0 )
            prior.var <- slope.use^2

          # The envelope's standard error, else the anchored fit's own. A
          # coefficient whose uncertainty cannot be judged earns no trust.
          se.use <- est$se
          if ( !( is.finite( se.use ) && se.use > 0 ) ) se.use <- est$se.fit

          w <- if ( is.finite( se.use ) && se.use > 0 )
            prior.var / ( prior.var + se.use^2 ) else 0

          w <- w * w.source

          if ( is.finite( est$disagreement ) &&
               est$disagreement > max.disagreement )
            w <- w * max.disagreement / est$disagreement

          if ( !is.null( joint ) ) {
            g <- joint[ source, channel ]
            glasso.agreement <- if ( g == 0 ) 0 else
              max( 0, 1 - abs( slope.use - g ) / ( abs( slope.use ) + abs( g ) ) )
            w <- w * glasso.agreement
          }

          marker.spillover[ source, channel ] <- slope.use
          candidates[[ length( candidates ) + 1L ]] <- list(
            source = source, channel = channel, row = length( pair.log ) + 1L )
        }

        trust[ source, channel ] <- w

        pair.log[[ length( pair.log ) + 1L ]] <- data.frame(
          source              = source,
          channel             = channel,
          slope               = if ( is.null( est ) ) NA_real_ else est$slope,
          slope.truncated     = if ( is.null( est ) ) NA_real_ else est$slope.truncated,
          slope.anchored      = if ( is.null( est ) ) NA_real_ else est$slope.anchored,
          slope.median        = if ( is.null( est ) ) NA_real_ else est$slope.alt,
          slope.glasso        = if ( is.null( joint ) ) NA_real_ else joint[ source, channel ],
          se                  = if ( is.null( est ) ) NA_real_ else est$se,
          se.fit              = if ( is.null( est ) ) NA_real_ else est$se.fit,
          glasso.agreement    = glasso.agreement,
          source.weight       = w.source,
          disagreement        = if ( is.null( est ) ) NA_real_ else est$disagreement,
          coverage            = if ( is.null( est ) ) NA_real_ else est$coverage,
          span                = if ( is.null( est ) ) NA_real_ else est$span,
          hypernegative.base  = NA_real_,
          hypernegative.after = NA_real_,
          hypernegative.new   = NA_real_,
          trust               = w,
          row.names           = NULL,
          stringsAsFactors    = FALSE )
      }
    }

    # The lasso fitted after the pair slopes, with them removed from its first
    # negativity mask, so a source whose positive residual spillover lifts its
    # events above the channel's boundary is not missing from its fit.
    if ( glasso.check && glasso.warm && length( candidates ) > 0L ) {

      joint <- lasso( anchored )

      for ( cand in candidates ) {
        s <- marker.spillover[ cand$source, cand$channel ]
        g <- joint[ cand$source, cand$channel ]
        agreement <- if ( g == 0 ) 0 else
          max( 0, 1 - abs( s - g ) / ( abs( s ) + abs( g ) ) )
        trust[ cand$source, cand$channel ] <- trust[ cand$source, cand$channel ] * agreement
        pair.log[[ cand$row ]]$slope.glasso     <- g
        pair.log[[ cand$row ]]$glasso.agreement <- agreement
        pair.log[[ cand$row ]]$trust            <- trust[ cand$source, cand$channel ]
      }
    }

    # A negative correlation between a source's abundance and a channel's
    # corrected value has no biological explanation, so trust falls linearly
    # from `min.hypernegative.frac` to zero at `max.hypernegative.frac` of the
    # source's clearly positive population left below the channel's negative
    # boundary. With `hypernegative.dominant` only the events the source
    # dominates are judged: an event another fluorophore dominates reads
    # positive for the source only through that fluorophore's error.
    frac.range <- max( max.hypernegative.frac - min.hypernegative.frac, eps )

    for ( cand in candidates ) {

      src <- cand$source
      ch  <- cand$channel

      source.positive <- unmixed.comp[ , src ] > t.clear[ , src ]
      if ( hypernegative.dominant && !is.null( dominant ) )
        source.positive <- source.positive & dominant == match( src, fluorophores )

      if ( sum( source.positive ) < min.hypernegative.events ) next

      boundary <- neg.threshold.matrix[ source.positive, ch ]
      values   <- unmixed.comp[ source.positive, ch ]

      hypernegative.base  <- mean( values < boundary )
      hypernegative.after <- mean(
        values - marker.spillover[ src, ch ] * unmixed.comp[ source.positive, src ] <
          boundary )

      trust[ src, ch ] <- trust[ src, ch ] *
        ( 1 - min( max( hypernegative.after - min.hypernegative.frac, 0 ) / frac.range, 1 ) )

      pair.log[[ cand$row ]]$hypernegative.base  <- hypernegative.base
      pair.log[[ cand$row ]]$hypernegative.after <- hypernegative.after
      pair.log[[ cand$row ]]$hypernegative.new   <- max( hypernegative.after - hypernegative.base, 0 )
      pair.log[[ cand$row ]]$trust               <- trust[ src, ch ]
    }

    slope.error <- marker.spillover - identity.matrix

    source.log <- do.call( rbind, source.log )

    # The change each source's back-solved row would take this step, and how
    # much it looks like autofluorescence. With `af.gate`, a source whose step
    # looks like AF takes none: marker-associated AF would otherwise be booked
    # as spillover.
    row.change <- ( ( trust * slope.error ) * off.mask ) %*%
      ( spillover.curr %*% spectra.fluor )

    gated <- character( 0 )

    for ( f in seq_len( fluorophore.n ) ) {

      like <- if ( is.null( af.gate.spectra ) ) NA_real_ else
        .fix.af.likeness( row.change[ f, ], af.gate.spectra )

      source.log$af.likeness[ f ] <- like

      if ( af.gate && is.finite( like ) && like > af.gate.max ) {
        trust[ f, off.mask[ f, ] ] <- 0
        source.log$af.gated[ f ] <- TRUE
        gated <- c( gated, fluorophores[ f ] )
      }
    }

    if ( length( gated ) > 0L )
      pair.log <- lapply( pair.log, function( r ) {
        if ( r$source %in% gated ) r$trust <- 0
        r
      } )

    # The next anchored fit starts from the slope this step leaves behind.
    if ( estimator == "anchored" ) slope.warm <- slope.warm * ( 1 - trust )

    # Under "settled", a coefficient below `stop.trust` does not hold phase one
    # open: a slope too weakly trusted ever to finish would otherwise floor the
    # stop test at its own size.
    judged.error <- if ( stop.rule == "settled" )
      slope.error * ( trust >= stop.trust ) else slope.error

    off.diagonal <- abs( judged.error[ off.mask ] )

    delta          <- stats::sd( judged.error )
    delta.quantile <- unname( stats::quantile(
      off.diagonal, probs = convergence.quantile, names = FALSE ) )
    delta.max      <- max( off.diagonal )
    n.fitted       <- sum( trust > 0 ) - fluorophore.n

    convergence.log[ nrow( convergence.log ) + 1L, ] <-
      list( iter, delta, delta.quantile, delta.max, n.fitted )

    if ( verbose )
      message( sprintf(
        "\033[34miter %d: delta %.5f, q%.0f %.5f, max %.5f, %d coefficient(s) fitted\033[0m",
        iter, delta, 100 * convergence.quantile, delta.quantile, delta.max,
        n.fitted ) )

    # Trust damps the step; it must not touch the estimate or the stop test, or
    # the fixed point becomes trust-weighted identity rather than identity.
    spillover.next <- spillover.curr + ( trust * slope.error ) %*% spillover.curr

    decay.factor <- 1 - coefficient.decay * ( 1 - trust )
    diag( decay.factor ) <- 1
    spillover.next <- spillover.next * decay.factor

    if ( keep.history )
      spillover.history[[ iter ]] <- list(
        marker.spillover = marker.spillover,
        trust            = trust,
        spillover.next   = spillover.next,
        coefficient.log  = do.call( rbind, pair.log ),
        source.log       = source.log )

    diag.next <- diag( spillover.next )

    if ( any( !is.finite( diag.next ) ) || any( diag.next <= 0 ) ) {
      stop.reason <- "diagonal left the positive range"
      warning( "Spillover diagonal left the positive range; stopping iteration.",
               call. = FALSE )
      break
    }

    spillover.next <- sweep( spillover.next, 1, diag.next, "/" )

    # The part of each back-solved row outside its plausibility box is taken
    # back out of its coefficients, so a coefficient carrying co-expression or
    # AF stops steering phase two's regressors. The box bounds rows, not single
    # coefficients: bounding those holds back real corrections.
    if ( plausible ) {
      projected <- .fix.plausible.project( spillover.next, spectra.fluor,
                                           reference.spectra, allow )
      spillover.next <- projected$spillover
      source.log$plausible.project <- unname( projected$share[ fluorophores ] )
    }

    spillover.curr <- spillover.next

    below.run <- if ( delta.quantile < convergence.threshold ) below.run + 1L else 0L

    if ( below.run >= stop.after ) {
      stop.reason <- "converged"
      if ( verbose ) message( "\033[32mConverged.\033[0m" )
      break
    }

    delta.history <- c( delta.history[ -1L ], delta )

    if ( all( is.finite( delta.history ) ) &&
         mean( diff( delta.history ) ) >= 0 ) {
      stop.reason <- "stalled"
      if ( verbose ) message( "\033[33mRefinement stalled, stopping.\033[0m" )
      break
    }

    if ( iter == max.iter && verbose )
      message( "\033[33mReached iteration limit.\033[0m" )
  }

  coefficient.log <- do.call( rbind, pair.log )

  compensation.curr <- solve( spillover.curr )

  # Reported over every event of the run, not only a downsample phase one
  # iterated on.
  unmixed.final  <- stained.full$abundance %*% compensation.curr
  residual.final <- stained.full$residual

  bounds.final         <- boundaries( unstained.fit$abundance %*% compensation.curr )
  thresholds.final     <- bounds.final$upper
  neg.thresholds.final <- bounds.final$lower

  threshold.matrix.final <- get.spread.thresholds(
    unmixed          = unmixed.final,
    thresholds       = thresholds.final,
    spillover.spread = spillover.spread,
    spread.kappa     = spread.kappa,
    verbose          = FALSE )

  neg.threshold.matrix.final <- get.spread.thresholds(
    unmixed          = unmixed.final,
    thresholds       = neg.thresholds.final,
    spillover.spread = spillover.spread,
    spread.kappa     = spread.kappa,
    side             = "lower",
    verbose          = FALSE )

  dominant.final <- .fix.dominance(
    unmixed.final, threshold.matrix.final,
    pmax( apply( unmixed.final, 2, stats::quantile, probs = 0.999,
                 names = FALSE ) - thresholds.final, eps ) )

  spectra.backsolved <- spillover.curr %*% spectra.fluor
  row.max <- apply( spectra.backsolved, 1, max )
  row.max[ !is.finite( row.max ) | row.max <= 0 ] <- 1
  spectra.backsolved <- spectra.backsolved / row.max

  # ---------------------------------------------------------------------------
  # Phase two: signatures re-measured in raw space
  # ---------------------------------------------------------------------------

  spectra.new   <- spectra.fluor
  signature.log <- NULL

  if ( update.spectra ) {

    unmixed.comp <- stained.fit$abundance %*% compensation.curr
    residual     <- stained.fit$residual
    thresholds   <- thresholds.final

    threshold.matrix <- get.spread.thresholds(
      unmixed          = unmixed.comp,
      thresholds       = thresholds,
      spillover.spread = spillover.spread,
      spread.kappa     = spread.kappa,
      verbose          = FALSE )

    dominant <- .fix.dominance(
      unmixed.comp, threshold.matrix,
      pmax( apply( unmixed.comp, 2, stats::quantile, probs = 0.999,
                   names = FALSE ) - thresholds, eps ) )

    background.idx <- which( dominant == 0L )

    condition.curr <- calculate.condition.number( design )

    # A row phase one left alone has no in-span change to re-measure.
    flagged <- apply( abs( spillover.curr - identity.matrix ), 1, max ) >=
      convergence.threshold
    names( flagged ) <- fluorophores

    source.weight <- if ( is.null( source.log ) )
      stats::setNames( rep( 0, fluorophore.n ), fluorophores ) else
        stats::setNames( source.log$source.weight, source.log$source )

    log.rows <- list()

    for ( f in seq_along( fluorophores ) ) {

      j   <- fluorophores[ f ]
      idx <- which( dominant == f )

      reject <- NA_character_

      if ( j %in% af.frozen ) reject <- "af.coupled"

      if ( is.na( reject ) && phase.two.rows == "flagged" && !flagged[ j ] )
        reject <- "not.flagged"

      candidate <- NULL

      if ( is.na( reject ) ) {

        # Only the dyes present in the row's own population are regressors of
        # its signature, so a population of a few hundred events still gets a
        # joint fit instead of the fallback that books every other row's error
        # against it.
        active <- fluorophores

        if ( phase.two.active == "present" && length( idx ) > 0L ) {

          u.pop  <- unmixed.comp[ idx, , drop = FALSE ]
          tm.pop <- threshold.matrix[ idx, , drop = FALSE ]

          frac.above <- colMeans( u.pop > tm.pop )
          corr <- suppressWarnings( stats::cor( u.pop[ , j ], u.pop ) )[ 1, ]
          corr[ !is.finite( corr ) ] <- 0
          mean.abs <- colMeans( abs( u.pop ) ) / pmax( colMeans( tm.pop ), eps )

          present <- frac.above >= phase.two.present[[ "frac" ]] |
            abs( corr ) >= phase.two.present[[ "corr" ]] |
            mean.abs >= phase.two.present[[ "mean" ]]
          present[ j ] <- TRUE

          active <- fluorophores[ present ]
        }

        candidate <- extract.raw.signature(
          raw.data       = residual[ idx, , drop = FALSE ],
          spectra        = spectra.new,
          abundance      = unmixed.comp[ idx, , drop = FALSE ],
          target         = j,
          active         = active,
          intercept      = intercept,
          multivariate   = multivariate,
          ridge          = ridge,
          n.levels       = n.levels,
          min.bin.events = min.bin.events,
          min.events     = min.negative.events,
          background.raw = if ( length( background.idx ) > 0 )
            residual[ background.idx, , drop = FALSE ] else NULL )

        if ( is.null( candidate ) ) reject <- "no.fit"
      }

      if ( is.na( reject ) ) {

        st <- candidate$stats

        # Ordered cheapest and most decisive first; the first to fire names the
        # reason. Under the joint fit the fluorophore's own term is a partial
        # slope, so `explained.total` is the one that must stay near one.
        reject <-
          if ( st$x.span <= min.span * abs( thresholds[ j ] ) ) "span" else
          if ( !is.finite( st$explained.total ) ||
               st$explained.total < min.explained ||
               st$explained.total > max.explained ) "explained" else
          if ( st$resid.rel > max.resid ) "fit" else
          if ( st$intercept.rel > max.intercept ) "offset" else
          if ( is.finite( st$bg.align ) && st$bg.align < min.bg.align )
            "background.confound" else
          if ( st$clamp.frac > max.clamp.frac ) "negative.mass" else
          if ( st$deg.change > max.angle ) "angle" else
          if ( is.finite( st$anchor.rel ) && st$anchor.rel > max.anchor )
            "background" else
          if ( st$vif.target > max.vif ) "collinear" else
          if ( st$peak.new != st$peak.curr &&
               st$peak.new.rel < peak.shift.min.rel ) "peak.shift" else
            NA_character_
      }

      # A population too small for the joint fit was re-measured with every
      # other row held at its spectrum, which books their errors against it.
      not.joint <- !is.null( candidate ) && multivariate && fluorophore.n > 1L &&
        isFALSE( candidate$stats$joint )

      if ( is.na( reject ) && not.joint && phase.two.fallback == "reject" )
        reject <- "not.joint"

      # The leakage gate sums the off-target medians of the row's population;
      # with "negative" only over the channels that population is negative for,
      # so a marker expressed on the same cells is not counted as leakage.
      leak.panel <- fluorophores

      if ( leakage.channels == "negative" && length( idx ) > 0L ) {
        share <- colMeans( unmixed.comp[ idx, , drop = FALSE ] >
                             threshold.matrix[ idx, , drop = FALSE ] )
        leak.panel <- fluorophores[ fluorophores == j | share < 0.5 ]
      }

      # The row gates on one candidate row: the plausibility box and the AF gate
      # (`box`), conditioning and leakage.
      judge <- function( target.row, box = TRUE ) {

        out <- list( reject = NA_character_, trial = NULL,
                     af.likeness = NA_real_, implausible = NA_real_,
                     condition.after = NA_real_, leak.before = NA_real_,
                     leak.after = NA_real_ )

        stepped <- pmax( spectra.new[ j, ] + step * ( target.row - spectra.new[ j, ] ), 0 )

        if ( max( stepped ) <= 0 ) {
          out$reject <- "no.fit"
          return( out )
        }

        trial <- spectra.new
        trial[ j, ] <- stepped / max( stepped )

        if ( plausible && box ) {
          clipped <- .fix.clip.to.plausible( trial[ j, ], reference.spectra[ j, ],
                                             allow[ j, ], plausibility.reset )
          trial[ j, ] <- clipped$row
          out$implausible <- clipped$frac
          if ( clipped$frac > plausibility.reject ) out$reject <- "implausible"
        }

        out$trial <- trial

        if ( !is.null( af.gate.spectra ) )
          out$af.likeness <- .fix.af.likeness( trial[ j, ] - spectra.new[ j, ],
                                               af.gate.spectra )

        if ( is.na( out$reject ) && box && af.gate &&
             is.finite( out$af.likeness ) && out$af.likeness > af.gate.max )
          out$reject <- "af.like"

        if ( !is.na( out$reject ) ) return( out )

        design.trial <- if ( is.null( background.basis ) ) trial else
          rbind( background.basis, trial )

        out$condition.after <- calculate.condition.number( design.trial )

        # Matching seeds make "before" and "after" score the same held-out
        # events, so a difference reflects the design change only.
        leak.seed <- seed + f * 1009L

        set.seed( leak.seed )
        out$leak.before <- .fix.leakage( residual[ idx, , drop = FALSE ],
                                         design, j, leak.panel )
        set.seed( leak.seed )
        out$leak.after  <- .fix.leakage( residual[ idx, , drop = FALSE ],
                                         design.trial, j, leak.panel )

        out$reject <-
          if ( out$condition.after > max.condition.increase * condition.curr )
            "conditioning" else
          if ( is.finite( out$leak.before ) && is.finite( out$leak.after ) &&
               out$leak.after >= ( 1 + leakage.margin ) * out$leak.before )
            "leakage" else NA_character_

        out
      }

      judged         <- NULL
      candidate.type <- if ( not.joint && phase.two.fallback == "backsolved" )
        "backsolved" else "remeasured"
      rescued.from   <- NA_character_

      if ( is.na( reject ) ) {
        judged <- judge( if ( candidate.type == "backsolved" )
          spectra.backsolved[ j, ] else candidate$signature )
        reject <- judged$reject
      }

      # Phase one's back-solved row as a second candidate, when the first failed
      # a gate that says nothing against the row itself and phase one's evidence
      # for the row is strong.
      if ( !is.na( reject ) && phase.two.rescue == "backsolved" &&
           reject %in% rescue.reasons && flagged[ j ] &&
           source.weight[ j ] >= rescue.min.weight &&
           .signature.row.angle( spectra.new[ j, , drop = FALSE ],
                                 spectra.backsolved[ j, , drop = FALSE ] ) <= max.angle ) {

        second <- judge( spectra.backsolved[ j, ] )

        if ( is.na( second$reject ) ) {
          rescued.from   <- reject
          candidate.type <- "backsolved"
          judged         <- second
          reject         <- NA_character_
        }
      }

      # A tandem whose other candidates failed: its starting row with only the
      # parent fraction changed by phase one's tandem-to-parent coefficient.
      # Breakdown adds parent signal where a tandem has little, so the parent
      # direction is the expected one and the plausibility box does not apply.
      parent <- if ( j %in% names( tandems ) ) tandems[[ j ]] else NA_character_
      tandem.coefficient <- if ( is.na( parent ) ) NA_real_ else
        spillover.curr[ j, parent ]

      if ( !is.na( reject ) && tandem.candidate && !is.na( parent ) &&
           flagged[ j ] && abs( tandem.coefficient ) >= tandem.min.coefficient ) {

        target.row <- spectra.fluor[ j, ] + tandem.coefficient * spectra.fluor[ parent, ]

        if ( max( target.row ) > 0 ) {

          target.row <- target.row / max( target.row )

          if ( .signature.row.angle( spectra.new[ j, , drop = FALSE ],
                                     matrix( target.row, nrow = 1 ) ) <= max.angle ) {

            third <- judge( target.row, box = FALSE )

            if ( is.na( third$reject ) ) {
              rescued.from   <- reject
              candidate.type <- "parent"
              judged         <- third
              reject         <- NA_character_
            }
          }
        }
      }

      if ( is.na( reject ) ) {
        spectra.new    <- judged$trial
        design         <- if ( is.null( background.basis ) ) spectra.new else
          rbind( background.basis, spectra.new )
        condition.curr <- judged$condition.after
      }

      # The statistics block is taken from `extract.raw.signature()` when there
      # is one and built from the same column set when there is not, so a change
      # to what that function reports cannot silently break the row bind.
      stats.block <- if ( is.null( candidate ) )
        .fix.empty.signature.stats() else
          candidate$stats[ , setdiff( names( candidate$stats ),
                                      c( "fluorophore", "n.events" ) ),
                           drop = FALSE ]

      log.rows[[ length( log.rows ) + 1L ]] <- data.frame(
        fluorophore        = j,
        n.events           = length( idx ),
        stats.block,
        candidate          = if ( is.na( reject ) ) candidate.type else NA_character_,
        rescued.from       = rescued.from,
        tandem.parent      = parent,
        tandem.coefficient = tandem.coefficient,
        af.likeness        = if ( is.null( judged ) ) NA_real_ else judged$af.likeness,
        implausible.frac   = if ( is.null( judged ) ) NA_real_ else judged$implausible,
        leak.before        = if ( is.null( judged ) ) NA_real_ else judged$leak.before,
        leak.after         = if ( is.null( judged ) ) NA_real_ else judged$leak.after,
        condition.after    = if ( is.null( judged ) ) NA_real_ else judged$condition.after,
        accepted           = is.na( reject ),
        reason             = if ( is.na( reject ) ) "accepted" else reject,
        row.names          = NULL,
        stringsAsFactors   = FALSE )
    }

    signature.log <- do.call( rbind, log.rows )

    if ( verbose )
      message( sprintf(
        "\033[32mAccepted signature updates for %d of %d fluorophore(s)%s.\033[0m",
        sum( signature.log$accepted ), fluorophore.n,
        if ( any( signature.log$accepted ) ) paste0(
          ": ", paste( signature.log$fluorophore[ signature.log$accepted ],
                       collapse = ", " ) ) else "" ) )
  }

  accepted <- stats::setNames( rep( FALSE, fluorophore.n ), fluorophores )
  if ( !is.null( signature.log ) )
    accepted[ signature.log$fluorophore ] <- signature.log$accepted

  list(
    spectra                    = spectra.new,
    spectra.backsolved         = spectra.backsolved,
    spillover                  = spillover.curr,
    compensation               = compensation.curr,
    trust                      = trust,
    coefficient.log            = coefficient.log,
    source.log                 = source.log,
    signature.log              = signature.log,
    convergence.log            = convergence.log,
    stop.reason                = stop.reason,
    smooth.log                 = NULL,
    spillover.history          = spillover.history,
    unmixed.final              = unmixed.final,
    residual.final             = residual.final,
    thresholds.final           = thresholds.final,
    threshold.matrix.final     = threshold.matrix.final,
    neg.thresholds.final       = neg.thresholds.final,
    neg.threshold.matrix.final = neg.threshold.matrix.final,
    dominant.final             = dominant.final,
    af.basis                   = background.basis,
    af.hotspot                 = af.hotspot,
    af.frozen                  = af.frozen,
    leakage.prior              = prior,
    accepted                   = accepted
  )
}
