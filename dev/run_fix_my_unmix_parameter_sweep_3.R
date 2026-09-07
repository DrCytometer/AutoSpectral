# run_fix_my_unmix_parameter_sweep.R
#
# A parameter sweep for fix.my.unmix() / the Stage R alternating-pass
# prototype, run against the real bead/cell control pair, scored against
# ground truth exactly as run_fix_my_unmix_diagnostics.R's Stage C and
# Stage R already do. The goal is the same three questions for every
# argument fix.my.unmix() exposes (plus the pass-loop arguments Stage R
# adds ahead of merging it into production):
#
#   1. Does the metric move at all across a plausible range? If not, the
#      parameter is a candidate to hard-code and drop from the signature.
#   2. Where does it move, and in which direction? That sets the
#      recommended default.
#   3. Where does it break (identity collapses, an estimator error inflates,
#      or the true impact gets worse instead of better)? That sets a cap.
#
# PREREQUISITE - this script does not stand alone. It reuses the objects and
# `.fx.*` helper functions that run_fix_my_unmix_diagnostics.R builds, so
# that a two-substrate real dataset only has to be read, gated, and
# calibrated once per session:
#
#   - Source or run Section 0 (paths/objects) through Stage L (null
#     calibration) of run_fix_my_unmix_diagnostics.R first, in this same R
#     session. That leaves `fx.prep`, `fx.null`, `fx.stage.l`, `fx.variants`,
#     `af.name`, `fx.lapply`, and the helpers `.fx.angle()`, `.fx.renorm()`,
#     `.fx3.dominant()`, `.fx.induced.spillover()`, `.fx.stratified.sample()`
#     in scope. Stage R itself does not need to have been run.
#   - The package must be loaded (`devtools::load_all()`), since this script
#     calls `.fix.envelope.slope()`, `.fix.leakage()`, `extract.raw.signature()`,
#     `deconvolve.af.background()`, `get.spread.thresholds()`,
#     `get.variant.leakage.prior()`, `calculate.hotspot.matrix()`, and
#     `calculate.condition.number()` directly.
#
# WHAT THIS DOES NOT COVER YET, AND WHY - a fidelity note before the numbers
#
# run_fix_my_unmix_diagnostics.R's own Stage R is a simplified
# reimplementation of the phase-one pair loop: it hardcodes
# `quantiles = c(0.05, 0.5)`, a coverage floor of 0.10, a span floor of 5,
# and a flat `delta.max < 0.005` stop, and it never passes
# `min.bin.negative`, `spread.addback`, `anchor.weight`, `max.mask.passes`,
# `source.dominant`, or a leakage prior into `.fix.envelope.slope()` at all,
# so those arguments have silently been running at `.fix.envelope.slope()`'s
# own defaults (which are not all the same as fix.my.unmix()'s - notably
# `spread.addback` defaults to TRUE inside `.fix.envelope.slope()` and to
# FALSE in fix.my.unmix()). Stage R's phase-two loop is similarly partial: it
# never checks `min.bg.align`, `max.vif`, `max.condition.increase`,
# `allow.peak.shift`, or `max.hotspot`/AF-coupling freezing at all.
#
# None of that is a bug in the diagnostics script - it was built to be a
# cheap first look, not a parameter-accurate stand-in - but it does mean a
# sweep run through Stage R as it stands would recommend defaults for
# behaviour that was never actually exercised. So this script does not call
# Stage R's `.fx.estimate.spillover()`. It defines two new functions,
# `.fps.estimate.spillover()` and `.fps.run.pass.loop()`, that mirror
# fix.my.unmix()'s own phase-one and phase-two blocks argument-for-argument
# (envelope quantiles, the coverage/span/rise floors, the masking loop, the
# leakage prior, and the full phase-two gate cascade including
# background-confound, VIF, condition-number, and leakage checks), extended
# with only the four arguments that exist in the prototype and not yet in
# fix.my.unmix()'s signature (`n.pass`, `pass.events`, `step.decay`,
# `total.max.angle`). Everything reported below is therefore a property of
# the actual gates and the actual estimator, not of the shortcuts Stage R
# took to stay fast.
#
# Three defaults also disagree between fix_my_unmix.R and the diagnostics
# script's globals: `unstained.threshold` (0.99 in fix.my.unmix(), 0.995 in
# fx.threshold.probs), `max.iter` (20L vs fx.max.iter's 15L), and
# `convergence.threshold` (0.01 in fix.my.unmix(), a hardcoded 0.005 inside
# Stage R's `.fx.estimate.spillover()`). This script's `fps.defaults` uses
# fix.my.unmix()'s shipped values throughout, since establishing defaults is
# the point of the exercise; where a sweep shows the diagnostics script's
# value was the better one, that is itself a finding to report, not
# something to have assumed going in.
#
# STAGING
#
#   Stage P1 (cheap): phase-one only, one iterated fit from the starting
#   spectra, scored against the induced spillover truth exactly as Stage C
#   does. No signature extraction, so this is the estimator's own arguments
#   in isolation: the coverage/span/rise/disagreement floors, the masking
#   loop, the envelope shape, the leakage prior, and the convergence rule.
#   Wide grids are affordable here.
#
#   Stage P2 (expensive): the full alternating pass loop, scored against the
#   true panel impact exactly as Stage R's `pass.table$recovered` is. This is
#   where the phase-two gates, the pass-loop meta-parameters, and the
#   signature-extraction arguments live. Grids here are deliberately
#   narrower, and the script prints a timing estimate and waits for
#   `fps.confirm.p2 <- TRUE` before running the full grid, since this stage
#   costs roughly one full Stage R run per grid point.
#
#   NOT swept here: `downsample` / `downsample.background.frac` /
#   `downsample.min.stratum`. Changing the event budget means re-gating and
#   rebuilding the AF basis (`.fx.prepare()`'s KDE gate and
#   `get.af.basis()`'s PCA), which is Stage A's cost, not Stage R's, and
#   isn't worth re-paying inside an inner sweep loop. Test 2-3 budgets by
#   editing `fx.max.events` / `fx.background.frac` / `fx.min.stratum` and
#   re-running Stage A of the diagnostics script, then re-run this script's
#   Stage P2 once per budget.
#
#   Also not swept: `ridge` (1e-6 numerical stabiliser for the joint
#   signature fit - a priori not worth a slot in the grid, included once at
#   its default plus two order-of-magnitude bracket points to confirm it
#   really is inert) and `figures` / `save` / `verbose` (output plumbing,
#   not model behaviour).
#
# Run top to bottom. Both stages write one summary CSV and, for Stage P2, a
# per-pass detail CSV, to `fps.output.dir`.

# ---------------------------------------------------------------------------
# 0. Preconditions
# ---------------------------------------------------------------------------

fps.required <- c( "fx.prep", "fx.null", "fx.stage.l", "fx.variants", "af.name" )
fps.missing  <- fps.required[ !vapply( fps.required, exists, logical( 1 ) ) ]

if ( length( fps.missing ) > 0 )
  stop( sprintf( paste0(
    "Missing %s. Run Section 0 through Stage L of ",
    "run_fix_my_unmix_diagnostics.R in this session first." ),
    paste( fps.missing, collapse = ", " ) ), call. = FALSE )

fps.output.dir <- file.path( fx.output.dir, "parameter_sweep" )
if ( !dir.exists( fps.output.dir ) ) dir.create( fps.output.dir, recursive = TRUE )

fps.timing <- data.frame( stage = character( 0 ), unit = character( 0 ),
                          elapsed = numeric( 0 ), stringsAsFactors = FALSE )

#' @noRd
fps.log.time <- function( stage, unit, elapsed ) {
  fps.timing <<- rbind( fps.timing, data.frame(
    stage = stage, unit = unit, elapsed = unname( elapsed ),
    stringsAsFactors = FALSE ) )
}

#' @noRd
fps.write.csv <- function( df, name ) {
  utils::write.csv( df, file.path( fps.output.dir, name ), row.names = FALSE )
}

# ---------------------------------------------------------------------------
# 1. Parameter registry - every fix.my.unmix() argument this sweep touches,
#    at its shipped default, plus the four pass-loop-only prototype
#    arguments (n.pass, pass.events, step.decay, total.max.angle).
# ---------------------------------------------------------------------------

fps.defaults <- list(

  # Updated from fix.my.unmix()'s shipped defaults after two Stage P1 rounds
  # and one Stage P2 round: spread.addback -> TRUE held up under P1's
  # p90/median re-check. max.truncated.events went 20000 -> 40000 -> 20000
  # -> 40000: P1's max.abs.err liked 40000, P1's rms.err/p90.abs.err didn't
  # (round two), but P2's real ground-truth impact clearly does (Cells
  # 0.310 vs 0.172 at pass 3, Beads flat) - this is exactly the case P2 was
  # built to arbitrate, so it wins. n.pass -> 3 and step -> 0.5 are new from
  # this P2 round: continuing the outer loop past pass 3 was net-harmful on
  # Cells in 97% of pass-4 transitions and 100% of pass-5 ones regardless of
  # which gate was being swept, and step size showed a clean, monotonic,
  # both-substrates-agreeing improvement from 1 down to 0.5 with no sign of
  # plateauing yet - worth a follow-up sweep below 0.5 before treating that
  # value as final.

  # positivity / thresholds
  unstained.threshold   = 0.99,
  unstained.margin       = 1.3,
  spread.kappa           = 2,

  # phase-one pair estimator
  estimator               = "truncated",
  envelope.quantiles      = c( 0.05, 0.5 ),
  min.negative.frac       = 0.10,
  max.disagreement        = 0.5,
  min.negative.events     = 200L,
  min.bin.negative        = 25L,
  min.span                = 5,
  min.rise                = 1,
  n.levels.pair           = 10L,
  max.truncated.events    = 40000L,
  max.mask.passes         = 3L,
  source.dominant         = TRUE,
  spread.addback          = TRUE,
  anchor.weight           = 1,
  max.coefficient         = 0.2,
  leakage.prior           = TRUE,
  span.fraction           = 0.6,
  max.iter                = 20L,
  convergence.threshold   = 0.01,
  convergence.quantile    = 0.95,

  # phase-two signature extraction and gates
  n.levels                = 60L,
  min.bin.events           = 50L,
  multivariate             = TRUE,
  ridge                    = 1e-6,
  intercept                = TRUE,
  min.explained            = 0.8,
  max.explained            = 1.2,
  max.clamp.frac           = 0.15,
  min.impact.ratio         = 2,
  max.angle                = 10,
  min.bg.align             = -0.9,
  max.anchor               = 0.10,
  max.vif                  = 500,
  max.condition.increase   = 1.05,
  peak.shift.min.rel       = 0.7,
  max.hotspot              = 5,
  step                     = 0.5,

  # pass-loop-only (prototype; not yet in fix.my.unmix()'s own signature)
  n.pass                   = 3L,
  pass.events               = 40000L,
  step.decay                = 0.95,
  total.max.angle            = 15
)

# ---------------------------------------------------------------------------
# 2. Parameterised phase-one estimator - mirrors fix_my_unmix.R lines
#    ~531-731 argument-for-argument, operating on already-deconvolved
#    abundance rather than raw + background basis, since the caller (either
#    the P1 evaluator below or `.fps.run.pass.loop()`) has already run
#    `deconvolve.af.background()` with the spectra this pass/evaluation
#    starts from. No "best iteration" tracking: like production, it returns
#    whatever the loop was holding when it stopped, not the iteration with
#    the smallest delta.
# ---------------------------------------------------------------------------

#' @noRd
.fps.estimate.spillover <- function( prep, abundance, unstained.abundance,
                                     params, variants = NULL, m.init = NULL ) {

  panel <- prep$panel
  n.fl  <- length( panel )

  prior <- NULL
  if ( isTRUE( params$leakage.prior ) && !is.null( variants ) ) {
    af.row <- prep$s.start[ setdiff( rownames( prep$s.start ), panel ),
                            , drop = FALSE ]
    prior <- tryCatch(
      get.variant.leakage.prior(
        spectra       = rbind( prep$s.start.panel, af.row ),
        variants      = variants,
        extra.rows    = prep$af.basis,
        af.name       = af.name,
        span.fraction = params$span.fraction,
        verbose       = FALSE ),
      error = function( e ) NULL )
  }

  m.hat <- if ( !is.null( m.init ) ) m.init else diag( n.fl )
  dimnames( m.hat ) <- list( panel, panel )

  delta.history   <- rep( NA_real_, 3L )
  convergence.log <- data.frame( iter = integer( 0 ), delta = numeric( 0 ),
                                 delta.quantile = numeric( 0 ),
                                 delta.max = numeric( 0 ) )
  final.log   <- NULL
  n.iter.used <- 0L

  for ( iter in seq_len( as.integer( params$max.iter ) ) ) {

    comp <- tryCatch( solve( m.hat ), error = function( e ) NULL )
    if ( is.null( comp ) ) break

    unmixed <- abundance           %*% comp
    unst    <- unstained.abundance %*% comp

    thr <- params$unstained.margin *
      apply( unst, 2, stats::quantile, probs = params$unstained.threshold,
            names = FALSE )
    names( thr ) <- panel

    neg.var <- apply( unst, 2, stats::mad )^2
    names( neg.var ) <- panel

    thr.mat <- get.spread.thresholds(
      unmixed = unmixed, thresholds = thr,
      spillover.spread = prep$spillover.spread,
      spread.kappa = params$spread.kappa, verbose = FALSE )

    dominant <- NULL
    if ( isTRUE( params$source.dominant ) ) {

      dyn.range <- pmax( apply( unmixed, 2, stats::quantile, probs = 0.999,
                                names = FALSE ) - thr, .Machine$double.eps )
      dominance.score <- sweep( pmax( unmixed - thr.mat, 0 ), 2, dyn.range, "/" )
      dominant <- max.col( dominance.score, ties.method = "first" )
      dominant[ dominance.score[ cbind( seq_along( dominant ), dominant ) ] <= 0 ] <- 0L
    }

    slope <- diag( n.fl ); trust <- diag( n.fl )
    dimnames( slope ) <- dimnames( trust ) <- list( panel, panel )

    pairs <- data.frame(
      src = rep( panel, each = n.fl - 1L ),
      ch  = unlist( lapply( panel, function( s ) setdiff( panel, s ) ) ),
      stringsAsFactors = FALSE )

    ests <- fx.lapply( seq_len( nrow( pairs ) ), function( i ) {

      src <- pairs$src[ i ]; ch <- pairs$ch[ i ]

      sv <- 0
      if ( !is.null( prep$spillover.spread ) &&
           src %in% rownames( prep$spillover.spread ) &&
           ch  %in% colnames( prep$spillover.spread ) ) {
        sv <- prep$spillover.spread[ src, ch ]
        if ( !is.finite( sv ) ) sv <- 0
        sv <- max( sv, 0 )
      }

      .fix.envelope.slope(
        x.source             = unmixed[ , src ],
        x.target             = unmixed[ , ch ],
        threshold.source     = thr.mat[ , src ],
        threshold.target     = thr.mat[ , ch ],
        spread.var           = sv,
        neg.var              = neg.var[ ch ],
        source.mask          = if ( is.null( dominant ) ) NULL else
          dominant == match( src, panel ),
        quantiles             = params$envelope.quantiles,
        n.levels              = params$n.levels.pair,
        min.events            = params$min.negative.events,
        min.bin.negative      = params$min.bin.negative,
        spread.addback        = params$spread.addback,
        anchor.weight         = params$anchor.weight,
        max.truncated.events  = params$max.truncated.events,
        max.coefficient       = params$max.coefficient,
        max.mask.passes       = params$max.mask.passes )
    } )

    rows <- vector( "list", nrow( pairs ) )

    for ( i in seq_len( nrow( pairs ) ) ) {

      src <- pairs$src[ i ]; ch <- pairs$ch[ i ]
      est <- ests[[ i ]]

      w <- 0
      slope.use <- if ( is.null( est ) ) NA_real_ else
        if ( identical( params$estimator, "truncated" ) )
          est$slope.truncated else est$slope

      if ( !is.null( est ) && is.finite( slope.use ) &&
           est$coverage >= params$min.negative.frac &&
           abs( slope.use ) <= params$max.coefficient &&
           abs( slope.use ) * est$span >= params$min.rise * est$noise &&
           est$span > params$min.span * abs( thr[ src ] ) ) {

        prior.var <- if ( is.null( prior ) ) slope.use^2 else
          prior$variance[ src, ch ]
        if ( !is.finite( prior.var ) || prior.var <= 0 ) prior.var <- slope.use^2

        w <- if ( is.finite( est$se ) && est$se > 0 )
          prior.var / ( prior.var + est$se^2 ) else 1

        if ( is.finite( est$disagreement ) &&
             est$disagreement > params$max.disagreement )
          w <- w * params$max.disagreement / est$disagreement

        slope[ src, ch ] <- slope.use
      }

      trust[ src, ch ] <- w

      rows[[ i ]] <- data.frame(
        source = src, channel = ch, trust = w,
        slope.truncated = if ( is.null( est ) ) NA_real_ else est$slope.truncated,
        se           = if ( is.null( est ) ) NA_real_ else est$se,
        disagreement = if ( is.null( est ) ) NA_real_ else est$disagreement,
        coverage     = if ( is.null( est ) ) NA_real_ else est$coverage,
        span         = if ( is.null( est ) ) NA_real_ else est$span,
        row.names = NULL, stringsAsFactors = FALSE )
    }

    slope.error  <- slope - diag( n.fl )
    off.diagonal <- abs( slope.error[ row( slope.error ) != col( slope.error ) ] )

    delta          <- stats::sd( slope.error )
    delta.quantile <- unname( stats::quantile(
      off.diagonal, probs = params$convergence.quantile, names = FALSE ) )
    delta.max      <- max( off.diagonal )

    convergence.log[ nrow( convergence.log ) + 1L, ] <-
      list( iter, delta, delta.quantile, delta.max )
    n.iter.used <- iter
    final.log   <- do.call( rbind, rows )

    if ( delta.quantile < params$convergence.threshold ) break

    delta.history <- c( delta.history[ -1L ], delta )
    if ( all( is.finite( delta.history ) ) && mean( diff( delta.history ) ) >= 0 )
      break

    m.next <- m.hat + ( trust * slope.error ) %*% m.hat
    dn <- diag( m.next )
    if ( any( !is.finite( dn ) ) || any( dn <= 0 ) ) break
    m.hat <- sweep( m.next, 1, dn, "/" )
  }

  list( spillover = m.hat, log = final.log, convergence.log = convergence.log,
        n.iter = n.iter.used,
        delta.quantile = if ( nrow( convergence.log ) > 0 )
          tail( convergence.log$delta.quantile, 1 ) else NA_real_,
        delta.max = if ( nrow( convergence.log ) > 0 )
          tail( convergence.log$delta.max, 1 ) else NA_real_ )
}

# ---------------------------------------------------------------------------
# 3. STAGE P1 - phase-one only, scored against the induced spillover truth
# ---------------------------------------------------------------------------

#' @noRd
.fps.run.phase1.eval <- function( tg, params ) {

  p <- fx.prep[[ tg ]]

  # Reset before every evaluation, not once at the top of the sweep, so that
  # two grid points differ only in `params` and not also in where each one
  # happened to land in the random stream. `.fix.envelope.slope()`'s bulk
  # subsampling is the only source of randomness this touches.
  if ( exists( "fx.seed" ) ) set.seed( fx.seed )

  dec  <- deconvolve.af.background( p$raw,       p$s.start.panel, p$af.basis,
                                    af.name = NULL )
  decu <- deconvolve.af.background( p$unstained, p$s.start.panel, p$af.basis,
                                    af.name = NULL )

  t0  <- Sys.time()
  res <- .fps.estimate.spillover( p, dec$abundance, decu$abundance,
                                  params = params, variants = fx.variants[[ tg ]] )
  elapsed <- as.numeric( Sys.time() - t0, units = "secs" )

  induced <- .fx.induced.spillover( p )
  shared  <- intersect( rownames( res$spillover ), rownames( induced ) )
  err     <- res$spillover[ shared, shared ] - induced[ shared, shared ]
  off     <- err[ row( err ) != col( err ) ]

  # max.abs.err is a single-worst-pair statistic and swings on whichever one
  # pair sits closest to a cap or a coverage floor; median.abs.err and
  # p90.abs.err report the typical and near-worst behaviour so a jump in the
  # max alone doesn't get read as a change in the estimator generally.
  data.frame(
    target         = tg,
    max.abs.err    = max( abs( off ) ),
    p90.abs.err    = unname( stats::quantile( abs( off ), 0.90, names = FALSE ) ),
    median.abs.err = stats::median( abs( off ) ),
    rms.err        = sqrt( mean( off^2 ) ),
    n.fitted       = if ( is.null( res$log ) ) 0L else sum( res$log$trust > 0 ),
    n.iter         = res$n.iter,
    delta.quantile.final = res$delta.quantile,
    elapsed.sec    = elapsed,
    stringsAsFactors = FALSE )
}

fps.grid.p1 <- list(
  unstained.threshold  = c( 0.95, 0.97, 0.99, 0.995, 0.999 ),
  unstained.margin     = c( 1.0, 1.15, 1.3, 1.5, 2.0 ),
  spread.kappa         = c( 0, 1, 2, 3, 4 ),
  estimator            = c( "truncated", "envelope" ),
  min.negative.frac    = c( 0.02, 0.05, 0.10, 0.20, 0.35 ),
  max.disagreement     = c( 0.25, 0.5, 0.75, 1, 2 ),
  min.negative.events  = c( 50L, 100L, 200L, 400L, 800L ),
  min.bin.negative     = c( 10L, 15L, 25L, 40L, 60L ),
  min.span             = c( 2, 3, 5, 8, 12 ),
  min.rise             = c( 0, 0.5, 1, 2, 4 ),
  n.levels.pair        = c( 5L, 8L, 10L, 15L, 20L ),
  max.truncated.events = c( 5000L, 10000L, 20000L, 40000L ),
  max.mask.passes      = c( 0L, 1L, 2L, 3L, 5L ),
  source.dominant      = c( TRUE, FALSE ),
  spread.addback       = c( TRUE, FALSE ),
  anchor.weight        = c( 0.25, 0.5, 1, 2, 4 ),
  max.coefficient      = c( 0.1, 0.15, 0.2, 0.3, 0.4 ),
  leakage.prior        = c( TRUE, FALSE ),
  span.fraction        = c( 0.3, 0.45, 0.6, 0.75, 0.9 ),
  max.iter             = c( 5L, 10L, 15L, 20L, 30L ),
  convergence.threshold = c( 0.002, 0.005, 0.01, 0.02, 0.05 ),
  convergence.quantile  = c( 0.75, 0.85, 0.95, 0.99, 1.0 )
)

# envelope.quantiles is a pair, not a scalar, so it gets its own small grid
# rather than living inside fps.grid.p1.
fps.grid.p1.envelope.quantiles <- list(
  c( 0.05, 0.5 ), c( 0.10, 0.5 ), c( 0.02, 0.5 ), c( 0.05, 0.25 ), c( 0.05, 0.75 ) )

cat( "\n===================== STAGE P1: phase-one sweep =====================\n" )

fps.p1.rows <- list()

for ( tg in names( fx.prep ) ) {
  t0 <- Sys.time()
  fps.p1.rows[[ paste( "baseline", tg ) ]] <- cbind(
    parameter = "baseline", value = "default",
    .fps.run.phase1.eval( tg, fps.defaults ) )
  fps.log.time( "P1", paste( "baseline", tg ),
               as.numeric( Sys.time() - t0, units = "secs" ) )
}

for ( pname in names( fps.grid.p1 ) ) {
  for ( val in fps.grid.p1[[ pname ]] ) {

    params <- fps.defaults
    params[[ pname ]] <- val

    for ( tg in names( fx.prep ) ) {
      t0  <- Sys.time()
      row <- .fps.run.phase1.eval( tg, params )
      fps.log.time( "P1", paste( pname, val, tg ),
                   as.numeric( Sys.time() - t0, units = "secs" ) )

      key <- paste( pname, val, tg )
      fps.p1.rows[[ key ]] <- cbind( parameter = pname, value = as.character( val ), row )

      cat( sprintf( "  [P1] %-22s = %-10s (%s): max.err %.4f  rms %.4f  n.fitted %d  n.iter %d\n",
                    pname, as.character( val ), tg,
                    row$max.abs.err, row$rms.err, row$n.fitted, row$n.iter ) )
    }
  }
}

for ( q in fps.grid.p1.envelope.quantiles ) {

  params <- fps.defaults
  params$envelope.quantiles <- q
  label  <- paste( q, collapse = "/" )

  for ( tg in names( fx.prep ) ) {
    t0  <- Sys.time()
    row <- .fps.run.phase1.eval( tg, params )
    fps.log.time( "P1", paste( "envelope.quantiles", label, tg ),
                 as.numeric( Sys.time() - t0, units = "secs" ) )

    key <- paste( "envelope.quantiles", label, tg )
    fps.p1.rows[[ key ]] <- cbind( parameter = "envelope.quantiles", value = label, row )

    cat( sprintf( "  [P1] %-22s = %-10s (%s): max.err %.4f  rms %.4f  n.fitted %d  n.iter %d\n",
                  "envelope.quantiles", label, tg,
                  row$max.abs.err, row$rms.err, row$n.fitted, row$n.iter ) )
  }
}

fps.p1.table <- do.call( rbind, fps.p1.rows )
rownames( fps.p1.table ) <- NULL
fps.write.csv( fps.p1.table, "stage_P1_phase1_sweep.csv" )

cat( "\n  wrote stage_P1_phase1_sweep.csv. For each parameter, compare every\n" )
cat( "  row's max.abs.err / rms.err against the \"baseline\" row for the same\n" )
cat( "  target. A parameter whose whole grid sits within noise of baseline on\n" )
cat( "  both targets is a candidate to hard-code at its default and drop from\n" )
cat( "  fix.my.unmix()'s signature. Watch n.fitted too: a value that improves\n" )
cat( "  max.abs.err by starving the estimator down to one or two coefficients\n" )
cat( "  is not really an improvement.\n" )

# ---------------------------------------------------------------------------
# 4. STAGE P2 - the full alternating pass loop, scored against true impact
# ---------------------------------------------------------------------------
# Mirrors fix_my_unmix.R's phase-two gate cascade (span, explained, fit,
# offset, background.confound, negative.mass, bias, angle, background,
# collinear, peak.shift, conditioning, leakage) in the order production
# checks them, layered under the prototype's own settled/plateau/drift
# retirement and AF-hotspot freezing, exactly as CONTEXT_fix_my_unmix.md's
# "alternating phase-one/phase-two loop with per-row plateau retirement"
# describes it.

#' @noRd
.fps.run.pass.loop <- function( tg, params ) {

  p      <- fx.prep[[ tg ]]

  # Same fix as .fps.run.phase1.eval(): reset before the loop so two grid
  # points differ only in `params`, not also in where each one landed in the
  # random stream via `sample.int()` below.
  if ( exists( "fx.seed" ) ) set.seed( fx.seed )

  shared <- intersect( rownames( fx.stage.l[[ tg ]]$s.null ), p$panel )

  bias <- fx.stage.l[[ tg ]]$s.null[ shared, p$detectors, drop = FALSE ] -
    fx.null[[ tg ]]$s.start.panel[ shared, p$detectors, drop = FALSE ]

  m.null.bias <- NULL
  if ( !is.null( fx.stage.l[[ tg ]]$m.spillover ) ) {
    m.null.bias <- fx.stage.l[[ tg ]]$m.spillover -
      diag( nrow( fx.stage.l[[ tg ]]$m.spillover ) )
    dimnames( m.null.bias ) <- dimnames( fx.stage.l[[ tg ]]$m.spillover )
  }

  reference     <- fx.stage.l[[ tg ]]$signature
  cut.resid     <- 3 * stats::median( reference$resid.rel,     na.rm = TRUE )
  cut.intercept <- 3 * stats::median( reference$intercept.rel, na.rm = TRUE )

  af.hotspot <- calculate.hotspot.matrix( rbind( p$af.basis, p$s.start.panel ) )
  coupling   <- apply( af.hotspot[ p$panel, rownames( p$af.basis ), drop = FALSE ],
                       1, max )
  af.frozen  <- names( coupling )[ coupling > params$max.hotspot ]

  s.curr         <- p$s.start.panel
  design.curr    <- rbind( p$af.basis, s.curr )
  condition.curr <- calculate.condition.number( design.curr )

  ever.accepted <- character( 0 )
  settled       <- character( 0 )
  step.prev     <- stats::setNames( rep( Inf, length( p$panel ) ), p$panel )
  m.prev        <- NULL

  pass.rows <- list()
  gate.rows <- list()

  for ( pass in seq_len( as.integer( params$n.pass ) ) ) {

    dec  <- deconvolve.af.background( p$raw,       s.curr, p$af.basis, af.name = NULL )
    decu <- deconvolve.af.background( p$unstained, s.curr, p$af.basis, af.name = NULL )

    fit.idx <- if ( nrow( dec$abundance ) > params$pass.events )
      sample.int( nrow( dec$abundance ), params$pass.events ) else
        seq_len( nrow( dec$abundance ) )

    ph <- .fps.estimate.spillover(
      p, dec$abundance[ fit.idx, , drop = FALSE ], decu$abundance,
      params = params, variants = fx.variants[[ tg ]], m.init = m.prev )
    m.prev <- ph$spillover

    unmixed <- dec$abundance %*% solve( ph$spillover )

    thr <- params$unstained.margin *
      apply( decu$abundance %*% solve( ph$spillover ), 2, stats::quantile,
            probs = params$unstained.threshold, names = FALSE )
    names( thr ) <- p$panel

    thr.mat <- get.spread.thresholds(
      unmixed = unmixed, thresholds = thr,
      spillover.spread = p$spillover.spread,
      spread.kappa = params$spread.kappa, verbose = FALSE )

    dom <- .fx3.dominant( unmixed, thr, thr.mat )

    x.high <- apply( unmixed, 2, stats::quantile, probs = 0.999, names = FALSE )
    names( x.high ) <- p$panel

    s.raw      <- s.curr
    accepted   <- character( 0 )
    step.total <- 0

    for ( j in p$panel ) {

      reason <- NA_character_

      if ( j %in% af.frozen ) reason <- "af.coupled"

      idx <- if ( is.na( reason ) ) which( dom == match( j, p$panel ) ) else integer( 0 )
      if ( is.na( reason ) && length( idx ) < params$min.negative.events )
        reason <- "no.fit"

      cand <- NULL
      if ( is.na( reason ) ) {
        cand <- extract.raw.signature(
          raw.data       = dec$residual[ idx, , drop = FALSE ],
          spectra        = s.curr,
          abundance      = unmixed[ idx, , drop = FALSE ],
          target         = j,
          active         = p$panel,
          intercept      = params$intercept,
          multivariate   = params$multivariate,
          ridge          = params$ridge,
          n.levels       = params$n.levels,
          min.bin.events = params$min.bin.events,
          min.events     = params$min.negative.events,
          background.raw = if ( any( dom == 0L ) )
            dec$residual[ dom == 0L, , drop = FALSE ] else NULL )
        if ( is.null( cand ) ) reason <- "no.fit"
      }

      step.impact <- NA_real_
      bias.impact <- NA_real_

      if ( is.na( reason ) ) {

        st <- cand$stats

        unmixing.curr <- solve( tcrossprod( design.curr ), design.curr )[
          p$panel, , drop = FALSE ]

        step.impact <- unname( x.high[ j ] ) * sqrt( sum(
          ( ( cand$signature - s.curr[ j, ] ) %*% t( unmixing.curr ) )^2 ) )

        bias.impact <- if ( j %in% shared )
          unname( x.high[ j ] ) *
          sqrt( sum( ( bias[ j, ] %*% t( unmixing.curr ) )^2 ) ) else 0

        drift <- .fx.angle( s.curr[ j, , drop = FALSE ],
                            p$s.start.panel[ j, , drop = FALSE ] )

        # Same order fix.my.unmix() checks its own gate cascade in, with the
        # prototype's settled / plateau / drift retirement appended after
        # every production gate has already passed.
        if ( j %in% settled )                                          reason <- "settled" else
        if ( st$x.span <= params$min.span * abs( thr[ j ] ) )          reason <- "span" else
        if ( !is.finite( st$explained.total ) ||
             st$explained.total < params$min.explained ||
             st$explained.total > params$max.explained )               reason <- "explained" else
        if ( st$resid.rel > cut.resid )                                reason <- "fit" else
        if ( st$intercept.rel > cut.intercept )                        reason <- "offset" else
        if ( is.finite( st$bg.align ) && st$bg.align < params$min.bg.align )
                                                                         reason <- "background.confound" else
        if ( st$clamp.frac > params$max.clamp.frac )                   reason <- "negative.mass" else
        if ( is.finite( step.impact ) && is.finite( bias.impact ) &&
             step.impact < params$min.impact.ratio * bias.impact )     reason <- "bias" else
        if ( st$deg.change > params$max.angle )                        reason <- "angle" else
        if ( is.finite( st$anchor.rel ) && st$anchor.rel > params$max.anchor )
                                                                         reason <- "background" else
        if ( st$vif.target > params$max.vif )                          reason <- "collinear" else
          if ( st$peak.new != st$peak.curr &&
               st$peak.new.rel < params$peak.shift.min.rel ) reason <- "peak.shift" else
        if ( step.impact >= params$step.decay * step.prev[ j ] )       reason <- "plateau" else
        if ( drift > params$total.max.angle )                          reason <- "drift" else
                                                                         reason <- "accepted"

        if ( identical( reason, "plateau" ) ) settled <- c( settled, j )
      }

      gate.rows[[ length( gate.rows ) + 1L ]] <- data.frame(
        pass = pass, fluorophore = j, reason = reason,
        step.impact = step.impact, bias.impact = bias.impact,
        row.names = NULL, stringsAsFactors = FALSE )

      if ( !identical( reason, "accepted" ) ) next

      proposed <- pmax( cand$signature +
                          params$step * ( cand$signature - s.curr[ j, ] ), 0 )
      if ( max( proposed ) <= 0 ) next
      trial <- s.raw
      trial[ j, ] <- proposed / max( proposed )

      design.trial   <- rbind( p$af.basis, trial )
      condition.after <- calculate.condition.number( design.trial )

      leak.before <- .fix.leakage( dec$residual[ idx, , drop = FALSE ],
                                   design.curr, j, p$panel )
      leak.after  <- .fix.leakage( dec$residual[ idx, , drop = FALSE ],
                                   design.trial, j, p$panel )

      if ( condition.after > params$max.condition.increase * condition.curr ) next
      if ( is.finite( leak.before ) && is.finite( leak.after ) &&
           leak.after >= leak.before ) next

      step.prev[ j ] <- step.impact
      step.total     <- step.total + step.impact
      s.raw[ j, ]    <- trial[ j, ]
      accepted       <- c( accepted, j )
      design.curr    <- design.trial
      condition.curr <- condition.after
    }

    # A pass that accepts nothing still has phase-one's spillover.pass folded
    # back into s.curr below (back.raw <- spillover.pass %*% s.curr), which
    # reapplies phase-one's residual-error estimate with no new phase-two
    # signal to justify it. Across the first full P2 grid, a pass with
    # n.accepted == 0 made recovered worse on that same pass 85% of the time
    # (mean change -0.074), against +0.061 on average for a pass that
    # accepted anything - so stop here instead of folding this pass in.
    if ( length( accepted ) == 0L ) break

    ever.accepted <- union( ever.accepted, accepted )
    spillover.pass <- ph$spillover

    if ( !is.null( m.null.bias ) ) {
      shared.m <- intersect( rownames( m.null.bias ), p$panel )
      spillover.pass[ shared.m, shared.m ] <-
        spillover.pass[ shared.m, shared.m ] -
        m.null.bias[ shared.m, shared.m, drop = FALSE ]
      diag( spillover.pass ) <- 1
    }

    frozen.now <- union( af.frozen, setdiff( p$panel, ever.accepted ) )
    if ( length( frozen.now ) > 0 ) {
      spillover.pass[ frozen.now, ] <- 0
      spillover.pass[ cbind( frozen.now, frozen.now ) ] <- 1
    }

    back.raw <- spillover.pass %*% s.curr

    design.split <- rbind( p$af.basis, s.curr )
    perp      <- .signature.span.split( s.raw - s.curr, design.split )$perpendicular
    bias.perp <- .signature.span.split( bias, design.split )$perpendicular

    s.next <- back.raw + perp

    debias <- intersect( shared, accepted )
    if ( length( debias ) > 0 )
      s.next[ debias, ] <- s.next[ debias, , drop = FALSE ] -
      bias.perp[ debias, , drop = FALSE ]

    s.curr <- .fx.renorm( pmax( s.next, 0 ) )
    design.curr    <- rbind( p$af.basis, s.curr )
    condition.curr <- calculate.condition.number( design.curr )

    pass.rows[[ pass ]] <- list( s = s.curr, n.accepted = length( accepted ),
                                 delta.quantile = ph$delta.quantile,
                                 step = step.total,
                                 angle = sum( .fx.angle( s.curr, p$s.true.panel ) ) )
  }

  design.q   <- rbind( p$af.basis, p$s.true.panel )
  unmixing.q <- solve( tcrossprod( design.q ), design.q )[ p$panel, , drop = FALSE ]

  dec.q <- deconvolve.af.background( p$raw, p$s.start.panel, p$af.basis, af.name = NULL )
  x.max <- apply( dec.q$abundance, 2, stats::quantile, probs = 0.999, names = FALSE )

  impact <- function( s ) sum( x.max * sqrt( rowSums(
    ( ( s[ p$panel, , drop = FALSE ] - p$s.true.panel ) %*% t( unmixing.q ) )^2 ) ) )

  impact.start <- impact( p$s.start.panel )

  pass.table <- data.frame(
    pass       = seq_along( pass.rows ),
    n.accepted = vapply( pass.rows, function( r ) r$n.accepted, integer( 1 ) ),
    delta.quantile = vapply( pass.rows, function( r ) r$delta.quantile, numeric( 1 ) ),
    step       = vapply( pass.rows, function( r ) r$step, numeric( 1 ) ),
    angle      = vapply( pass.rows, function( r ) r$angle, numeric( 1 ) ),
    impact     = vapply( pass.rows, function( r ) impact( r$s ), numeric( 1 ) ),
    row.names = NULL )

  pass.table$recovered <- 1 - pass.table$impact / impact.start

  list( pass.table = pass.table, gate.detail = do.call( rbind, gate.rows ),
        impact.start = impact.start, n.frozen = length( af.frozen ) )
}

# ---------------------------------------------------------------------------
# 5. STAGE P2 driver - narrower grids, timing estimate + confirmation gate
# ---------------------------------------------------------------------------

fps.grid.p2 <- list(
  n.pass             = c( 2L, 3L, 5L, 8L ),
  pass.events        = c( 15000L, 25000L, 40000L, 60000L ),
  step.decay         = c( 0.85, 0.90, 0.95, 0.99 ),
  total.max.angle    = c( 8, 12, 15, 20 ),
  min.impact.ratio   = c( 1, 1.5, 2, 3 ),
  min.explained      = c( 0.6, 0.7, 0.8, 0.9 ),
  max.explained      = c( 1.1, 1.2, 1.4, 1.6 ),
  max.clamp.frac     = c( 0.05, 0.10, 0.15, 0.25 ),
  max.angle          = c( 5, 10, 15, 20 ),
  min.bg.align       = c( -0.99, -0.95, -0.9, -0.7 ),
  max.anchor         = c( 0.05, 0.10, 0.15, 0.25 ),
  max.vif            = c( 50, 200, 500, 1000 ),
  max.condition.increase = c( 1.02, 1.05, 1.10, 1.20 ),
  peak.shift.min.rel = c( Inf, 0.7, 0.5 ),
  max.hotspot        = c( 3, 5, 8, 12 ),
  step               = c( 0.5, 0.75, 1 ),
  n.levels           = c( 30L, 45L, 60L, 90L ),
  min.bin.events     = c( 25L, 50L, 75L, 100L ),
  multivariate       = c( TRUE, FALSE ),
  intercept          = c( TRUE, FALSE ),
  ridge              = c( 1e-8, 1e-6, 1e-4 ),

  # Phase-one parameters that survived Stage P1's p90.abs.err / rms.err
  # re-check (round two), not just the max.abs.err pass from round one.
  # max.mask.passes, source.dominant, span.fraction, and max.coefficient
  # were all in this list after round one and were dropped after round two
  # showed their round-one signal was a single volatile pair, not a real
  # effect - see the sweep_timing notes / conversation log for that
  # withdrawal. min.rise is the one promotion: a real, substrate-divergent
  # rms.err effect with a large coverage swing on Beads (39 -> 240 fitted
  # pairs), which P1's induced-spillover metric can't fully arbitrate since
  # many of those extra pairs may have a near-zero true coefficient either
  # way - P2's abundance-weighted impact can.
  min.rise             = c( 0, 0.5, 1 ),
  max.truncated.events = c( 20000L, 40000L ),
  spread.kappa         = c( 1, 2, 3 )
)

fps.n.p2.runs <- sum( vapply( fps.grid.p2, length, integer( 1 ) ) ) *
  length( names( fx.prep ) ) + length( names( fx.prep ) )

cat( "\n===================== STAGE P2: full pass-loop sweep =====================\n" )
cat( sprintf( "  grid covers %d parameters, %d total runs (including baseline) per target.\n",
             length( fps.grid.p2 ), fps.n.p2.runs / length( fx.prep ) ) )

t0 <- Sys.time()
fps.p2.baseline <- lapply( names( fx.prep ), function( tg )
  .fps.run.pass.loop( tg, fps.defaults ) )
names( fps.p2.baseline ) <- names( fx.prep )
fps.baseline.elapsed <- as.numeric( Sys.time() - t0, units = "secs" )

cat( sprintf( paste0(
  "  one baseline pass loop over both targets took %.1fs. At that rate the\n",
  "  full grid above is roughly %.1f minutes. Edit fps.grid.p2 to trim it, or\n",
  "  set fps.confirm.p2 <- TRUE and re-run this section to proceed as-is.\n" ),
  fps.baseline.elapsed,
  fps.baseline.elapsed / length( fx.prep ) * fps.n.p2.runs / 60 ) )

if ( !exists( "fps.confirm.p2" ) ) fps.confirm.p2 <- FALSE

fps.p2.summary.rows <- list()
fps.p2.detail.rows  <- list()

#' @noRd
.fps.p2.record <- function( pname, value, tg, result ) {

  pt <- result$pass.table
  best <- which.min( pt$impact )

  fps.p2.summary.rows[[ length( fps.p2.summary.rows ) + 1L ]] <<- data.frame(
    parameter        = pname, value = as.character( value ), target = tg,
    n.pass.run       = nrow( pt ),
    best.pass        = pt$pass[ best ],
    best.recovered   = pt$recovered[ best ],
    final.recovered  = pt$recovered[ nrow( pt ) ],
    n.frozen         = result$n.frozen,
    impact.start     = result$impact.start,
    row.names = NULL, stringsAsFactors = FALSE )

  fps.p2.detail.rows[[ length( fps.p2.detail.rows ) + 1L ]] <<- cbind(
    parameter = pname, value = as.character( value ), target = tg, pt )
}

for ( tg in names( fx.prep ) )
  .fps.p2.record( "baseline", "default", tg, fps.p2.baseline[[ tg ]] )

if ( fps.confirm.p2 ) {

  for ( pname in names( fps.grid.p2 ) ) {
    for ( val in fps.grid.p2[[ pname ]] ) {

      params <- fps.defaults
      params[[ pname ]] <- val

      for ( tg in names( fx.prep ) ) {

        t0     <- Sys.time()
        result <- .fps.run.pass.loop( tg, params )
        fps.log.time( "P2", paste( pname, val, tg ),
                     as.numeric( Sys.time() - t0, units = "secs" ) )

        .fps.p2.record( pname, val, tg, result )

        pt <- result$pass.table
        cat( sprintf(
          "  [P2] %-22s = %-10s (%s): final recovered %.3f, best %.3f at pass %d, %d frozen\n",
          pname, as.character( val ), tg,
          tail( pt$recovered, 1 ), max( pt$recovered ),
          pt$pass[ which.max( pt$recovered ) ], result$n.frozen ) )
      }
    }
  }

  fps.p2.summary <- do.call( rbind, fps.p2.summary.rows )
  fps.p2.detail  <- do.call( rbind, fps.p2.detail.rows )

  fps.write.csv( fps.p2.summary, "stage_P2_full_sweep_summary.csv" )
  fps.write.csv( fps.p2.detail,  "stage_P2_full_sweep_detail.csv" )

  cat( "\n  wrote stage_P2_full_sweep_summary.csv and stage_P2_full_sweep_detail.csv.\n" )
  cat( "  Compare best.recovered / final.recovered against the \"baseline\" row per\n" )
  cat( "  target. Where best.pass < n.pass.run consistently, step.decay or\n" )
  cat( "  total.max.angle is retiring rows before they should stop, or a gate is\n" )
  cat( "  intermittently re-admitting a row that should have settled - check\n" )
  cat( "  stage_P2_full_sweep_detail.csv's per-pass angle column for that value.\n" )

} else {

  fps.write.csv( do.call( rbind, fps.p2.summary.rows ), "stage_P2_baseline_only.csv" )
  cat( "\n  fps.confirm.p2 is FALSE: only the baseline pass loop ran, written to\n" )
  cat( "  stage_P2_baseline_only.csv. Set fps.confirm.p2 <- TRUE and re-run this\n" )
  cat( "  section (Stage P2 driver) to run the full grid above.\n" )
}

# ---------------------------------------------------------------------------
# Timing summary
# ---------------------------------------------------------------------------

fps.write.csv( fps.timing, "sweep_timing.csv" )
cat( sprintf( "\nAll sweep tables written to %s\n", fps.output.dir ) )
