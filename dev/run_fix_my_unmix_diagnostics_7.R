# run_fix_my_unmix_diagnostics.R
#
# Unified staged diagnostics for the in-span signature correction, run against
# the real bead / cell single-stained control pair, in both directions. The
# bead-versus-cell spectral difference is the error being corrected, and both
# sets of spectra were extracted independently, so the target is measured
# rather than imposed.
#
# This replaces run_fix_my_unmix_diagnostics_2_1.R and _3_1.R. What moved:
#
#   - Stage K (envelope-vs-truncated-vs-four-ingredient arm sweep) is gone.
#     Settled: the truncated estimator beats the envelope on both substrates,
#     confirmed twice (CONTEXT_fix_my_unmix.md 2.1). `.fx3.slope()`,
#     `fx3.arms`, `.fx3.sweep()`, `.fx3.summarise()`, `fx3.dim.dyes`, and
#     `.fx.truncated.slope()` all existed only to run or feed that sweep and
#     are not reproduced here.
#   - Stage D's three-way raw/back/hybrid comparison is gone. Settled: the
#     hybrid is the right combination (CONTEXT_fix_my_unmix.md, section 3).
#     Stage R computes the hybrid directly every pass; there is no separate
#     single-pass version to compare it against any more.
#   - Stage E (gate behaviour by arm) and Stage P (proposed gate stack, scored)
#     are gone. The gate stack they were used to tune is settled and is what
#     Stage R applies every pass (CONTEXT_fix_my_unmix.md, section 4). Per-row
#     gate reasons are still printed, from Stage R's own `gate.rows`.
#   - Stage M (deg.start/deg.bias ratio) and Stage N (bias-subtraction vs
#     gated-rule sweep, k = 1/2/3/5) are gone. Settled: bias subtraction beats
#     every gated rule tried (CONTEXT_fix_my_unmix.md 2.2); it is what Stage R
#     applies to both the spillover matrix and the signatures, not something
#     re-derived per run.
#   - Stage O (perp.frac as a gate) is gone. Explicitly rejected
#     (CONTEXT_fix_my_unmix.md 2.2): it predicts how much subtraction will
#     help, not whether a correction should be trusted.
#   - Stage Q (single-pass impact-weighted error) is gone. Stage R's
#     `pass.table` reports impact per pass, including pass 1, which is a
#     superset of what Stage Q showed.
#   - Stage S (autofluorescence coefficients as extraction regressors) is
#     gone. Settled and closed, rejected at the bin count that supports it
#     (CONTEXT_fix_my_unmix_session_summary.md, section 5).
#
# What stayed, and why:
#
#   - Stage A (error structure) and Stage B (AF basis diagnostics): cheap,
#     single-pass, and the right first look at a new dataset.
#   - Stage C (spillover vs induced truth): simplified to report the
#     estimator's own iterated result, not a separate single-pass fit,
#     removing a second, larger implementation of the same pair loop.
#   - Stage L (null calibration): kept in full. It is where `bias` and
#     `m.spillover` come from for Stage R's debiasing, and it is the one
#     stage explicitly flagged as needing to transfer to a second dataset
#     (CONTEXT_fix_my_unmix_session_summary.md, section 6, item 1).
#   - Stage R (alternating passes): the cost centre and the current
#     candidate algorithm. Kept in full, now fed by the stratified sample
#     from `.fx.prepare()` rather than a flat downsample.
#   - Stage T (spillover row audit) and Stage U (what is left, oblique
#     split): kept, now pointed at Stage R's final pass rather than the
#     retired Stage D.
#   - Stage F (spectral trace plots): kept, off by default (`fx.write.traces`),
#     simplified to start / final / truth since there are no separate raw and
#     back arms to plot any more.
#
# A bug fix carried over from the old harness: `.fx.estimate.spillover()`'s
# per-pair log recorded `slope` (the envelope estimate) but not
# `slope.truncated`, even though `fx.estimator = "truncated"` is what the
# matrix update actually uses. Production `fix.my.unmix()`'s own pair log
# already carries both. Fixed here so Stage C and Stage L's "largest
# coefficients invented on the null" table report the estimate the matrix
# was actually built from.
#
# Downsampling is now stratified by dominant fluorophore rather than a flat
# random draw (see `.fx.stratified.sample()`). In a concatenated single-stain
# pool only one tube in sixteen is positive for any one channel, so a uniform
# sample carries forward the same lopsided mix the raw file had: a great deal
# of "unstained for this channel" bulk from the other fifteen tubes, and only
# a thin slice of the population that actually carries the signal the pair
# estimator and the signature fit need. Stratifying gives every fluorophore's
# own positive population a floor, caps the redundant bulk separately, and
# lets `fx.max.events` come down substantially (200,000 to 80,000 below) for
# the same or better coverage of the populations that matter, which pays for
# itself every one of the five passes in Stage R re-derives its populations.
#
# Run top to bottom. Every stage prints a table; those tables are the output
# to report back.

# ---------------------------------------------------------------------------
# 0. Setup - EDIT THIS SECTION
# ---------------------------------------------------------------------------

asp <- get.autospectral.param( cytometer = "aurora" )

fx.spectra <- list(
  Cells = cell.spectra,
  Beads = bd.spectra
)

fx.raw <- list(
  Cells = concat.cells,
  Beads = concat.beads
)

fx.unstained <- list(
  Cells = readFCS( file.path( cell.dir, "A10 unstained_010_Cells.fcs" ) ),
  Beads = readFCS( file.path( bd.dir,   "A10 unstained_010_Beads_1.fcs" ) )
)

# Variant lists from get.spectral.variants(), one per particle type. Used for
# spillover.spread.
fx.variants <- list(
  Cells = readRDS( "./figure_spectral_variants/Spectral_variants_cells.rds" ),
  Beads = readRDS( "./figure_spectral_variants/Spectral_variants_beads.rds" )
)

af.name <- "AF"

fx.threshold.probs  <- 0.995
fx.threshold.margin <- 1.3
fx.spread.kappa     <- 2
fx.seed             <- 42L
fx.max.iter         <- 15L
fx.n.levels.pair    <- 10L
fx.min.events       <- 200L
fx.estimator        <- "truncated"
fx.max.coefficient  <- 0.2
fx.max.truncated    <- 5000L

# Stratified downsampling. fx.max.events is the total event budget out of
# .fx.prepare(); fx.background.frac is the share of that budget reserved for
# events dominant for nothing; fx.min.stratum is the floor below which a
# fluorophore's own positive population is kept whole rather than thinned.
# See `.fx.stratified.sample()`.
fx.max.events       <- 80000L
fx.background.frac  <- 0.3
fx.min.stratum      <- 3000L

# Stage R. fx.pass.events is a secondary safety cap on top of the stratified
# pool, applied only to the per-pass spillover fit, since a spillover
# coefficient is a population statistic and does not need every event.
fx.n.pass           <- 5L
fx.pass.events       <- 40000L
fx.impact.ratio      <- 2
fx.step.decay        <- 0.95
fx.total.max.angle   <- 15

fx.audit.n       <- 4L
fx.write.traces  <- TRUE
trace.dir        <- "./figure_fix_my_unmix"

# ---------------------------------------------------------------------------
# Parallel backend for the pair loop
# ---------------------------------------------------------------------------
# `.fx.estimate.spillover()`'s inner loop is every (source, channel) pair,
# every iteration, every pass, and each pair is independent of every other:
# nothing is shared but read-only inputs. That makes it an easy target for
# fork-based parallelism on macOS and Linux, which `.fix.envelope.slope()`
# now supports cleanly (see the fit.trace patch in fix_my_unmix.R) because
# it no longer makes any compiled linear-algebra call, so there is nothing
# for a forked worker to inherit mid-computation from a threaded BLAS. This
# is a local, prescheduled wrapper rather than the package's own
# `create.parallel.lapply()`: that helper defaults to one fork per task
# (`mc.preschedule = FALSE`), which suits its usual job of a handful of
# large, unevenly-sized per-file tasks, but here there are 200-1500 small,
# near-uniform-cost tasks per call, where prescheduling into `fx.mc.cores`
# chunks avoids paying fork overhead once per pair.
#
# Set fx.parallel <- FALSE to force sequential execution (also the automatic
# fallback on Windows, where fork-based mclapply is unavailable).

fx.parallel  <- TRUE
fx.mc.cores  <- if ( fx.parallel && .Platform$OS.type == "unix" &&
                     parallel::detectCores() > 1L )
  max( 1L, parallel::detectCores() - 1L ) else 1L

#' @noRd
fx.lapply <- function( x, fun, ... ) {

  if ( fx.mc.cores <= 1L ) return( lapply( x, fun, ... ) )

  # fit.truncated()'s call into MASS::rlm() does real, repeated linear
  # algebra (an IRLS fit re-solves a weighted least squares system every
  # iteration), not the single closed-form calculation fit.trace() was
  # reduced to. Left alone, mc.cores forked children each potentially
  # running multi-threaded Accelerate underneath that oversubscribes
  # whatever cores are available, which shows up as real but underwhelming
  # speedup rather than a crash or a clean failure. Pinning each child to a
  # single BLAS thread before it runs anything is what
  # create.parallel.lapply()'s own mclapply branch already does for this
  # reason; this wrapper needs the same guard since it does not go through
  # that function.
  wrapper.fun <- function( ... ) {
    if ( requireNamespace( "RhpcBLASctl", quietly = TRUE ) ) {
      RhpcBLASctl::blas_set_num_threads( 1L )
      RhpcBLASctl::omp_set_num_threads( 1L )
    }
    fun( ... )
  }

  parallel::mclapply( x, wrapper.fun, mc.cores = fx.mc.cores,
                      mc.preschedule = TRUE )
}

# ---------------------------------------------------------------------------
# Timing and CSV output
# ---------------------------------------------------------------------------
# Every stage's per-substrate for loop, and Stage R's per-pass loop, is timed
# with system.time() and logged below rather than left to be read off the
# console. Each stage's key table is written to fx.output.dir as its own csv
# so results can be handed over as files rather than copy-pasted output;
# fx.timing itself is written last, once every stage has logged into it.

fx.output.dir <- "./figure_fix_my_unmix_output"
if ( !dir.exists( fx.output.dir ) ) dir.create( fx.output.dir, recursive = TRUE )

fx.timing <- data.frame( stage = character( 0 ), unit = character( 0 ),
                         elapsed = numeric( 0 ), stringsAsFactors = FALSE )

#' @noRd
fx.log.time <- function( stage, unit, elapsed ) {
  fx.timing <<- rbind( fx.timing, data.frame(
    stage = stage, unit = unit, elapsed = unname( elapsed ),
    stringsAsFactors = FALSE ) )
  cat( sprintf( "  [%s / %s: %.1fs]\n", stage, unit, elapsed ) )
}

#' @noRd
fx.write.csv <- function( df, name ) {
  utils::write.csv( df, file.path( fx.output.dir, name ), row.names = FALSE )
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

#' @noRd
.fx.renorm <- function( s ) {
  row.max <- apply( s, 1, max )
  row.max[ !is.finite( row.max ) | row.max <= 0 ] <- 1
  s / row.max
}

#' @noRd
.fx.angle <- function( a, b ) {
  a <- as.matrix( a ); b <- as.matrix( b )
  cs <- rowSums( a * b ) / ( sqrt( rowSums( a^2 ) ) * sqrt( rowSums( b^2 ) ) )
  180 / pi * acos( pmin( 1, pmax( -1, cs ) ) )
}

#' Main-population scatter gate, density relative so it needs no landmarks. On
#' bead files this removes the noise and debris fraction, whose distinct
#' background profile otherwise enters the dominance populations and the
#' evaluation medians.
#' @noRd
.fx.main.gate <- function( scatter, gate.level = 0.1, grid.n = 128L,
                           max.events = 100000L ) {
  sc <- as.matrix( scatter[ , 1:2, drop = FALSE ] )
  fit.idx <- seq_len( nrow( sc ) )
  if ( length( fit.idx ) > max.events )
    fit.idx <- sample( fit.idx, max.events )
  bw <- vapply( 1:2, function( k ) {
    b <- tryCatch( MASS::bandwidth.nrd( sc[ fit.idx, k ] ),
                   error = function( e ) 0 )
    if ( !is.finite( b ) || b <= 0 )
      b <- 4 * 1.06 * stats::sd( sc[ fit.idx, k ] ) * length( fit.idx )^( -0.2 )
    max( b, .Machine$double.eps )
  }, numeric( 1 ) )
  kde  <- MASS::kde2d( sc[ fit.idx, 1 ], sc[ fit.idx, 2 ], h = bw, n = grid.n )
  ix   <- findInterval( sc[ , 1 ], kde$x, all.inside = TRUE )
  iy   <- findInterval( sc[ , 2 ], kde$y, all.inside = TRUE )
  kde$z[ cbind( ix, iy ) ] >= gate.level * max( kde$z )
}

#' Stratified subsample of a gated raw pool by per-fluorophore positivity. A
#' uniform downsample keeps every population in the same proportion the file
#' already had it in; in a concatenated single-stain pool that means the
#' fifteen tubes that read as background in any one channel set how hard the
#' sixteenth tube's own positive population gets thinned as well, and in a
#' fully stained sample the same thing happens because most panels have one
#' or two large lineage populations that dwarf everything else. Either way
#' that is backwards: the positive population is exactly what the pair
#' estimator and the signature fit need the most events from.
#'
#' Membership is by positivity, not by a single winner-take-all dominance
#' label: an event above threshold for two fluorophores at once belongs to
#' both strata, not whichever scores higher. That matters specifically for a
#' fully stained sample, where co-expression is real and is exactly the
#' population the `coverage` identifiability gate is checking; a
#' mutually-exclusive assignment would credit a co-positive event to only one
#' of the two markers and let a large stratum's subsampling discard it from
#' the other's population by chance. On a concatenated single-stain pool the
#' distinction is close to moot, since the spread-scaled boundary already
#' keeps one tube's events out of another's positive population except right
#' at the edge, so this reduces to the same partition as before there.
#'
#' Every fluorophore's stratum is sampled against its own quota: a floor that
#' keeps a dim dye's whole positive population intact whenever it is smaller
#' than the floor, and a proportional share of what is left of the budget for
#' strata that can use more. The chosen events are the union across
#' fluorophores, so a co-positive event drawn by more than one stratum is
#' kept once and effectively has a higher retention chance, which is the
#' right direction. Events positive for nothing are capped separately as a
#' flat fraction of the total, since additional events there mostly refine a
#' threshold that is already well determined.
#' @noRd
.fx.stratified.sample <- function( abundance, thresholds, threshold.matrix,
                                   n.total, background.frac = 0.3,
                                   min.stratum = 2000L ) {

  n <- nrow( abundance )
  if ( n <= n.total ) return( seq_len( n ) )

  positive <- ( abundance - threshold.matrix ) > 0

  background.idx <- which( rowSums( positive ) == 0L )
  positive.idx    <- lapply( seq_len( ncol( abundance ) ),
                             function( j ) which( positive[ , j ] ) )

  n.background <- min( length( background.idx ),
                       round( background.frac * n.total ) )
  budget.positive <- n.total - n.background

  sizes   <- vapply( positive.idx, length, integer( 1 ) )
  floor.n <- pmin( sizes, as.integer( min.stratum ) )

  # A panel with many fluorophores can demand more floor than the budget
  # holds; when it does, every stratum's floor is scaled down by the same
  # factor rather than starving whichever strata come later in column order.
  quota <- if ( sum( floor.n ) > budget.positive ) {
    pmin( sizes, floor( floor.n * budget.positive / max( sum( floor.n ), 1 ) ) )
  } else {
    spare           <- pmax( sizes - floor.n, 0 )
    leftover.budget <- budget.positive - sum( floor.n )
    extra <- if ( sum( spare ) > 0 )
      leftover.budget * spare / sum( spare ) else rep( 0, length( sizes ) )
    pmin( sizes, floor.n + round( extra ) )
  }

  chosen.positive <- unique( unlist( Map( function( idx, k )
    if ( length( idx ) <= k ) idx else sample( idx, k ),
    positive.idx, quota ), use.names = FALSE ) )

  chosen.background <- if ( n.background >= length( background.idx ) )
    background.idx else sample( background.idx, n.background )

  c( chosen.positive, chosen.background )
}

#' Assemble everything one direction of the test needs: the gated raw data
#' (stratified down to `fx.max.events`), the gated unstained, the starting
#' and target spectra restricted to shared rows and detectors, and the
#' autofluorescence basis. The scatter gate now runs on the full file before
#' any sampling, so the event budget is not spent on debris that would have
#' been gated out anyway.
#' @noRd
.fx.prepare <- function( target, start ) {

  shared.fluor <- intersect( rownames( fx.spectra[[ target ]] ),
                             rownames( fx.spectra[[ start ]] ) )
  shared.det   <- intersect( colnames( fx.spectra[[ target ]] ),
                             colnames( fx.spectra[[ start ]] ) )

  s.start <- as.matrix( fx.spectra[[ start  ]][ shared.fluor, shared.det, drop = FALSE ] )
  s.true  <- as.matrix( fx.spectra[[ target ]][ shared.fluor, shared.det, drop = FALSE ] )

  panel <- setdiff( shared.fluor, af.name )
  s.start.panel <- s.start[ panel, , drop = FALSE ]

  scatter.det <- grep( "^(FSC|SSC)", colnames( fx.raw[[ target ]] ), value = TRUE )

  gate.keep <- .fx.main.gate( fx.raw[[ target ]][ , scatter.det, drop = FALSE ] )
  gated.idx <- which( gate.keep )

  unst         <- fx.unstained[[ target ]][ , shared.det,  drop = FALSE ]
  unst.scatter <- fx.unstained[[ target ]][ , scatter.det, drop = FALSE ]
  unst.keep    <- .fx.main.gate( unst.scatter )
  unst         <- unst[ unst.keep, , drop = FALSE ]

  af.basis <- get.af.basis( unst, n.pc = "auto", verbose = FALSE )

  raw.gated    <- fx.raw[[ target ]][ gated.idx, shared.det, drop = FALSE ]
  design.rough <- rbind( af.basis, s.start.panel )

  abundance.full <- unmix.ols.fast( raw.gated, design.rough )
  colnames( abundance.full ) <- rownames( design.rough )
  abundance.rough <- abundance.full[ , panel, drop = FALSE ]

  unst.full <- unmix.ols.fast( unst, design.rough )
  colnames( unst.full ) <- rownames( design.rough )
  unst.abundance.rough <- unst.full[ , panel, drop = FALSE ]

  thr <- fx.threshold.margin * apply( unst.abundance.rough, 2, stats::quantile,
                                      probs = fx.threshold.probs, names = FALSE )
  names( thr ) <- panel

  thr.mat <- get.spread.thresholds(
    unmixed = abundance.rough, thresholds = thr,
    spillover.spread = fx.variants[[ target ]]$spillover.spread,
    spread.kappa = fx.spread.kappa, verbose = FALSE )

  set.seed( fx.seed )
  sample.idx <- .fx.stratified.sample(
    abundance = abundance.rough, thresholds = thr, threshold.matrix = thr.mat,
    n.total = fx.max.events, background.frac = fx.background.frac,
    min.stratum = fx.min.stratum )

  raw <- raw.gated[ sample.idx, , drop = FALSE ]

  list( target = target, start = start,
        panel = panel, detectors = shared.det,
        s.start = s.start, s.true = s.true,
        s.start.panel = s.start.panel,
        s.true.panel  = s.true[ panel, , drop = FALSE ],
        raw = raw, unstained = unst, af.basis = af.basis,
        gate.frac = mean( gate.keep ),
        spillover.spread = fx.variants[[ target ]]$spillover.spread )
}

#' The spillover actually induced by the difference between the two measured
#' spectra sets, under the design the data are unmixed with. This, not any
#' nominal mixing matrix, is what phase one can recover: it is the projection
#' of the real spectral error onto the row space, in unmixed coordinates.
#' @noRd
.fx.induced.spillover <- function( prep ) {

  design   <- rbind( prep$af.basis, prep$s.start.panel )
  unmixing <- solve( tcrossprod( design ), design )[ prep$panel, , drop = FALSE ]

  m <- diag( length( prep$panel ) ) +
    ( prep$s.true.panel - prep$s.start.panel ) %*% t( unmixing )
  dimnames( m ) <- list( prep$panel, prep$panel )

  d <- diag( m )
  d[ !is.finite( d ) | d == 0 ] <- 1
  m / d
}

#' Phase one: the residual spillover matrix, iterated to its best delta.
#' `m.init` lets a caller warm-start from a previous call's result rather
#' than the identity; the fixed point being solved for is the same either
#' way, so this only changes how many iterations it takes to get there. Stage
#' R uses it, since each pass's true matrix is close to the last pass's.
#' @noRd
.fx.estimate.spillover <- function( prep, abundance, unstained.abundance,
                                    max.iter = fx.max.iter, verbose = TRUE,
                                    m.init = NULL ) {

  panel <- prep$panel
  n.fl  <- length( panel )

  m.hat <- if ( !is.null( m.init ) ) m.init else diag( n.fl )
  dimnames( m.hat ) <- list( panel, panel )
  m.best <- m.hat; delta.best <- Inf; worst <- NA_character_
  history <- numeric( 0 ); final.log <- NULL
  
  # Warm-starts each pair's Huber fit from its own slope on the previous
  # inner iteration, the same change made to .fix.envelope.slope()'s caller
  # in fix_my_unmix.R -- the fixed point is unchanged, only how many IRLS
  # iterations it takes to reach it.
  slope.warm <- matrix( 0, n.fl, n.fl, dimnames = list( panel, panel ) )

  for ( iter in seq_len( max.iter ) ) {

    comp <- tryCatch( solve( m.hat ), error = function( e ) NULL )
    if ( is.null( comp ) ) break

    unmixed <- abundance           %*% comp
    unst    <- unstained.abundance %*% comp

    thr <- fx.threshold.margin *
      apply( unst, 2, stats::quantile, probs = fx.threshold.probs, names = FALSE )
    names( thr ) <- panel

    neg.var <- apply( unst, 2, stats::mad )^2
    names( neg.var ) <- panel

    thr.mat <- get.spread.thresholds(
      unmixed = unmixed, thresholds = thr,
      spillover.spread = prep$spillover.spread,
      spread.kappa = fx.spread.kappa, verbose = FALSE )

    slope <- diag( n.fl ); trust <- diag( n.fl )
    dimnames( slope ) <- dimnames( trust ) <- list( panel, panel )

    # Every (source, channel) pair reads only these already-built inputs and
    # writes nothing shared, so the pair loop is flattened into one task list
    # and handed to fx.lapply. Row order matches the original nested loop
    # (source outer, channel inner) so a log diffs the same way it always did.
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
        quantiles            = c( 0.05, 0.5 ),
        n.levels             = fx.n.levels.pair,
        min.events           = fx.min.events,
        max.truncated.events = fx.max.truncated,
        max.coefficient      = fx.max.coefficient,
        start.slope          = slope.warm[ src, ch ] )
    } )

    rows <- vector( "list", nrow( pairs ) )
    
    for ( i in seq_len( nrow( pairs ) ) ) {
      
      src <- pairs$src[ i ]; ch <- pairs$ch[ i ]
      est <- ests[[ i ]]
      
      if ( !is.null( est ) && is.finite( est$slope.truncated ) )
        slope.warm[ src, ch ] <- est$slope.truncated
      
      # Both estimates come back from the one call, so the robust fit is run
      # once per pair per iteration rather than twice.
      slope.use <- if ( is.null( est ) ) NA_real_ else
        if ( identical( fx.estimator, "truncated" ) ) est$slope.truncated else
          est$slope

      w <- 0
      if ( !is.null( est ) && is.finite( slope.use ) &&
           abs( slope.use ) <= fx.max.coefficient &&
           est$coverage >= 0.10 &&
           est$span > 5 * abs( thr[ src ] ) ) {
        w <- 1
        slope[ src, ch ] <- slope.use
      }
      trust[ src, ch ] <- w

      rows[[ i ]] <- data.frame(
        source = src, channel = ch,
        slope           = if ( is.null( est ) ) NA_real_ else est$slope,
        slope.truncated = if ( is.null( est ) ) NA_real_ else
          est$slope.truncated,
        slope.median = if ( is.null( est ) ) NA_real_ else est$slope.alt,
        se           = if ( is.null( est ) ) NA_real_ else est$se,
        disagreement = if ( is.null( est ) ) NA_real_ else est$disagreement,
        coverage     = if ( is.null( est ) ) NA_real_ else est$coverage,
        trust = w, row.names = NULL, stringsAsFactors = FALSE )
    }

    slope.error <- slope - diag( n.fl )
    delta.max   <- max( abs( slope.error ) )
    history     <- c( history, delta.max )

    if ( delta.max < delta.best ) {
      delta.best <- delta.max
      m.best     <- m.hat
      final.log  <- do.call( rbind, rows )
      hit        <- which( abs( slope.error ) == delta.max, arr.ind = TRUE )[ 1, ]
      worst      <- paste( panel[ hit[ 1 ] ], "->", panel[ hit[ 2 ] ] )
    }

    if ( verbose )
      cat( sprintf( "    iter %2d: delta.max %.5f%s\n", iter, delta.max,
                    if ( delta.max <= delta.best ) "   (best)" else "" ) )

    if ( delta.max < 0.005 ) break

    m.next <- m.hat + ( trust * slope.error ) %*% m.hat
    dn <- diag( m.next )
    if ( any( !is.finite( dn ) ) || any( dn <= 0 ) ) break
    m.hat <- sweep( m.next, 1, dn, "/" )
  }

  list( spillover = m.best, delta = delta.best, worst = worst,
        history = history, log = final.log )
}

#' Dominance assignment on an abundance matrix, the same score the signature
#' phase uses to build its populations. Each abundance is scored against its
#' own event's positivity boundary when a spread threshold matrix is
#' supplied, because a channel a bright dye spills into has a wider negative
#' population there and a flat cut donates that dye's mid-range events to
#' whichever narrower neighbour it spreads into.
#' @noRd
.fx3.dominant <- function( abundance, thresholds, threshold.matrix = NULL ) {

  flat <- thresholds[ colnames( abundance ) ]

  boundary <- if ( is.null( threshold.matrix ) )
    matrix( flat, nrow( abundance ), ncol( abundance ), byrow = TRUE ) else
      threshold.matrix[ , colnames( abundance ), drop = FALSE ]

  dyn <- pmax( apply( abundance, 2, stats::quantile, probs = 0.999,
                      names = FALSE ) - flat, .Machine$double.eps )

  score <- sweep( pmax( abundance - boundary, 0 ), 2, dyn, "/" )
  dom   <- max.col( score, ties.method = "first" )
  dom[ score[ cbind( seq_along( dom ), dom ) ] <= 0 ] <- 0L
  dom
}

# ---------------------------------------------------------------------------
# STAGE A - the real error, and the geometry it has to be measured in
# ---------------------------------------------------------------------------

cat( "\n===================== STAGE A: error structure =====================\n" )

fx.prep <- list(
  Cells = .fx.prepare( "Cells", "Beads" ),
  Beads = .fx.prepare( "Beads", "Cells" )
)

for ( tg in names( fx.prep ) ) {

  t.elapsed <- system.time( {

    p        <- fx.prep[[ tg ]]
    design   <- rbind( p$af.basis, p$s.start.panel )
    hot      <- calculate.hotspot.matrix( design )
    pair.hot <- hot[ p$panel, p$panel ]

    cat( sprintf( "\n-- %s raw data, starting from %s spectra --\n", p$target, p$start ) )
    cat( sprintf( "  gate.frac %.3f, events after stratified sample: %d\n",
                  p$gate.frac, nrow( p$raw ) ) )

    err <- p$s.true.panel - p$s.start.panel

    split        <- .signature.span.split( err, p$s.start.panel )
    split.design <- .signature.span.split( err, design )

    a.table <- data.frame(
      fluorophore    = p$panel,
      deg.start      = .fx.angle( p$s.start.panel, p$s.true.panel ),
      err.norm       = sqrt( rowSums( err^2 ) ),
      frac.in.span   = sqrt( rowSums( split$parallel^2 ) ) / sqrt( rowSums( err^2 ) ),
      frac.in.design = sqrt( rowSums( split.design$parallel^2 ) ) /
        sqrt( rowSums( err^2 ) ),
      row.names      = NULL )

    print( a.table, digits = 4 )
    cat( sprintf( "  panel: ||E|| %.4f, in-span fraction %.3f\n",
                  norm( err, "F" ),
                  norm( split$parallel, "F" ) / norm( err, "F" ) ) )

    a.table <- cbind( target = tg, a.table, stringsAsFactors = FALSE )

  } )

  fx.log.time( "A", tg, t.elapsed[ "elapsed" ] )

  if ( tg == names( fx.prep )[ 1 ] ) fx.csv.a <- list()
  fx.csv.a[[ tg ]] <- a.table
}

fx.write.csv( do.call( rbind, fx.csv.a ), "stage_A_error_structure.csv" )

cat( "\n  INTERPRETATION: frac.in.span is the share of the real error the\n" )
cat( "  spillover matrix can carry. What is left is what only a detector-space\n" )
cat( "  estimate can see, and is the entire justification for the hybrid.\n" )

# ---------------------------------------------------------------------------
# STAGE B - autofluorescence basis on the real unstained controls
# ---------------------------------------------------------------------------

cat( "\n===================== STAGE B: autofluorescence basis =====================\n" )

for ( tg in names( fx.prep ) ) {

  t.elapsed <- system.time( {

    p  <- fx.prep[[ tg ]]
    sv <- attr( p$af.basis, "singular.values" )

    cat( sprintf( "\n-- %s --\n", tg ) )
    cat( sprintf( "  components retained: %d\n", nrow( p$af.basis ) ) )
    cat( sprintf( "  singular values: %s\n",
                  paste( sprintf( "%.3g", sv ), collapse = " " ) ) )
    cat( sprintf( "  ratio to the next:  %s\n",
                  paste( sprintf( "%.2f", sv[ -length( sv ) ] / sv[ -1 ] ),
                         collapse = " " ) ) )
    permuted.threshold <- svd( apply( p$unstained, 2, sample ),
                               nu = 0L, nv = 0L )$d[ 2 ]
    cat( sprintf( "  permuted retention threshold: %.3g\n", permuted.threshold ) )

    hot <- calculate.hotspot.matrix( rbind( p$af.basis, p$s.start.panel ) )
    coupling <- apply( hot[ p$panel, rownames( p$af.basis ), drop = FALSE ], 1, max )
    cat( "  hotspot coupling to the AF basis, worst five:\n" )
    print( round( head( sort( coupling, decreasing = TRUE ), 5 ), 2 ) )

    b.table <- data.frame(
      target = tg, fluorophore = names( coupling ),
      hotspot.coupling = unname( coupling ),
      permuted.threshold = permuted.threshold,
      n.components = nrow( p$af.basis ),
      row.names = NULL, stringsAsFactors = FALSE )

  } )

  fx.log.time( "B", tg, t.elapsed[ "elapsed" ] )

  if ( tg == names( fx.prep )[ 1 ] ) fx.csv.b <- list()
  fx.csv.b[[ tg ]] <- b.table
}

fx.write.csv( do.call( rbind, fx.csv.b ), "stage_B_af_basis.csv" )

cat( "\n  INTERPRETATION: a ratio near 1.0 between consecutive singular values is\n" )
cat( "  the noise plateau; components inside it are being fitted to noise. Any\n" )
cat( "  dye with coupling above about 5 cannot be separated from the background\n" )
cat( "  and should be frozen rather than corrected.\n" )

# ---------------------------------------------------------------------------
# STAGE C - spillover estimation against the induced truth
# ---------------------------------------------------------------------------
# The target is the spillover the measured spectral difference actually
# induces under this design, not a nominal mixing matrix. This is now the
# same iterated estimator Stage L and Stage R use, run once from the starting
# spectra, so it reports what the mechanism itself achieves before any null
# debiasing rather than a separately-implemented single-pass fit.

cat( "\n===================== STAGE C: spillover coefficients =====================\n" )

fx.stage.c <- list()

for ( tg in names( fx.prep ) ) {

  t.elapsed <- system.time( {

    p <- fx.prep[[ tg ]]

    dec  <- deconvolve.af.background( p$raw,       p$s.start.panel, p$af.basis,
                                      af.name = NULL )
    decu <- deconvolve.af.background( p$unstained, p$s.start.panel, p$af.basis,
                                      af.name = NULL )

    m.induced <- .fx.induced.spillover( p )
    ph        <- .fx.estimate.spillover( p, dec$abundance, decu$abundance,
                                         verbose = FALSE )

    c.table <- data.frame(
      source    = ph$log$source,
      channel   = ph$log$channel,
      induced   = m.induced[ cbind( match( ph$log$source,  p$panel ),
                                    match( ph$log$channel, p$panel ) ) ],
      estimate  = ph$log$slope.truncated,
      coverage  = ph$log$coverage,
      trust     = ph$log$trust,
      row.names = NULL )

    c.table$err <- c.table$estimate - c.table$induced

    cat( sprintf( "\n-- %s --\n", tg ) )
    cat( "  worst eight pairs, iterated estimate vs induced truth:\n" )
    print( head( c.table[ order( -abs( c.table$err ) ), ], 8 ), digits = 3 )
    cat( sprintf(
      "  median |induced| %.4f, median |error| %.4f, recovered fraction %.3f\n",
      stats::median( abs( c.table$induced ) ),
      stats::median( abs( c.table$err ), na.rm = TRUE ),
      1 - stats::median( abs( c.table$err ), na.rm = TRUE ) /
        stats::median( abs( c.table$induced ) ) ) )

    fx.stage.c[[ tg ]] <- list( dec = dec, decu = decu, m.induced = m.induced,
                                table = c.table )

    c.table <- cbind( target = tg, c.table, stringsAsFactors = FALSE )

  } )

  fx.log.time( "C", tg, t.elapsed[ "elapsed" ] )

  if ( tg == names( fx.prep )[ 1 ] ) fx.csv.c <- list()
  fx.csv.c[[ tg ]] <- c.table
}

fx.write.csv( do.call( rbind, fx.csv.c ), "stage_C_spillover_vs_induced.csv" )

cat( "\n  READ FIRST: this is the estimator's own result before null debiasing,\n" )
cat( "  kept as a sanity check against ground truth on this substrate. Stage R\n" )
cat( "  is the one to trust in production, where there is no ground truth; if\n" )
cat( "  this table and Stage R disagree sharply in direction, look at Stage L\n" )
cat( "  first.\n" )

# ---------------------------------------------------------------------------
# STAGE L - null calibration
# ---------------------------------------------------------------------------
# Same particle type on both sides. The induced spillover is the identity and
# the true signature error is zero, so everything the estimator returns is
# its own bias. This run is available in production from the control data the
# reference spectra were extracted from, via `null.fit`.

cat( "\n===================== STAGE L: null calibration =====================\n" )

fx.null <- list(
  Cells = .fx.prepare( "Cells", "Cells" ),
  Beads = .fx.prepare( "Beads", "Beads" )
)

fx.stage.l <- list()

for ( tg in names( fx.null ) ) {

  t.elapsed <- system.time( {

    p <- fx.null[[ tg ]]

    cat( sprintf( "\n-- %s --\n", tg ) )

    dec  <- deconvolve.af.background( p$raw,       p$s.start.panel, p$af.basis,
                                      af.name = NULL )
    decu <- deconvolve.af.background( p$unstained, p$s.start.panel, p$af.basis,
                                      af.name = NULL )

    # The iterated phase-one estimator on its own control. Everything off the
    # diagonal is the artefact Stage R must subtract from the cross-run matrix.
    null.fit.idx <- if ( nrow( dec$abundance ) > fx.pass.events )
      sample.int( nrow( dec$abundance ), fx.pass.events ) else
        seq_len( nrow( dec$abundance ) )

    ph.null <- .fx.estimate.spillover(
      p, dec$abundance[ null.fit.idx, , drop = FALSE ], decu$abundance,
      verbose = FALSE )

    cat( "  largest coefficients invented on the null:\n" )
    ph.null$log$passes.gate <-
      abs( ph.null$log$slope.truncated ) <= fx.max.coefficient &
      ph.null$log$coverage >= 0.10
    print( head( ph.null$log[ order( -abs( ph.null$log$slope.truncated ) ),
                              c( "source", "channel", "slope.truncated",
                                 "coverage", "trust", "passes.gate" ) ], 8 ),
           digits = 3 )
    cat( sprintf( "  pairs above the coefficient cap: %d of %d\n",
                  sum( abs( ph.null$log$slope.truncated ) > fx.max.coefficient,
                      na.rm = TRUE ),
                  nrow( ph.null$log ) ) )

    # Signature phase on the null. Any angular change is pure estimator bias,
    # and its size per row is the floor below which a real correction cannot be
    # believed.

    thr <- fx.threshold.margin *
      apply( decu$abundance, 2, stats::quantile,
            probs = fx.threshold.probs, names = FALSE )
    names( thr ) <- p$panel

    thr.mat.l <- get.spread.thresholds(
      unmixed = dec$abundance, thresholds = thr,
      spillover.spread = p$spillover.spread,
      spread.kappa = fx.spread.kappa, verbose = FALSE )

    dom    <- .fx3.dominant( dec$abundance, thr, thr.mat.l )
    s.null <- p$s.start.panel
    n.rows <- list()

    for ( j in p$panel ) {

      idx <- which( dom == match( j, p$panel ) )
      if ( length( idx ) < fx.min.events ) next

      cand <- extract.raw.signature(
        raw.data       = dec$residual[ idx, , drop = FALSE ],
        spectra        = p$s.start.panel,
        abundance      = dec$abundance[ idx, , drop = FALSE ],
        target         = j,
        active         = p$panel,
        multivariate   = TRUE,
        ridge          = 1e-6,
        n.levels       = 60L,
        min.events     = fx.min.events,
        background.raw = if ( any( dom == 0L ) )
          dec$residual[ dom == 0L, , drop = FALSE ] else NULL )

      if ( is.null( cand ) ) next
      s.null[ j, ] <- cand$signature
      n.rows[[ length( n.rows ) + 1L ]] <- cand$stats
    }

    n.stats <- do.call( rbind, n.rows )

    n.table <- data.frame(
      fluorophore = p$panel,
      deg.bias    = .fx.angle( s.null, p$s.start.panel ),
      row.names   = NULL )

    m <- match( n.table$fluorophore, n.stats$fluorophore )
    n.table$explained       <- n.stats$explained[ m ]
    n.table$explained.total <- n.stats$explained.total[ m ]
    n.table$resid.rel     <- n.stats$resid.rel[ m ]
    n.table$intercept.rel <- n.stats$intercept.rel[ m ]
    n.table$clamp.frac    <- n.stats$clamp.frac[ m ]
    n.table$vif.target    <- n.stats$vif.target[ m ]

    cat( "\n  angular bias of the signature phase with nothing to correct:\n" )
    print( n.table, digits = 3 )

    fx.stage.l[[ tg ]] <- list( signature = n.table, s.null = s.null,
                                m.spillover = ph.null$spillover )

    l.coef.table <- cbind( target = tg, ph.null$log, stringsAsFactors = FALSE )
    l.bias.table <- cbind( target = tg, n.table, stringsAsFactors = FALSE )

  } )

  fx.log.time( "L", tg, t.elapsed[ "elapsed" ] )

  if ( tg == names( fx.null )[ 1 ] ) {
    fx.csv.l.coef <- list()
    fx.csv.l.bias <- list()
  }
  fx.csv.l.coef[[ tg ]] <- l.coef.table
  fx.csv.l.bias[[ tg ]] <- l.bias.table
}

fx.write.csv( do.call( rbind, fx.csv.l.coef ), "stage_L_null_coefficients.csv" )
fx.write.csv( do.call( rbind, fx.csv.l.bias ), "stage_L_null_signature_bias.csv" )

cat( "\n  READ FIRST: deg.bias is the floor. On a new dataset, watch whether it\n" )
cat( "  is small and stable relative to Stage A's deg.start; if it is not, the\n" )
cat( "  gate constants tuned on the Aurora panel (fx.impact.ratio, fx.step.decay,\n" )
cat( "  min.bin.events, max.mask.passes) are the ones to re-check first.\n" )

# ---------------------------------------------------------------------------
# STAGE R - alternating spillover and signature passes
# ---------------------------------------------------------------------------
# Each pass re-estimates the spillover matrix and the signatures against the
# spectra the previous pass produced, so the nuisance removal for a collinear
# dye sees its neighbour's corrected row rather than the row it started with.
# This is the cost centre and the current candidate algorithm.

cat( "\n===================== STAGE R: alternating passes =====================\n" )

fx.stage.r <- list()

for ( tg in names( fx.prep ) ) {

  t.elapsed.tg <- system.time( {

  p      <- fx.prep[[ tg ]]
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

  s.curr <- p$s.start.panel
  track  <- data.frame( fluorophore = p$panel,
                        pass0 = .fx.angle( s.curr, p$s.true.panel ),
                        row.names = NULL )

  ph <- NULL
  fx.csv.r.gate.tg <- list()

  for ( pass in seq_len( fx.n.pass ) ) {

    t.elapsed.pass <- system.time( {

    dec  <- deconvolve.af.background( p$raw,       s.curr, p$af.basis,
                                      af.name = NULL )
    decu <- deconvolve.af.background( p$unstained, s.curr, p$af.basis,
                                      af.name = NULL )

    fit.idx <- if ( nrow( dec$abundance ) > fx.pass.events )
      sample.int( nrow( dec$abundance ), fx.pass.events ) else
        seq_len( nrow( dec$abundance ) )

    # Warm-started from the previous pass's own result rather than the
    # identity: the true matrix moves only as much as this pass's accepted
    # rows moved it, so starting from last pass's answer converges the
    # inner iteration in a handful of steps instead of retracing the same
    # path from scratch every pass.
    ph <- .fx.estimate.spillover(
      p, dec$abundance[ fit.idx, , drop = FALSE ], decu$abundance,
      m.init = if ( pass > 1L ) ph$spillover else NULL,
      verbose = FALSE )

    unmixed <- dec$abundance %*% solve( ph$spillover )

    thr <- fx.threshold.margin *
      apply( decu$abundance %*% solve( ph$spillover ), 2, stats::quantile,
            probs = fx.threshold.probs, names = FALSE )
    names( thr ) <- p$panel

    thr.mat <- get.spread.thresholds(
      unmixed = unmixed, thresholds = thr,
      spillover.spread = p$spillover.spread,
      spread.kappa = fx.spread.kappa, verbose = FALSE )

    dom <- .fx3.dominant( unmixed, thr, thr.mat )

    # The unmixing operator of the current design, so a proposed step is
    # scored by the abundance error it would move rather than by its angle.
    # Its panel rows annihilate anything in the span of the autofluorescence
    # basis, so a step confined to that subspace scores zero, which is the
    # correct weight for a change that cannot reach the abundances.
    design.curr   <- rbind( p$af.basis, s.curr )
    unmixing.curr <- solve( tcrossprod( design.curr ),
                            design.curr )[ p$panel, , drop = FALSE ]

    x.high <- apply( unmixed, 2, stats::quantile, probs = 0.999, names = FALSE )
    names( x.high ) <- p$panel

    s.raw      <- s.curr
    accepted   <- character( 0 )
    gate.rows  <- list()
    step.total <- 0

    if ( pass == 1L ) {
      step.prev     <- stats::setNames( rep( Inf, length( p$panel ) ), p$panel )
      settled       <- character( 0 )
      ever.accepted <- character( 0 )
    }

    for ( j in p$panel ) {

      idx <- which( dom == match( j, p$panel ) )
      if ( length( idx ) < fx.min.events ) next

      cand <- extract.raw.signature(
        raw.data       = dec$residual[ idx, , drop = FALSE ],
        spectra        = s.curr,
        abundance      = unmixed[ idx, , drop = FALSE ],
        target         = j,
        active         = p$panel,
        multivariate   = TRUE,
        ridge          = 1e-6,
        n.levels       = 60L,
        min.events     = fx.min.events,
        background.raw = if ( any( dom == 0L ) )
          dec$residual[ dom == 0L, , drop = FALSE ] else NULL )

      if ( is.null( cand ) ) next

      st <- cand$stats
      if ( st$explained.total < 0.8 || st$explained.total > 1.2 ) next
      if ( st$resid.rel     > cut.resid )     next
      if ( st$intercept.rel > cut.intercept ) next
      step.impact <- unname( x.high[ j ] ) *
        sqrt( sum( ( ( cand$signature - s.curr[ j, ] ) %*%
                       t( unmixing.curr ) )^2 ) )

      bias.impact <- if ( j %in% shared )
        unname( x.high[ j ] ) *
        sqrt( sum( ( bias[ j, ] %*% t( unmixing.curr ) )^2 ) ) else 0

      drift <- .fx.angle( s.curr[ j, , drop = FALSE ],
                          p$s.start.panel[ j, , drop = FALSE ] )

      reason <- if ( !isTRUE( st$joint ) ) "univariate" else "accepted"

      if ( j %in% settled )                     reason <- "settled" else
        if ( st$explained.total < 0.8 || st$explained.total > 1.2 )
          reason <- "explained" else
            if ( st$resid.rel     > cut.resid )     reason <- "fit"     else
              if ( st$intercept.rel > cut.intercept ) reason <- "offset"  else
                if ( st$clamp.frac    > 0.15 )          reason <- "negative.mass" else
                  if ( step.impact < fx.impact.ratio * bias.impact )
                    reason <- "bias" else
                      if ( step.impact >= fx.step.decay * step.prev[ j ] )
                        reason <- "plateau" else
                          if ( drift > fx.total.max.angle )       reason <- "drift"

      gate.rows[[ length( gate.rows ) + 1L ]] <- data.frame(
        fluorophore = j, step.impact = step.impact,
        bias.impact = bias.impact,
        ratio       = step.impact / pmax( bias.impact, .Machine$double.eps ),
        step.prev   = unname( step.prev[ j ] ),
        joint       = isTRUE( st$joint ),
        reason      = reason,
        row.names   = NULL, stringsAsFactors = FALSE )

      # A row whose step has stopped shrinking has taken everything this
      # control can tell it. Every further pass swaps one draw of estimator
      # noise for another and the row wanders, so it is retired rather than
      # damped.
      if ( identical( reason, "plateau" ) ) settled <- c( settled, j )

      if ( ! reason %in% c( "accepted", "univariate" ) ) next

      step.prev[ j ] <- step.impact
      step.total     <- step.total + step.impact

      s.raw[ j, ] <- cand$signature
      accepted    <- c( accepted, j )
    }

    ever.accepted <- union( ever.accepted, accepted )

    # A row is frozen only if it has never once cleared the signature gate,
    # in this pass or any earlier one: no usable population, or a fit that
    # never passed explained/resid/intercept/clamp. A row that is merely
    # `settled` or `plateau` this pass converged in an earlier one and keeps
    # the spillover coefficients and last accepted signature it already
    # earned; `accepted` alone would freeze it too, on precisely the pass
    # where everything has converged and nothing is newly accepted, which
    # zeroes the whole matrix rather than the rows that actually never fit.
    spillover.pass <- ph$spillover

    if ( !is.null( m.null.bias ) ) {
      shared.m <- intersect( rownames( m.null.bias ), p$panel )
      spillover.pass[ shared.m, shared.m ] <-
        spillover.pass[ shared.m, shared.m ] -
        m.null.bias[ shared.m, shared.m, drop = FALSE ]
      diag( spillover.pass ) <- 1
    }

    frozen <- setdiff( p$panel, ever.accepted )

    if ( length( frozen ) > 0 ) {
      spillover.pass[ frozen, ] <- 0
      spillover.pass[ cbind( frozen, frozen ) ] <- 1
    }

    back.raw <- spillover.pass %*% s.curr

    # The split must use the same design the unmixing sees. Error in the span
    # of the AF basis shifts per-event AF coefficients and nothing else, so it
    # is left uncorrected on purpose, and only the bias component the perp
    # channel actually re-applies each pass is subtracted; the in-span bias is
    # handled in the matrix above.
    design.split <- rbind( p$af.basis, s.curr )

    perp      <- .signature.span.split( s.raw - s.curr,
                                        design.split )$perpendicular
    bias.perp <- .signature.span.split( bias, design.split )$perpendicular

    s.next <- back.raw + perp

    debias <- intersect( shared, accepted )

    if ( length( debias ) > 0 )
      s.next[ debias, ] <- s.next[ debias, , drop = FALSE ] -
      bias.perp[ debias, , drop = FALSE ]

    s.curr <- .fx.renorm( pmax( s.next, 0 ) )

    track[[ paste0( "pass", pass ) ]] <- .fx.angle( s.curr, p$s.true.panel )

    if ( pass == 1L ) {
      s.by.pass     <- list()
      delta.by.pass <- numeric( 0 )
      step.by.pass  <- numeric( 0 )
    }

    s.by.pass[[ pass ]]   <- s.curr
    delta.by.pass[ pass ] <- ph$delta
    step.by.pass[ pass ]  <- step.total

    cat( sprintf( "  %s pass %d: %d row(s) accepted, spillover delta %.5f\n",
                  tg, pass, length( accepted ), ph$delta ) )

    if ( length( gate.rows ) > 0 ) {
      cat( sprintf( "    gate detail, pass %d:\n", pass ) )
      gate.table <- do.call( rbind, gate.rows )
      print( gate.table, digits = 3 )
      fx.csv.r.gate.tg[[ pass ]] <- cbind(
        target = tg, pass = pass, gate.table, stringsAsFactors = FALSE )
    }

    } )

    fx.log.time( "R", sprintf( "%s pass %d", tg, pass ), t.elapsed.pass[ "elapsed" ] )
  }

  cat( sprintf( "\n-- %s --\n", tg ) )
  print( track, digits = 3 )
  design.q   <- rbind( p$af.basis, p$s.true.panel )
  unmixing.q <- solve( tcrossprod( design.q ), design.q )[ p$panel, ,
                                                           drop = FALSE ]

  dec.q <- deconvolve.af.background( p$raw, p$s.start.panel, p$af.basis,
                                     af.name = NULL )
  x.max <- apply( dec.q$abundance, 2, stats::quantile, probs = 0.999,
                  names = FALSE )

  impact <- function( s ) sum( x.max * sqrt( rowSums(
    ( ( s[ p$panel, , drop = FALSE ] - p$s.true.panel ) %*%
        t( unmixing.q ) )^2 ) ) )

  impact.start <- impact( p$s.start.panel )

  # The residual spillover the estimator still finds after a pass is the only
  # convergence signal that survives into production, where there is no
  # truth to score against. If it selects the same pass the impact selects,
  # the stopping rule transfers unchanged.
  pass.table <- data.frame(
    pass   = seq_along( s.by.pass ),
    delta  = delta.by.pass,
    step   = step.by.pass,
    angle  = vapply( s.by.pass, function( s )
      sum( .fx.angle( s, p$s.true.panel ) ), numeric( 1 ) ),
    impact = vapply( s.by.pass, impact, numeric( 1 ) ),
    row.names = NULL )

  pass.table$recovered <- 1 - pass.table$impact / impact.start

  best.delta  <- which.min( pass.table$delta )
  best.step   <- which.min( pass.table$step )
  best.impact <- which.min( pass.table$impact )

  cat( sprintf( "  start: angle %.2f deg, impact %.4g\n",
                sum( track$pass0 ), impact.start ) )
  print( pass.table, digits = 3 )
  cat( sprintf( paste0( "  smallest residual spillover at pass %d, smallest ",
                        "accepted step at pass %d, smallest true impact at ",
                        "pass %d\n" ),
                best.delta, best.step, best.impact ) )

  fx.stage.r[[ tg ]] <- list(
    s.final         = s.by.pass[[ length( s.by.pass ) ]],
    s.best          = s.by.pass[[ best.step ]],
    pass.table      = pass.table,
    x.max           = x.max,
    unmixing.q      = unmixing.q,
    spillover.final = spillover.pass,
    log.final       = ph$log )

  if ( tg == names( fx.prep )[ 1 ] ) {
    fx.csv.r.pass  <- list()
    fx.csv.r.track <- list()
    fx.csv.r.gate  <- list()
  }
  fx.csv.r.pass[[ tg ]]  <- cbind( target = tg, pass.table, stringsAsFactors = FALSE )
  fx.csv.r.track[[ tg ]] <- cbind( target = tg, track, stringsAsFactors = FALSE )
  fx.csv.r.gate[[ tg ]]  <- do.call( rbind, fx.csv.r.gate.tg )

  } )

  fx.log.time( "R", paste0( tg, " total" ), t.elapsed.tg[ "elapsed" ] )
}

fx.write.csv( do.call( rbind, fx.csv.r.pass ),  "stage_R_pass_table.csv" )
fx.write.csv( do.call( rbind, fx.csv.r.track ), "stage_R_angle_track.csv" )
fx.write.csv( do.call( rbind, fx.csv.r.gate ),  "stage_R_gate_detail.csv" )

cat( "\n  READ FIRST: the row to watch is PerCP-eFluor 710 on cells. If a second\n" )
cat( "  pass improves it once BV711's row has been corrected, the collinear\n" )
cat( "  nuisance removal was the limit and the alternating loop is the fix. If\n" )
cat( "  it plateaus at pass one, the limit is elsewhere and the co-expression\n" )
cat( "  machinery is the next place to look.\n" )

# ---------------------------------------------------------------------------
# STAGE T - the spillover row, channel by channel, for the top-impact dyes
# ---------------------------------------------------------------------------
# Stage R's final pass, one row at a time, against the induced truth, so a
# failure is attributed to named channels rather than to a dye.

cat( "\n===================== STAGE T: spillover row audit =====================\n" )

for ( tg in names( fx.prep ) ) {

  if ( is.null( fx.stage.r[[ tg ]] ) ) next

  t.elapsed <- system.time( {

  p         <- fx.prep[[ tg ]]
  r         <- fx.stage.r[[ tg ]]
  m.induced <- fx.stage.c[[ tg ]]$m.induced

  design.t   <- rbind( p$af.basis, p$s.true.panel )
  unmixing.t <- solve( tcrossprod( design.t ), design.t )[ p$panel, ,
                                                           drop = FALSE ]

  x.max <- apply( fx.stage.c[[ tg ]]$dec$abundance, 2, stats::quantile,
                  probs = 0.999, names = FALSE )

  impact <- x.max * sqrt( rowSums(
    ( ( p$s.start.panel - p$s.true.panel ) %*% t( unmixing.t ) )^2 ) )

  top <- p$panel[ order( -impact ) ][ seq_len( min( fx.audit.n,
                                                    length( p$panel ) ) ) ]

  cat( sprintf( "\n-- %s --\n", tg ) )

  fx.csv.t.tg <- list()

  for ( j in top ) {

    ch <- setdiff( p$panel, j )

    t.table <- data.frame(
      channel   = ch,
      induced   = m.induced[ j, ch ],
      estimate  = r$spillover.final[ j, ch ],
      row.names = NULL )

    t.table$err <- t.table$estimate - t.table$induced

    key <- match( paste( j, ch ),
                  paste( r$log.final$source, r$log.final$channel ) )
    t.table$fitted <- r$log.final$trust[ key ] > 0

    cat( sprintf(
      "\n  source %s: sum|induced| %.4f, sum|err| %.4f, cos %.3f\n", j,
      sum( abs( t.table$induced ) ), sum( abs( t.table$err ) ),
      sum( t.table$induced * t.table$estimate ) /
        max( sqrt( sum( t.table$induced^2 ) ) *
               sqrt( sum( t.table$estimate^2 ) ), .Machine$double.eps ) ) )

    print( head( t.table[ order( -abs( t.table$err ) ), ], 6 ), digits = 3 )

    fx.csv.t.tg[[ j ]] <- cbind( target = tg, source = j, t.table,
                                 stringsAsFactors = FALSE )
  }

  } )

  fx.log.time( "T", tg, t.elapsed[ "elapsed" ] )

  if ( tg == names( fx.prep )[ 1 ] ) fx.csv.t <- list()
  fx.csv.t[[ tg ]] <- do.call( rbind, fx.csv.t.tg )
}

fx.write.csv( do.call( rbind, fx.csv.t ), "stage_T_row_audit.csv" )

cat( "\n  READ FIRST: a cosine near zero or negative with sum|err| comparable to\n" )
cat( "  sum|induced| means the row is being estimated from pairs that carry no\n" )
cat( "  usable information, and the row should be frozen rather than corrected.\n" )
cat( "  If the error concentrates on a few named channels, those pairs are the\n" )
cat( "  identifiability problem and the fix belongs in the pair estimator.\n" )

# ---------------------------------------------------------------------------
# STAGE U - the anatomy of what is left
# ---------------------------------------------------------------------------
# The remaining error of each row, against Stage R's final pass, split three
# ways against the final design: the part inside the panel span, which the
# back-solve can still reach; the part inside the autofluorescence span but
# orthogonal to the panel, which shifts per-event AF coefficients and nothing
# else and contributes exactly zero to the panel impact; and the part outside
# the design entirely, which only the perpendicular channel can reach. The
# impact column is the honest headroom.

cat( "\n===================== STAGE U: what is left =====================\n" )

for ( tg in names( fx.stage.r ) ) {

  t.elapsed <- system.time( {

  p       <- fx.prep[[ tg ]]
  s.final <- fx.stage.r[[ tg ]]$s.final
  x.max   <- fx.stage.r[[ tg ]]$x.max
  w.panel <- fx.stage.r[[ tg ]]$unmixing.q

  # The split must be oblique, not orthogonal. The unmixing operator resolves
  # a spectrum into autofluorescence and panel coefficients that are not
  # orthogonal to each other, and it is those coefficients the abundances
  # inherit. An orthogonal split cuts the error in a place the operator does
  # not recognise, and the pieces then carry impacts that do not sum to the
  # whole.
  design.u <- rbind( p$af.basis, p$s.true.panel )
  resolve  <- t( solve( tcrossprod( design.u ), design.u ) )
  n.af     <- nrow( p$af.basis )

  err  <- s.final[ p$panel, , drop = FALSE ] - p$s.true.panel
  coef <- err %*% resolve

  part.af    <- coef[ ,  seq_len( n.af ), drop = FALSE ] %*% p$af.basis
  part.panel <- coef[ , -seq_len( n.af ), drop = FALSE ] %*% p$s.true.panel
  part.out   <- err - part.af - part.panel

  norm.row <- function( m ) sqrt( rowSums( m^2 ) )

  err.start  <- p$s.start.panel - p$s.true.panel
  coef.start <- err.start %*% resolve

  u.table <- data.frame(
    fluorophore = p$panel,
    err.start   = norm.row( err.start ),
    err.final   = norm.row( err ),
    af.part     = norm.row( part.af ),
    panel.part  = norm.row( part.panel ),
    out.part    = norm.row( part.out ),
    row.names   = NULL )

  # Only the panel coefficients reach the abundances. Error resolved onto the
  # autofluorescence basis moves per-event AF coefficients and nothing else,
  # and error the design cannot represent lands in the residual, so both are
  # invisible to the panel to first order.
  u.table$impact.start <- x.max *
    norm.row( coef.start[ , -seq_len( n.af ), drop = FALSE ] )
  u.table$impact.final <- x.max *
    norm.row( coef[ , -seq_len( n.af ), drop = FALSE ] )
  u.table$recovered <- 1 - u.table$impact.final / u.table$impact.start

  cat( sprintf( "\n-- %s --\n", tg ) )
  print( u.table[ order( -u.table$impact.final ), ], digits = 3 )

  cat( sprintf( paste0( "  remaining spectral error %.4g of %.4g; of that, AF ",
                        "%.4g, panel %.4g, outside the design %.4g\n" ),
                sum( u.table$err.final ), sum( u.table$err.start ),
                sum( u.table$af.part ), sum( u.table$panel.part ),
                sum( u.table$out.part ) ) )
  cat( sprintf( "  remaining panel impact %.4g of %.4g, recovered %.3f\n",
                sum( u.table$impact.final ), sum( u.table$impact.start ),
                1 - sum( u.table$impact.final ) /
                  sum( u.table$impact.start ) ) )

  u.table <- cbind( target = tg, u.table, stringsAsFactors = FALSE )

  } )

  fx.log.time( "U", tg, t.elapsed[ "elapsed" ] )

  if ( tg == names( fx.stage.r )[ 1 ] ) fx.csv.u <- list()
  fx.csv.u[[ tg ]] <- u.table
}

fx.write.csv( do.call( rbind, fx.csv.u ), "stage_U_remaining_error.csv" )

cat( "\n  READ FIRST: `panel.part` is the only component the abundances see. A\n" )
cat( "  row with a large `err.final` but a small `panel.part` is spectrally\n" )
cat( "  wrong and numerically harmless and should not be chased. A row where\n" )
cat( "  `panel.part` is still a large share of `err.final` is where the\n" )
cat( "  spillover phase has signal left on the table.\n" )

# ---------------------------------------------------------------------------
# STAGE F (optional) - spectral traces against ground truth
# ---------------------------------------------------------------------------
# Off by default; set fx.write.traces <- TRUE to render. Simplified to start /
# final / truth, since there are no separate raw and back arms any more.

if ( fx.write.traces ) {

  cat( "\n===================== STAGE F: spectral traces =====================\n" )

  if ( !dir.exists( trace.dir ) ) dir.create( trace.dir, recursive = TRUE )

  for ( tg in names( fx.prep ) ) {

    p       <- fx.prep[[ tg ]]
    s.final <- fx.stage.r[[ tg ]]$s.final

    for ( j in p$panel ) {

      cmp <- rbind( p$s.start.panel[ j, ], s.final[ j, ], p$s.true.panel[ j, ] )
      colnames( cmp ) <- p$detectors
      rownames( cmp ) <- c(
        sprintf( "start (%.2f deg)", .fx.angle( p$s.start.panel[ j, , drop = FALSE ],
                                                 p$s.true.panel[ j, , drop = FALSE ] ) ),
        sprintf( "final (%.2f deg)", .fx.angle( s.final[ j, , drop = FALSE ],
                                                 p$s.true.panel[ j, , drop = FALSE ] ) ),
        "truth" )

      spectral.trace(
        spectral.matrix = cmp, asp = asp,
        title    = paste0( tolower( tg ), "_signature_", make.names( j ) ),
        plot.dir = trace.dir, split.lasers = FALSE,
        color.palette = "viridis", save = TRUE )
    }
  }

  cat( sprintf( "  wrote traces to %s\n", trace.dir ) )
}

# ---------------------------------------------------------------------------
# Timing summary
# ---------------------------------------------------------------------------

fx.write.csv( fx.timing, "stage_timing.csv" )
cat( sprintf( "\nAll stage tables and %s written to %s\n",
              "stage_timing.csv", fx.output.dir ) )
print( fx.timing, digits = 3 )
