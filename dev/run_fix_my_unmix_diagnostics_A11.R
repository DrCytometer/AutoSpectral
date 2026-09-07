# run_fix_my_unmix_diagnostics_A11.R
#
# Consolidated re-run of the bead / cell benchmark with every correction from
# the A10 round applied, in one script, so a single execution produces every
# result needed for the next determination. Replaces re-running A8 and A9.
#
# What A10 established, and what this script does about it:
#
#   1. The A8 harness never passed `source.mask` into the pair estimator,
#      while production fix.my.unmix() passes it by default. With the mask,
#      cor(estimate, induced) on Cells goes from -0.21 to +0.74 and every
#      magnitude band recovers 0.87-0.94 (stage_AC). Every A8/A9 cell-side
#      result was therefore measured on an estimator production does not run.
#      -> The estimator here computes dominance every inner iteration and
#         passes the mask (fx.source.dominant).
#
#   2. With the mask on, the largest true coefficients come back at 0.24-0.35,
#      above fx.max.coefficient = 0.2, and a pair refused by that cap is
#      invisible to both the matrix update and the convergence criterion
#      (its slope entry stays at the identity), so it can never be corrected.
#      -> fx.max.coefficient = 0.5, with the count of pairs above the cap
#         logged every run.
#
#   3. The phase-two gate refused exactly the rows with the largest errors,
#      because resid.rel / intercept.rel / explained.total have a per-dye
#      floor the null run measures and the cascade compared raw values
#      against one panel-median cut. Stage AD: the per-row rule releases the
#      four non-recovering cell dyes (BUV395, BV510, BV570, PerCP-Cy5.5) and,
#      with the floor at the full panel cut, newly refuses nothing on either
#      substrate.
#      -> Stage R gates on per-row null-referenced cuts,
#         cut_j = max( fx.row.ratio * null_j, panel cut ).
#
#   4. `plateau` / `settled` rows were retired from the signature phase but
#      not from the spillover fold-back, so converged rows kept drifting on
#      estimator noise (AF647 on Cells: 4.3 -> 12.7 degrees after retirement).
#      -> Retired rows are added to the frozen set.
#
#   5. `coverage` fails open (returns 1 when no bright events survive).
#      -> Every gate here tests with isTRUE(), so a future NA-returning
#         coverage and the current fail-open value are both handled.
#
# The five information sources of CONTEXT_information_sources.md, adjudicated
# by A9/A10 and carried here:
#
#   Idea 5 (local topology). Mode trace refuted on both substrates (26-30%
#     reduction on refused pairs, below the 40% bar; worse than the mask on
#     fitted pairs). The within-cluster shear PASSED at high cluster count:
#     n=1000/min=50 gave a 72% error reduction on refused-and-real cell pairs
#     (stage_AE). But that refusal set came from the mask-off estimator, so
#     STAGE Z here re-scores the shear against the corrected estimator's own
#     refusals, at n = 1000 and 2000. This is section 10 item 6's condition,
#     re-run against the right baseline.
#   Idea 2 (structural completion). frac.oracle passed (0.71 / 0.94), blanket
#     licence failed (loo.mae.pred > loo.mae.base), donor axis failed
#     (frac.donor == frac.wrong except BUV805). Verdict stands: no database
#     change; completion only with per-row licensing; not re-run here because
#     its inputs (variant deltas) are mask-independent. Section 10 item 5: no.
#   Idea 1 (hypernegativity). Sign-dependent trust term refuted (asymmetry
#     inconsistent across substrates; upper tail the stronger predictor on
#     Cells). The symmetric tail-excess DID detect damaging-unfitted pairs
#     (AUC 0.75 / 0.87). STAGE W re-scores that detector against the
#     corrected estimator's refusals. Section 10 item 7: gate rescue only,
#     no trust term.
#   Idea 3 (spread as variance / VIF as trust). Refuted as weighting (rms
#     worse on both substrates); sd(z) 8.8 vs 43 across substrates so no
#     analytic error works; the null coefficient spread fallback also failed
#     (null k = 0.18 / 0.00 vs real 6.3 / 29.5). Section 10 item 4: no.
#     Not re-run.
#   Idea 4 (SSI). spearman(ssi, |induced|) = 0.07; `after` nowhere near
#     `truth`; 31% false-positive floor on beads. Section 10 item 8: not
#     usable as the scoreboard on this data. Not re-run.
#
# Two further questions this script answers:
#
#   STAGE G. The bead-versus-cell gap. The concrete substrate differences are
#     brightness and scatter, so Stage G audits the scatter gate (kept
#     fraction at three gate levels, and its stability) and the positivity
#     thresholds (flat vs spread-scaled inflation per dye, and each dye's
#     surviving bright population), per substrate. If the bead main
#     population is being cut by the scatter gate, or the spread term is
#     inflating bead thresholds until the positive populations vanish, it
#     shows up here as numbers rather than a guess.
#
#   STAGE L. The null run's validity (the ruler problem). The null "bias"
#     is e_est - e_lib: the difference between the estimator's artefact and
#     the library's own error, and subtracting it assumes the library is
#     truth. Stage L therefore adds a ground-truth-free held-out check: fit
#     the null signature on half the events, then ask whether the library or
#     the null-corrected signature unmixes the OTHER half more cleanly
#     (detector-space reconstruction residual, and off-target abundance mass
#     per dominance population). Neither metric references an assumed-true
#     spectrum. If the null-corrected signature wins on unseen events, the
#     "bias" contains real library error and Stage R's subtraction is
#     discarding signal; if it loses or ties, the subtraction stands. Both
#     split directions are run and averaged.
#
# PREREQUISITES. Same session inputs as A8 (a8.cell.spectra, a8.bd.spectra,
# a8.concat.cells, a8.concat.beads, a8.cell.dir, a8.bd.dir, a8.cell.variants,
# a8.bead.variants) and the AutoSpectral package loaded via
# devtools::load_all(). NO package edits are required to run this script;
# the production-side edits remain a separate step after these results.
#
# Run top to bottom. Report back the CSVs listed at the end.

# ---------------------------------------------------------------------------
# 0. Setup - EDIT THIS SECTION
# ---------------------------------------------------------------------------

asp <- get.autospectral.param( cytometer = "discover" )
asp$scatter.data.max.x <- 5e7
asp$scatter.data.max.y <- 5e7

fx.spectra <- list(
  Cells = a8.cell.spectra,
  Beads = a8.bd.spectra
)

fx.raw <- list(
  Cells = a8.concat.cells,
  Beads = a8.concat.beads
)

fx.unstained <- list(
  Cells = readFCS( file.path( a8.cell.dir, "H01_Negative_001.fcs" ) ),
  Beads = readFCS( file.path( a8.bd.dir,   "H01_Negative_001.fcs" ) )
)

fx.variants <- list(
  Cells = a8.cell.variants,
  Beads = a8.bead.variants
)

af.name <- "AF"

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------

fx.threshold.probs  <- 0.995
fx.threshold.margin <- 1.3
fx.spread.kappa     <- 2
fx.seed             <- 42L
fx.max.iter         <- 15L
fx.n.levels.pair    <- 10L
fx.min.events       <- 100L
fx.estimator        <- "truncated"

# Raised from 0.2: with the dominance mask the largest true coefficients on
# this benchmark are 0.24-0.35, and a pair refused by the cap contributes
# nothing to the update or the convergence criterion, so it stays wrong
# silently. The count of pairs above the cap is logged so a future dataset
# that needs more room announces itself.
fx.max.coefficient  <- 0.5

# Stage C's own cap comparison found 20,000 beats 400,000 on both substrates
# once the dominance mask is on - decisively on Beads, where 400,000 gave a
# negative correlation with the induced truth. Matches production's
# max.truncated.events default; the harness was carrying a stale, mask-off
# value.
fx.max.truncated    <- 20000L

# Restrict a source's bright events to those the source is dominant for, as
# production fix.my.unmix( source.dominant = TRUE ) does. Off, the estimator's
# correlation with the induced truth on this panel was -0.21; on, +0.74.
fx.source.dominant  <- TRUE

# Per-row null-referenced gate cuts: cut_j = max( ratio * null_j, panel cut ).
# With the floor at the full panel cut the rule is a pure loosening; verified
# on this dataset to release the four non-recovering cell dyes and newly
# refuse nothing on either substrate. The explained gate becomes a band
# around each row's own null value.
fx.row.ratio           <- 3
fx.max.explained.drift <- 0.2

fx.max.events       <- 400000L
fx.background.frac  <- 0.3
fx.min.stratum      <- 3000L

fx.n.pass           <- 3L
fx.pass.events      <- 400000L
fx.impact.ratio     <- 1.0
fx.step.decay       <- 0.95
fx.total.max.angle  <- 15

fx.min.span         <- 5
fx.min.bg.align     <- -0.9
fx.max.angle        <- 16
fx.max.anchor       <- 0.10
fx.max.vif           <- 500
fx.peak.shift.min.rel <- 0.7
fx.max.hotspot       <- 5

# Stage G scatter gate audit levels; the middle one is what .fx.prepare uses.
fx.gate.levels      <- c( 0.05, 0.10, 0.20 )

# Stage Z and W event subsample and cluster settings. n = 1000 was the best
# arm in stage_AE and had not plateaued, so 2000 brackets it from above.
fx.stat.events      <- 100000L
fx.cofactor         <- 500
fx.cluster.grid     <- c( 1000L, 2000L )
fx.min.cluster      <- 50L

# Stage W tail-excess quantiles, as A9 Stage W.
fx.hyper.q          <- 0.02
fx.upper.q          <- 0.98
fx.bright.frac      <- 1 / 3

# ---------------------------------------------------------------------------
# Parallel backend
# ---------------------------------------------------------------------------

fx.parallel     <- TRUE
fx.mc.cores.max <- 4L
fx.mc.cores     <- if ( fx.parallel && .Platform$OS.type == "unix" &&
                        parallel::detectCores() > 1L )
  max( 1L, min( fx.mc.cores.max, parallel::detectCores() - 1L ) ) else 1L

#' @noRd
fx.lapply <- function( x, fun, ... ) {

  if ( fx.mc.cores <= 1L ) return( lapply( x, fun, ... ) )

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

fx.output.dir <- "./figure_fix_my_unmix_output_A11"
if ( dir.exists( fx.output.dir ) ) unlink( fx.output.dir, recursive = TRUE )
dir.create( fx.output.dir, recursive = TRUE )

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
# Helpers (carried from A8 unchanged unless noted)
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

#' @noRd
.fx.stratified.sample <- function( abundance, thresholds, threshold.matrix,
                                   n.total, background.frac = 0.3,
                                   min.stratum = 2000L ) {

  n <- nrow( abundance )
  if ( n <= n.total ) return( seq_len( n ) )

  positive <- ( abundance - threshold.matrix ) > 0

  background.idx <- which( rowSums( positive ) == 0L )
  positive.idx   <- lapply( seq_len( ncol( abundance ) ),
                            function( j ) which( positive[ , j ] ) )

  n.background <- min( length( background.idx ),
                       round( background.frac * n.total ) )
  budget.positive <- n.total - n.background

  sizes   <- vapply( positive.idx, length, integer( 1 ) )
  floor.n <- pmin( sizes, as.integer( min.stratum ) )

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
        gate.frac = mean( gate.keep ), unst.gate.frac = mean( unst.keep ),
        spillover.spread = fx.variants[[ target ]]$spillover.spread )
}

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

#' Phase one, corrected: dominance is assigned every inner iteration and
#' passed as `source.mask`, matching production fix.my.unmix() under
#' `source.dominant = TRUE`; the trust gate tests coverage with isTRUE() so
#' a pair with no surviving bright population cannot pass on a fail-open
#' default; and the count of coefficients above the cap is logged.
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

    dom.iter <- if ( isTRUE( fx.source.dominant ) )
      .fx3.dominant( unmixed, thr, thr.mat ) else NULL

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
        source.mask          = if ( is.null( dom.iter ) ) NULL else
          dom.iter == match( src, panel ),
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

      slope.use <- if ( is.null( est ) ) NA_real_ else
        if ( identical( fx.estimator, "truncated" ) ) est$slope.truncated else
          est$slope

      w <- 0
      if ( !is.null( est ) && is.finite( slope.use ) &&
           abs( slope.use ) <= fx.max.coefficient &&
           isTRUE( est$coverage >= 0.10 ) &&
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

    if ( verbose ) {
      n.capped <- sum( vapply( ests, function( e )
        !is.null( e ) && is.finite( e$slope.truncated ) &&
          abs( e$slope.truncated ) > fx.max.coefficient, logical( 1 ) ) )
      cat( sprintf( "    iter %2d: delta.max %.5f, %d pair(s) above the cap%s\n",
                    iter, delta.max, n.capped,
                    if ( delta.max <= delta.best ) "   (best)" else "" ) )
    }

    if ( delta.max < 0.005 ) break

    m.next <- m.hat + ( trust * slope.error ) %*% m.hat
    dn <- diag( m.next )
    if ( any( !is.finite( dn ) ) || any( dn <= 0 ) ) break
    m.hat <- sweep( m.next, 1, dn, "/" )
  }

  list( spillover = m.best, delta = delta.best, worst = worst,
        history = history, log = final.log )
}

#' Signature-phase null fit on one raw pool: unmix against the pool's own
#' library, take dominance populations, and re-derive each row. Returns the
#' re-derived spectra and the per-row statistics whose values are the gate
#' floors. Factored out of Stage L so the held-out validity check can run it
#' on half-pools without duplicating the body.
#' @noRd
.fx11.null.signature <- function( p, raw.pool ) {

  dec  <- deconvolve.af.background( raw.pool,    p$s.start.panel, p$af.basis,
                                    af.name = NULL )
  decu <- deconvolve.af.background( p$unstained, p$s.start.panel, p$af.basis,
                                    af.name = NULL )

  thr <- fx.threshold.margin *
    apply( decu$abundance, 2, stats::quantile,
           probs = fx.threshold.probs, names = FALSE )
  names( thr ) <- p$panel

  thr.mat <- get.spread.thresholds(
    unmixed = dec$abundance, thresholds = thr,
    spillover.spread = p$spillover.spread,
    spread.kappa = fx.spread.kappa, verbose = FALSE )

  dom    <- .fx3.dominant( dec$abundance, thr, thr.mat )
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

  n.stats <- if ( length( n.rows ) > 0 ) do.call( rbind, n.rows ) else NULL

  n.table <- data.frame(
    fluorophore = p$panel,
    deg.bias    = .fx.angle( s.null, p$s.start.panel ),
    row.names   = NULL )

  if ( !is.null( n.stats ) ) {
    m <- match( n.table$fluorophore, n.stats$fluorophore )
    n.table$explained       <- n.stats$explained[ m ]
    n.table$explained.total <- n.stats$explained.total[ m ]
    n.table$resid.rel       <- n.stats$resid.rel[ m ]
    n.table$intercept.rel   <- n.stats$intercept.rel[ m ]
    n.table$clamp.frac      <- n.stats$clamp.frac[ m ]
    n.table$vif.target      <- n.stats$vif.target[ m ]
  }

  list( s.null = s.null, table = n.table )
}

#' Held-out fit quality for a signature matrix, with no assumed-true spectrum
#' anywhere in it. `resid` is the mean detector-space reconstruction error on
#' events the signature was not fitted to; `leak` is the off-target abundance
#' mass carried by each dye's own dominance population, relative to its own
#' median abundance.
#' @noRd
.fx11.holdout.quality <- function( raw.holdout, p, spectra ) {

  design <- rbind( p$af.basis, spectra[ p$panel, , drop = FALSE ] )

  x <- unmix.ols.fast( raw.holdout,  design )
  u <- unmix.ols.fast( p$unstained,  design )
  colnames( x ) <- rownames( design )
  colnames( u ) <- rownames( design )

  resid <- sqrt( rowSums( ( raw.holdout - x %*% design )^2 ) )

  xp <- x[ , p$panel, drop = FALSE ]
  up <- u[ , p$panel, drop = FALSE ]

  thr <- fx.threshold.margin *
    apply( up, 2, stats::quantile, probs = fx.threshold.probs, names = FALSE )
  names( thr ) <- p$panel

  thr.mat <- get.spread.thresholds(
    unmixed = xp, thresholds = thr,
    spillover.spread = p$spillover.spread,
    spread.kappa = fx.spread.kappa, verbose = FALSE )

  dom <- .fx3.dominant( xp, thr, thr.mat )

  leak <- vapply( seq_along( p$panel ), function( k ) {

    idx <- which( dom == k )
    if ( length( idx ) < fx.min.events ) return( NA_real_ )

    own <- stats::median( xp[ idx, k ] )
    if ( !is.finite( own ) || own <= 0 ) return( NA_real_ )

    sum( abs( apply( xp[ idx, -k, drop = FALSE ], 2, stats::median ) ) ) / own
  }, numeric( 1 ) )

  names( leak ) <- p$panel

  list( resid = mean( resid ), leak = leak )
}

#' Huber slope with the residual scale taken from a nominated subset. With
#' `scale.index` at the source-bright events, the negative bulk still pins
#' the intercept but no longer sets the outlier scale that decides how hard
#' the signal-carrying population is downweighted.
#' @noRd
.fx11.huber.slope <- function( x, y, k = 1.345, max.iter = 100L, tol = 1e-4,
                               scale.index = NULL ) {

  n <- length( x )
  if ( n < 3L ) return( NA_real_ )

  if ( is.null( scale.index ) || length( scale.index ) < 3L )
    scale.index <- seq_len( n )

  w    <- rep( 1, n )
  coef <- c( 0, 0 )

  for ( iter in seq_len( max.iter ) ) {

    x.mean <- sum( w * x ) / sum( w )
    y.mean <- sum( w * y ) / sum( w )
    sxx    <- sum( w * ( x - x.mean )^2 )

    if ( !is.finite( sxx ) || sxx <= 0 ) return( NA_real_ )

    slope     <- sum( w * ( x - x.mean ) * ( y - y.mean ) ) / sxx
    intercept <- y.mean - slope * x.mean
    coef.next <- c( intercept, slope )

    resid <- y - ( intercept + slope * x )

    ref   <- resid[ scale.index ]
    scale <- stats::median( abs( ref - stats::median( ref ) ) ) / 0.6745
    scale <- max( scale, .Machine$double.eps )

    w <- pmin( 1, k / pmax( abs( resid / scale ), .Machine$double.eps ) )

    moved <- max( abs( coef.next - coef ) ) < tol * max( abs( coef.next ), 1 )
    coef  <- coef.next

    if ( moved ) break
  }

  unname( coef[ 2 ] )
}

#' Micro-clusters on an asinh-stabilised copy of the abundances.
#' @noRd
.fx11.cluster <- function( u, n.clusters, min.size = 50L, seed = 42L ) {

  space <- asinh( u / fx.cofactor )
  set.seed( seed )

  km <- tryCatch(
    stats::kmeans( space, centers = min( n.clusters, nrow( space ) - 1L ),
                   iter.max = 50L, algorithm = "MacQueen" ),
    error = function( e ) NULL )

  if ( is.null( km ) ) return( rep( NA_integer_, nrow( u ) ) )

  id   <- as.integer( km$cluster )
  size <- table( id )
  id[ id %in% as.integer( names( size )[ size < min.size ] ) ] <- NA_integer_
  id
}

#' Within-cluster shear: one slope, cluster-specific intercepts, restricted
#' to source-positive events. Between-cluster variation is absorbed by the
#' intercepts, so only within-cluster covariance identifies the slope.
#' @noRd
.fx11.cluster.shear <- function( x.source, y.target, threshold.source,
                                 cluster.id, min.cluster = 50L ) {

  empty <- c( m.within = NA_real_, m.floor = NA_real_, n.clusters = 0 )

  index <- which( x.source > threshold.source & !is.na( cluster.id ) )
  if ( length( index ) < min.cluster ) return( empty )

  x <- x.source[ index ]
  y <- y.target[ index ]
  g <- cluster.id[ index ]

  size <- table( g )
  use  <- names( size )[ size >= min.cluster ]
  if ( length( use ) == 0L ) return( empty )

  num <- 0
  den <- 0
  per <- rep( NA_real_, length( use ) )

  for ( i in seq_along( use ) ) {

    take <- g == use[ i ]
    xc   <- x[ take ] - mean( x[ take ] )
    yc   <- y[ take ] - mean( y[ take ] )

    num <- num + sum( xc * yc )
    den <- den + sum( xc^2 )

    if ( sum( xc^2 ) > 0 )
      per[ i ] <- .fx11.huber.slope( x[ take ], y[ take ] )
  }

  ok <- is.finite( per )

  c( m.within   = if ( den > 0 ) num / den else NA_real_,
     m.floor    = if ( any( ok ) ) min( per[ ok ] ) else NA_real_,
     n.clusters = sum( ok ) )
}

#' Every ordered pair of a panel.
#' @noRd
.fx11.pairs <- function( panel ) {
  data.frame(
    source  = rep( panel, each = length( panel ) - 1L ),
    channel = unlist( lapply( panel, function( s ) setdiff( panel, s ) ) ),
    stringsAsFactors = FALSE )
}

# ---------------------------------------------------------------------------
# Preparation
# ---------------------------------------------------------------------------

cat( "\n===================== A11: preparation =====================\n" )

fx.prep <- list(
  Cells = .fx.prepare( "Cells", "Beads" ),
  Beads = .fx.prepare( "Beads", "Cells" )
)

for ( tg in names( fx.prep ) )
  cat( sprintf( "  %s: gate.frac %.3f (unstained %.3f), events %d\n",
                tg, fx.prep[[ tg ]]$gate.frac, fx.prep[[ tg ]]$unst.gate.frac,
                nrow( fx.prep[[ tg ]]$raw ) ) )

# ---------------------------------------------------------------------------
# STAGE G - scatter gate and positivity threshold audit
# ---------------------------------------------------------------------------
# The two concrete bead-versus-cell differences are brightness and scatter.
# The gate half asks whether the density gate is stable on beads: the kept
# fraction at three gate levels, on stained and unstained files. A main
# population should hold a similar fraction across levels; a kept fraction
# that swings hard with the level, or sits far below the cell file's, means
# the gate is cutting into the population and every downstream population
# statistic inherits that.
#
# The threshold half asks whether bead brightness inflates the spread-scaled
# thresholds until the positive populations vanish: per dye, the flat
# threshold, the median spread inflation over its own dominant events, the
# dominant population size, and how many of its dominant events sit above
# each boundary. This is where "beads are brighter" becomes measurable: the
# spread term is a sum over bright abundances, so a substrate where every
# tube is bright to begin with raises every boundary at once.

cat( "\n===================== STAGE G: gating and thresholds =====================\n" )

fx.csv.g.gate <- list()
fx.csv.g.thr  <- list()

for ( tg in names( fx.prep ) ) {

  t.elapsed <- system.time( {

    p <- fx.prep[[ tg ]]

    scatter.det <- grep( "^(FSC|SSC)", colnames( fx.raw[[ tg ]] ), value = TRUE )

    gate.rows <- lapply( fx.gate.levels, function( lev ) {
      set.seed( fx.seed )
      keep.raw <- .fx.main.gate(
        fx.raw[[ tg ]][ , scatter.det, drop = FALSE ], gate.level = lev )
      set.seed( fx.seed )
      keep.unst <- .fx.main.gate(
        fx.unstained[[ tg ]][ , scatter.det, drop = FALSE ], gate.level = lev )
      data.frame( target = tg, gate.level = lev,
                  kept.stained = mean( keep.raw ),
                  kept.unstained = mean( keep.unst ),
                  row.names = NULL, stringsAsFactors = FALSE )
    } )

    g.gate <- do.call( rbind, gate.rows )
    cat( sprintf( "\n-- %s scatter gate --\n", tg ) )
    print( g.gate, digits = 3 )

    dec  <- deconvolve.af.background( p$raw,       p$s.start.panel, p$af.basis,
                                      af.name = NULL )
    decu <- deconvolve.af.background( p$unstained, p$s.start.panel, p$af.basis,
                                      af.name = NULL )

    thr <- fx.threshold.margin *
      apply( decu$abundance, 2, stats::quantile,
             probs = fx.threshold.probs, names = FALSE )
    names( thr ) <- p$panel

    thr.mat <- get.spread.thresholds(
      unmixed = dec$abundance, thresholds = thr,
      spillover.spread = p$spillover.spread,
      spread.kappa = fx.spread.kappa, verbose = FALSE )

    dom <- .fx3.dominant( dec$abundance, thr, thr.mat )

    g.thr <- do.call( rbind, lapply( seq_along( p$panel ), function( k ) {

      j   <- p$panel[ k ]
      own <- which( dom == k )

      inflation <- if ( length( own ) > 0 )
        stats::median( thr.mat[ own, j ] ) / max( abs( thr[ j ] ),
                                                  .Machine$double.eps ) else
          NA_real_

      x <- dec$abundance[ , j ]

      data.frame(
        target        = tg,
        fluorophore   = j,
        thr.flat      = unname( thr[ j ] ),
        spread.inflation = inflation,
        n.dominant    = length( own ),
        n.above.flat  = sum( x > thr[ j ] ),
        n.above.spread = sum( x > thr.mat[ , j ] ),
        x.high        = stats::quantile( x, 0.999, names = FALSE ),
        row.names = NULL, stringsAsFactors = FALSE )
    } ) )

    cat( sprintf( "\n-- %s thresholds --\n", tg ) )
    print( g.thr, digits = 3 )
    cat( sprintf( paste0( "  median spread inflation %.2f; dyes with fewer ",
                          "than 2400 events above the spread boundary: %d ",
                          "of %d\n" ),
                  stats::median( g.thr$spread.inflation, na.rm = TRUE ),
                  sum( g.thr$n.above.spread < 2400L ),
                  nrow( g.thr ) ) )

    fx.csv.g.gate[[ tg ]] <- g.gate
    fx.csv.g.thr[[ tg ]]  <- g.thr

  } )

  fx.log.time( "G", tg, t.elapsed[ "elapsed" ] )
}

fx.write.csv( do.call( rbind, fx.csv.g.gate ), "stage_G_scatter_gate.csv" )
fx.write.csv( do.call( rbind, fx.csv.g.thr ),  "stage_G_thresholds.csv" )
gc()

# ---------------------------------------------------------------------------
# STAGE C - corrected estimator against the induced truth
# ---------------------------------------------------------------------------
# The iterated fixed point with the mask, the raised cap and the guarded
# coverage gate, compared against the induced spillover. Also: a single-pass
# run at the production truncation cap (20,000), to settle whether the cap
# matters once the mask is on, and the bright-scale Huber arm from stage_AB
# on the pairs above |induced| 0.05, to settle whether the scale fix still
# buys anything on top of the mask.

cat( "\n===================== STAGE C: corrected estimator =====================\n" )

fx.stage.c <- list()
fx.csv.c   <- list()
fx.csv.c.arms <- list()

for ( tg in names( fx.prep ) ) {

  t.elapsed <- system.time( {

    p <- fx.prep[[ tg ]]

    dec  <- deconvolve.af.background( p$raw,       p$s.start.panel, p$af.basis,
                                      af.name = NULL )
    decu <- deconvolve.af.background( p$unstained, p$s.start.panel, p$af.basis,
                                      af.name = NULL )

    cat( sprintf( "\n-- %s --\n", tg ) )

    ph <- .fx.estimate.spillover( p, dec$abundance, decu$abundance )

    m.induced <- .fx.induced.spillover( p )

    log <- ph$log
    ij  <- cbind( match( log$source, p$panel ), match( log$channel, p$panel ) )
    log$estimate <- ph$spillover[ ij ]
    log$induced  <- m.induced[ ij ]
    log$err      <- log$estimate - log$induced

    bands <- cut( abs( log$induced ),
                  breaks = c( 0, 0.005, 0.02, 0.05, 0.1, 0.2, Inf ),
                  include.lowest = TRUE )

    for ( b in levels( bands ) ) {
      take <- !is.na( bands ) & bands == b
      if ( sum( take ) == 0L ) next
      mi <- stats::median( abs( log$induced[ take ] ) )
      me <- stats::median( abs( log$err[ take ] ) )
      cat( sprintf( "  band %-14s n %4d  med |induced| %.5f  med |err| %.5f  recovered %.3f\n",
                    b, sum( take ), mi, me, 1 - me / mi ) )
    }

    cat( sprintf( "  cor(estimate, induced) %.3f; trust == 0 on %d of %d pairs\n",
                  suppressWarnings( stats::cor( log$estimate, log$induced,
                                                use = "complete.obs" ) ),
                  sum( log$trust == 0 ), nrow( log ) ) )

    # Single-pass cap comparison, mask on.
    cap.rows <- lapply( c( 20000L, 400000L ), function( cap ) {
      cap.hold <- fx.max.truncated
      fx.max.truncated <<- cap
      one <- .fx.estimate.spillover( p, dec$abundance, decu$abundance,
                                     max.iter = 1L, verbose = FALSE )
      fx.max.truncated <<- cap.hold
      one.log <- one$log
      ij1 <- cbind( match( one.log$source, p$panel ),
                    match( one.log$channel, p$panel ) )
      err1 <- one.log$slope.truncated - m.induced[ ij1 ]
      big  <- abs( m.induced[ ij1 ] ) > 0.05
      data.frame( target = tg, cap = cap,
                  med.err.all = stats::median( abs( err1 ), na.rm = TRUE ),
                  med.err.big = stats::median( abs( err1[ big ] ), na.rm = TRUE ),
                  cor.big = suppressWarnings( stats::cor(
                    one.log$slope.truncated[ big ], m.induced[ ij1 ][ big ],
                    use = "complete.obs" ) ),
                  row.names = NULL, stringsAsFactors = FALSE )
    } )
    cap.table <- do.call( rbind, cap.rows )
    cat( "  single-pass truncation cap comparison, mask on:\n" )
    print( cap.table, digits = 3 )

    # Bright-scale Huber arm on the big pairs, mask-consistent: an event is
    # kept when the source is dominant for it or it is below the source
    # boundary, the same rule .fix.envelope.slope() applies.
    unmixed.c <- dec$abundance %*% solve( ph$spillover )
    colnames( unmixed.c ) <- p$panel
    unst.c <- decu$abundance %*% solve( ph$spillover )
    colnames( unst.c ) <- p$panel

    thr.c <- fx.threshold.margin *
      apply( unst.c, 2, stats::quantile, probs = fx.threshold.probs,
             names = FALSE )
    names( thr.c ) <- p$panel

    thr.mat.c <- get.spread.thresholds(
      unmixed = unmixed.c, thresholds = thr.c,
      spillover.spread = p$spillover.spread,
      spread.kappa = fx.spread.kappa, verbose = FALSE )

    dom.c <- .fx3.dominant( unmixed.c, thr.c, thr.mat.c )

    pairs.big <- .fx11.pairs( p$panel )
    ij.big    <- cbind( match( pairs.big$source, p$panel ),
                        match( pairs.big$channel, p$panel ) )
    pairs.big$induced <- m.induced[ ij.big ]
    pairs.big <- pairs.big[ abs( pairs.big$induced ) > 0.05, , drop = FALSE ]

    hub <- fx.lapply( seq_len( nrow( pairs.big ) ), function( i ) {

      src <- pairs.big$source[ i ]
      ch  <- pairs.big$channel[ i ]
      k   <- match( src, p$panel )

      keep <- dom.c == k | unmixed.c[ , src ] <= thr.mat.c[ , src ]

      x <- unmixed.c[ keep, src ]
      y <- unmixed.c[ keep, ch ]

      bright <- which( x > thr.mat.c[ keep, src ] )

      if ( length( bright ) < fx.min.events ) return( NA_real_ )

      .fx11.huber.slope( x, y, scale.index = bright )
    } )

    pairs.big$bright.hub <- unlist( hub )
    key <- match( paste( pairs.big$source, pairs.big$channel ),
                  paste( log$source, log$channel ) )
    pairs.big$estimate <- log$estimate[ key ]
    pairs.big$err.est  <- pairs.big$estimate - pairs.big$induced
    pairs.big$err.hub  <- pairs.big$bright.hub - pairs.big$induced

    cat( sprintf( paste0( "  big pairs (n %d): med |err| estimator %.4f, ",
                          "bright-scale huber %.4f\n" ),
                  nrow( pairs.big ),
                  stats::median( abs( pairs.big$err.est ), na.rm = TRUE ),
                  stats::median( abs( pairs.big$err.hub ), na.rm = TRUE ) ) )

    fx.stage.c[[ tg ]] <- list( spillover = ph$spillover, log = log,
                                m.induced = m.induced,
                                dec = list( abundance = dec$abundance ) )

    fx.csv.c[[ tg ]]      <- cbind( target = tg, log, stringsAsFactors = FALSE )
    fx.csv.c.arms[[ tg ]] <- list( caps = cap.table,
                                   hub = cbind( target = tg, pairs.big,
                                                stringsAsFactors = FALSE ) )

  } )

  fx.log.time( "C", tg, t.elapsed[ "elapsed" ] )
  gc()
}

fx.write.csv( do.call( rbind, fx.csv.c ), "stage_C_corrected_estimator.csv" )
fx.write.csv( do.call( rbind, lapply( fx.csv.c.arms, `[[`, "caps" ) ),
              "stage_C_cap_comparison.csv" )
fx.write.csv( do.call( rbind, lapply( fx.csv.c.arms, `[[`, "hub" ) ),
              "stage_C_bright_huber.csv" )
gc()

# ---------------------------------------------------------------------------
# STAGE L - null calibration with the corrected estimator, and its validity
# ---------------------------------------------------------------------------
# Three parts. The full-pool null run supplies m.null.bias and the per-row
# gate floors Stage R consumes, exactly as before but through the corrected
# estimator - if the old null artefact was largely mask-off contamination,
# the "largest coefficients invented on the null" table shrinks here and
# that is a result in itself. Then the held-out validity check: the null
# signature refitted on each half-pool, both signatures scored on the half
# they never saw. The subtraction Stage R applies is only sound if the
# null-corrected signature does NOT unmix unseen events better than the
# library.

cat( "\n===================== STAGE L: null calibration =====================\n" )

fx.null <- list(
  Cells = .fx.prepare( "Cells", "Cells" ),
  Beads = .fx.prepare( "Beads", "Beads" )
)

fx.stage.l  <- list()
fx.csv.l    <- list()
fx.csv.l.ho <- list()

for ( tg in names( fx.null ) ) {

  t.elapsed <- system.time( {

    p <- fx.null[[ tg ]]

    cat( sprintf( "\n-- %s --\n", tg ) )

    dec  <- deconvolve.af.background( p$raw,       p$s.start.panel, p$af.basis,
                                      af.name = NULL )
    decu <- deconvolve.af.background( p$unstained, p$s.start.panel, p$af.basis,
                                      af.name = NULL )

    ph.null <- .fx.estimate.spillover( p, dec$abundance, decu$abundance,
                                       verbose = FALSE )

    cat( "  largest coefficients invented on the null:\n" )
    print( head( ph.null$log[ order( -abs( ph.null$log$slope.truncated ) ),
                              c( "source", "channel", "slope.truncated",
                                 "coverage", "trust" ) ], 8 ),
           digits = 3 )
    cat( sprintf( "  pairs above the coefficient cap on the null: %d of %d\n",
                  sum( abs( ph.null$log$slope.truncated ) > fx.max.coefficient,
                       na.rm = TRUE ),
                  nrow( ph.null$log ) ) )

    full <- .fx11.null.signature( p, p$raw )

    cat( "\n  angular bias of the signature phase with nothing to correct:\n" )
    print( full$table, digits = 3 )

    # Held-out validity, both directions.
    set.seed( fx.seed )
    half <- sample.int( nrow( p$raw ), floor( nrow( p$raw ) / 2 ) )

    ho.rows <- lapply( 1:2, function( d ) {

      fit.idx  <- if ( d == 1L ) half else setdiff( seq_len( nrow( p$raw ) ),
                                                    half )
      hold.idx <- setdiff( seq_len( nrow( p$raw ) ), fit.idx )

      s.half <- .fx11.null.signature( p, p$raw[ fit.idx, , drop = FALSE ] )

      q.lib  <- .fx11.holdout.quality( p$raw[ hold.idx, , drop = FALSE ], p,
                                       p$s.start.panel )
      q.null <- .fx11.holdout.quality( p$raw[ hold.idx, , drop = FALSE ], p,
                                       s.half$s.null )

      data.frame(
        target          = tg,
        direction       = d,
        resid.library   = q.lib$resid,
        resid.null      = q.null$resid,
        leak.library    = sum( q.lib$leak,  na.rm = TRUE ),
        leak.null       = sum( q.null$leak, na.rm = TRUE ),
        row.names = NULL, stringsAsFactors = FALSE )
    } )

    ho <- do.call( rbind, ho.rows )

    cat( "\n  held-out validity (mean of both split directions):\n" )
    cat( sprintf( paste0( "    reconstruction residual: library %.5g, ",
                          "null-corrected %.5g\n" ),
                  mean( ho$resid.library ), mean( ho$resid.null ) ) )
    cat( sprintf( "    off-target leak:         library %.5g, null-corrected %.5g\n",
                  mean( ho$leak.library ), mean( ho$leak.null ) ) )
    cat( if ( mean( ho$leak.null ) < mean( ho$leak.library ) &&
              mean( ho$resid.null ) < mean( ho$resid.library ) )
      paste0( "    -> the null-corrected signature fits unseen events better ",
              "than the library.\n",
              "       The null 'bias' contains real library error; treat Stage ",
              "R's subtraction\n       and the bias gate with caution on this ",
              "substrate.\n" ) else
      paste0( "    -> no improvement on unseen events; consistent with pure ",
              "estimator\n       artefact, and the subtraction stands.\n" ) )

    fx.stage.l[[ tg ]] <- list( signature = full$table, s.null = full$s.null,
                                m.spillover = ph.null$spillover )

    fx.csv.l[[ tg ]]    <- cbind( target = tg, full$table,
                                  stringsAsFactors = FALSE )
    fx.csv.l.ho[[ tg ]] <- ho

  } )

  fx.log.time( "L", tg, t.elapsed[ "elapsed" ] )
  gc()
}

fx.write.csv( do.call( rbind, fx.csv.l ),    "stage_L_null_signature_bias.csv" )
fx.write.csv( do.call( rbind, fx.csv.l.ho ), "stage_L_holdout_validity.csv" )
gc()

# ---------------------------------------------------------------------------
# STAGE R - the corrected alternating loop
# ---------------------------------------------------------------------------
# A8's Stage R with the four corrections applied: the masked estimator, the
# raised cap, per-row null-referenced gate cuts, and retirement that freezes
# a converged row against the spillover fold-back as well as the signature
# phase.

cat( "\n===================== STAGE R: alternating passes =====================\n" )

fx.stage.r  <- list()
fx.csv.r.gate <- list()
fx.csv.r.pass <- list()
fx.csv.r.track <- list()

for ( tg in names( fx.prep ) ) {

  p <- fx.prep[[ tg ]]

  shared <- intersect( rownames( fx.stage.l[[ tg ]]$s.null ), p$panel )

  bias <- fx.stage.l[[ tg ]]$s.null[ shared, p$detectors, drop = FALSE ] -
    p$s.start.panel[ shared, , drop = FALSE ]

  m.null.bias <- NULL
  if ( !is.null( fx.stage.l[[ tg ]]$m.spillover ) ) {
    m.null.bias <- fx.stage.l[[ tg ]]$m.spillover -
      diag( nrow( fx.stage.l[[ tg ]]$m.spillover ) )
    dimnames( m.null.bias ) <- dimnames( fx.stage.l[[ tg ]]$m.spillover )
  }

  reference     <- fx.stage.l[[ tg ]]$signature
  cut.resid     <- 3 * stats::median( reference$resid.rel,     na.rm = TRUE )
  cut.intercept <- 3 * stats::median( reference$intercept.rel, na.rm = TRUE )

  # Per-row cuts: each row scored against its own null floor, never tighter
  # than the panel rule.
  null.row <- function( column, cut.panel ) {
    v <- stats::setNames( rep( NA_real_, length( p$panel ) ), p$panel )
    k <- match( p$panel, reference$fluorophore )
    v[ !is.na( k ) ] <- reference[[ column ]][ k[ !is.na( k ) ] ]
    v[ !is.finite( v ) ] <- 0
    pmax( fx.row.ratio * v, cut.panel )
  }

  cut.resid.row     <- null.row( "resid.rel",     cut.resid )
  cut.intercept.row <- null.row( "intercept.rel", cut.intercept )

  null.explained <- stats::setNames( rep( 1, length( p$panel ) ), p$panel )
  k.exp <- match( p$panel, reference$fluorophore )
  null.explained[ !is.na( k.exp ) ] <-
    reference$explained.total[ k.exp[ !is.na( k.exp ) ] ]
  null.explained[ !is.finite( null.explained ) ] <- 1

  design.hot <- rbind( p$af.basis, p$s.start.panel )
  hot        <- calculate.hotspot.matrix( design.hot )
  coupling   <- apply( hot[ p$panel, rownames( p$af.basis ), drop = FALSE ],
                       1, max )
  af.frozen  <- names( coupling )[ coupling > fx.max.hotspot ]

  if ( length( af.frozen ) > 0 )
    cat( sprintf( "  %s: %d fluorophore(s) frozen for AF coupling above %.1f: %s\n",
                  tg, length( af.frozen ), fx.max.hotspot,
                  paste( af.frozen, collapse = ", " ) ) )

  s.curr <- p$s.start.panel
  track  <- data.frame( fluorophore = p$panel,
                        pass0 = .fx.angle( s.curr, p$s.true.panel ),
                        row.names = NULL )
  
  ph <- NULL
  
  # Initialised once per target here, not inside the pass loop gated on
  # pass == 1L: a pass that accepts nothing breaks the loop before that gate
  # is reached, which previously left these completely unset for this
  # target and let a same-session leftover from an unrelated run stand in
  # silently instead of failing loudly.
  s.by.pass     <- list()
  delta.by.pass <- numeric( 0 )
  step.by.pass  <- numeric( 0 )
  n.acc.by.pass <- integer( 0 )
  fx.csv.r.gate.tg <- list()

  for ( pass in seq_len( fx.n.pass ) ) {

    t.elapsed.pass <- system.time( {

    dec  <- deconvolve.af.background( p$raw,       s.curr, p$af.basis,
                                      af.name = NULL )
    decu <- deconvolve.af.background( p$unstained, s.curr, p$af.basis,
                                      af.name = NULL )

    set.seed( fx.seed + pass )
    fit.idx <- if ( nrow( dec$abundance ) > fx.pass.events )
      sample.int( nrow( dec$abundance ), fx.pass.events ) else
        seq_len( nrow( dec$abundance ) )

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

      if ( j %in% af.frozen ) next

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
      step.impact <- unname( x.high[ j ] ) *
        sqrt( sum( ( ( cand$signature - s.curr[ j, ] ) %*%
                       t( unmixing.curr ) )^2 ) )

      bias.impact <- if ( j %in% shared )
        unname( x.high[ j ] ) *
        sqrt( sum( ( bias[ j, ] %*% t( unmixing.curr ) )^2 ) ) else 0

      drift <- .fx.angle( s.curr[ j, , drop = FALSE ],
                          p$s.start.panel[ j, , drop = FALSE ] )

      reason <- if ( !isTRUE( st$joint ) ) "univariate" else "accepted"

      if ( j %in% settled ) {
        reason <- "settled"
      } else if ( st$x.span <= fx.min.span * abs( thr[ j ] ) ) {
        reason <- "span"
      } else if ( !is.finite( st$explained.total ) ||
                  abs( st$explained.total - null.explained[[ j ]] ) >
                  fx.max.explained.drift ) {
        reason <- "explained"
      } else if ( st$resid.rel > cut.resid.row[[ j ]] ) {
        reason <- "fit"
      } else if ( st$intercept.rel > cut.intercept.row[[ j ]] ) {
        reason <- "offset"
      } else if ( is.finite( st$bg.align ) && st$bg.align < fx.min.bg.align ) {
        reason <- "background.confound"
      } else if ( st$clamp.frac > 0.15 ) {
        reason <- "negative.mass"
      } else if ( step.impact < fx.impact.ratio * bias.impact ) {
        reason <- "bias"
      } else if ( st$deg.change > fx.max.angle ) {
        reason <- "angle"
      } else if ( is.finite( st$anchor.rel ) && st$anchor.rel > fx.max.anchor ) {
        reason <- "background"
      } else if ( st$vif.target > fx.max.vif ) {
        reason <- "collinear"
      } else if ( st$peak.new != st$peak.curr &&
                  st$peak.new.rel < fx.peak.shift.min.rel ) {
        reason <- "peak.shift"
      } else if ( step.impact >= fx.step.decay * step.prev[ j ] ) {
        reason <- "plateau"
      } else if ( drift > fx.total.max.angle ) {
        reason <- "drift"
      }

      gate.rows[[ length( gate.rows ) + 1L ]] <- data.frame(
        fluorophore = j, step.impact = step.impact,
        bias.impact = bias.impact,
        ratio       = step.impact / pmax( bias.impact, .Machine$double.eps ),
        step.prev   = unname( step.prev[ j ] ),
        joint       = isTRUE( st$joint ),
        reason      = reason,
        row.names   = NULL, stringsAsFactors = FALSE )

      # A row whose step has stopped shrinking has taken everything this
      # control can tell it, and is retired from the signature phase and the
      # spillover fold-back together: a retired row left out of `frozen`
      # would keep having a freshly estimated spillover row folded into its
      # signature every pass, replacing convergence with a random walk.
      if ( reason %in% c( "plateau", "settled" ) )
        settled <- union( settled, j )

      if ( ! reason %in% c( "accepted", "univariate" ) ) next

      step.prev[ j ] <- step.impact
      step.total     <- step.total + step.impact

      s.raw[ j, ] <- cand$signature
      accepted    <- c( accepted, j )
    }

    ever.accepted <- union( ever.accepted, accepted )

    cat( sprintf( "  %s pass %d: %d row(s) accepted, spillover delta %.5f\n",
                  tg, pass, length( accepted ), ph$delta ) )

    unstable.this.pass <- character( 0 )

    if ( length( gate.rows ) > 0 ) {
      gate.table <- do.call( rbind, gate.rows )
      print( gate.table, digits = 3 )
      fx.csv.r.gate.tg[[ pass ]] <- cbind(
        target = tg, pass = pass, gate.table, stringsAsFactors = FALSE )

      unstable.this.pass <- gate.table$fluorophore[
        ! gate.table$reason %in% c( "accepted", "univariate",
                                    "plateau", "settled" ) ]
    }

    if ( length( accepted ) == 0L ) {
      cat( sprintf( "  %s pass %d: nothing accepted, stopping here\n", tg, pass ) )
      break
    }

    spillover.pass <- ph$spillover

    if ( !is.null( m.null.bias ) ) {
      shared.m <- intersect( rownames( m.null.bias ), p$panel )
      spillover.pass[ shared.m, shared.m ] <-
        spillover.pass[ shared.m, shared.m ] -
        m.null.bias[ shared.m, shared.m, drop = FALSE ]
      diag( spillover.pass ) <- 1
    }

    frozen <- union( union( union( af.frozen,
                                   setdiff( p$panel, ever.accepted ) ),
                            unstable.this.pass ),
                     settled )

    if ( length( frozen ) > 0 ) {
      spillover.pass[ frozen, ] <- 0
      spillover.pass[ cbind( frozen, frozen ) ] <- 1
    }

    back.raw <- spillover.pass %*% s.curr

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

    s.by.pass[[ pass ]]   <- s.curr
    delta.by.pass[ pass ] <- ph$delta
    step.by.pass[ pass ]  <- step.total
    n.acc.by.pass[ pass ] <- length( accepted )

    } )

    fx.log.time( "R", sprintf( "%s pass %d", tg, pass ),
                 t.elapsed.pass[ "elapsed" ] )
    gc()
  }

  if ( length( s.by.pass ) == 0L ) {
    cat( sprintf( paste0( "  %s: nothing accepted on pass 1; no signature-",
                          "phase correction ran. Reporting the starting ",
                          "reference unchanged.\n" ), tg ) )
    s.by.pass[[ 1L ]]   <- s.curr
    delta.by.pass[ 1L ] <- NA_real_
    step.by.pass[ 1L ]  <- 0
    n.acc.by.pass[ 1L ] <- 0L
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

  impact.row <- function( s ) x.max * sqrt( rowSums(
    ( ( s[ p$panel, , drop = FALSE ] - p$s.true.panel ) %*%
        t( unmixing.q ) )^2 ) )

  impact.start <- impact.row( p$s.start.panel )

  pass.table <- data.frame(
    pass       = seq_along( s.by.pass ),
    n.accepted = n.acc.by.pass,
    delta      = delta.by.pass,
    step       = step.by.pass,
    angle      = vapply( s.by.pass, function( s )
      sum( .fx.angle( s, p$s.true.panel ) ), numeric( 1 ) ),
    impact     = vapply( s.by.pass, function( s ) sum( impact.row( s ) ),
                         numeric( 1 ) ),
    row.names = NULL )

  pass.table$recovered <- 1 - pass.table$impact / sum( impact.start )

  cat( sprintf( "  start: angle %.2f deg, impact %.4g\n",
                sum( track$pass0 ), sum( impact.start ) ) )
  print( pass.table, digits = 3 )

  best.pass    <- which.min( pass.table$impact )
  impact.final <- impact.row( s.by.pass[[ best.pass ]] )

  u.table <- data.frame(
    fluorophore   = p$panel,
    deg.start     = track$pass0,
    deg.final     = .fx.angle( s.by.pass[[ best.pass ]], p$s.true.panel ),
    impact.start  = unname( impact.start ),
    impact.final  = unname( impact.final ),
    recovered     = 1 - unname( impact.final ) / unname( impact.start ),
    ever.accepted = p$panel %in% ever.accepted,
    af.frozen     = p$panel %in% af.frozen,
    row.names = NULL )

  cat( "\n  per-dye recovery at the best pass:\n" )
  print( u.table[ order( -u.table$impact.final ), ], digits = 3 )

  fx.stage.r[[ tg ]] <- list( s.final = s.by.pass[[ best.pass ]],
                              u.table = u.table )

  fx.csv.r.gate[[ tg ]]  <- do.call( rbind, fx.csv.r.gate.tg )
  fx.csv.r.pass[[ tg ]]  <- cbind( target = tg, pass.table,
                                   stringsAsFactors = FALSE )
  fx.csv.r.track[[ tg ]] <- cbind( target = tg, u.table,
                                   stringsAsFactors = FALSE )
}

fx.write.csv( do.call( rbind, fx.csv.r.gate ),  "stage_R_gate_detail.csv" )
fx.write.csv( do.call( rbind, fx.csv.r.pass ),  "stage_R_pass_table.csv" )
fx.write.csv( do.call( rbind, fx.csv.r.track ), "stage_U_remaining_error.csv" )
gc()

# ---------------------------------------------------------------------------
# STAGE G2 - threshold audit against the corrected reference
# ---------------------------------------------------------------------------
# Stage G's inflation is measured against p$s.start.panel, the reference
# Stage R exists to correct. .align.spillover.spread() zeroes the diagonal,
# so a dye's own on-target abundance never inflates its own threshold - only
# cross-talk from every other channel does, read off whatever reference is
# currently in use. If Beads' inflation is mostly the starting reference's
# own unresolved spillover showing up as apparent cross-talk, it should
# shrink once a better reference is in play. This re-runs Stage G's
# threshold half against Stage R's corrected signature
# (fx.stage.r[[tg]]$s.final) to test that directly.

cat( "\n===================== STAGE G2: thresholds after correction =====================\n" )

fx.csv.g2.thr <- list()

for ( tg in names( fx.prep ) ) {
  
  t.elapsed <- system.time( {
    
    p      <- fx.prep[[ tg ]]
    s.corr <- fx.stage.r[[ tg ]]$s.final
    
    dec  <- deconvolve.af.background( p$raw,       s.corr, p$af.basis,
                                      af.name = NULL )
    decu <- deconvolve.af.background( p$unstained, s.corr, p$af.basis,
                                      af.name = NULL )
    
    thr <- fx.threshold.margin *
      apply( decu$abundance, 2, stats::quantile,
             probs = fx.threshold.probs, names = FALSE )
    names( thr ) <- p$panel
    
    thr.mat <- get.spread.thresholds(
      unmixed = dec$abundance, thresholds = thr,
      spillover.spread = p$spillover.spread,
      spread.kappa = fx.spread.kappa, verbose = FALSE )
    
    dom <- .fx3.dominant( dec$abundance, thr, thr.mat )
    
    g2.thr <- do.call( rbind, lapply( seq_along( p$panel ), function( k ) {
      
      j   <- p$panel[ k ]
      own <- which( dom == k )
      
      inflation <- if ( length( own ) > 0 )
        stats::median( thr.mat[ own, j ] ) / max( abs( thr[ j ] ),
                                                  .Machine$double.eps ) else
                                                    NA_real_
      
      x <- dec$abundance[ , j ]
      
      data.frame(
        target           = tg,
        fluorophore      = j,
        thr.flat         = unname( thr[ j ] ),
        spread.inflation = inflation,
        n.dominant       = length( own ),
        n.above.flat     = sum( x > thr[ j ] ),
        n.above.spread   = sum( x > thr.mat[ , j ] ),
        row.names = NULL, stringsAsFactors = FALSE )
    } ) )
    
    cat( sprintf( "\n-- %s thresholds, corrected reference --\n", tg ) )
    print( g2.thr, digits = 3 )
    cat( sprintf( "  median spread inflation %.2f (compare stage_G_thresholds.csv)\n",
                  stats::median( g2.thr$spread.inflation, na.rm = TRUE ) ) )
    
    fx.csv.g2.thr[[ tg ]] <- g2.thr
    
  } )
  
  fx.log.time( "G2", tg, t.elapsed[ "elapsed" ] )
}

fx.write.csv( do.call( rbind, fx.csv.g2.thr ), "stage_G2_thresholds_corrected.csv" )
gc()
cat( "\n  READ FIRST: compare spread.inflation and n.above.spread here against\n" )
cat( "  stage_G_thresholds.csv row for row. If Beads' inflation drops toward\n" )
cat( "  Cells' range once the corrected reference is in use, the threshold\n" )
cat( "  collapse is mostly a starting-reference artefact and the fix belongs\n" )
cat( "  in the pass loop, not in spread.kappa. If it stays high, treat it as\n" )
cat( "  a real substrate property and spread.kappa needs its own fix.\n" )

# ---------------------------------------------------------------------------
# STAGE Z - the within-cluster shear as the rescue path (idea 5)
# ---------------------------------------------------------------------------
# Section 10 item 6 of CONTEXT_information_sources.md, re-run against the
# corrected estimator's own refusals. The shear earned this re-test by a 72%
# error reduction on refused-and-real pairs at n = 1000 in stage_AE; whether
# it still earns a place depends on what the corrected estimator leaves
# refused. Scored per cluster setting on: all pairs, fitted pairs, and
# refused-and-real pairs, plus the combined rule (estimator where trusted,
# shear where refused) that would actually ship.

cat( "\n===================== STAGE Z: within-cluster shear =====================\n" )

fx.csv.z <- list()

for ( tg in names( fx.prep ) ) {

  t.elapsed <- system.time( {

    p         <- fx.prep[[ tg ]]
    m.induced <- fx.stage.c[[ tg ]]$m.induced
    est.log   <- fx.stage.c[[ tg ]]$log

    u.full <- fx.stage.c[[ tg ]]$dec$abundance %*%
      solve( fx.stage.c[[ tg ]]$spillover )
    colnames( u.full ) <- p$panel

    set.seed( fx.seed )
    sub <- if ( nrow( u.full ) > fx.stat.events )
      sort( sample.int( nrow( u.full ), fx.stat.events ) ) else
        seq_len( nrow( u.full ) )

    u <- u.full[ sub, , drop = FALSE ]

    decu <- deconvolve.af.background( p$unstained, p$s.start.panel, p$af.basis,
                                      af.name = NULL )
    uu <- decu$abundance %*% solve( fx.stage.c[[ tg ]]$spillover )
    colnames( uu ) <- p$panel

    thr <- fx.threshold.margin *
      apply( uu, 2, stats::quantile, probs = fx.threshold.probs, names = FALSE )
    names( thr ) <- p$panel

    thr.mat <- get.spread.thresholds(
      unmixed = u, thresholds = thr,
      spillover.spread = p$spillover.spread,
      spread.kappa = fx.spread.kappa, verbose = FALSE )

    pairs <- .fx11.pairs( p$panel )
    ij    <- cbind( match( pairs$source, p$panel ),
                    match( pairs$channel, p$panel ) )
    pairs$induced  <- m.induced[ ij ]
    pairs$estimate <- fx.stage.c[[ tg ]]$spillover[ ij ]

    key <- match( paste( pairs$source, pairs$channel ),
                  paste( est.log$source, est.log$channel ) )
    pairs$fitted <- !is.na( key ) & est.log$trust[ key ] > 0
    pairs$fitted[ is.na( pairs$fitted ) ] <- FALSE

    med.ind <- stats::median( abs( pairs$induced ), na.rm = TRUE )
    refused <- !pairs$fitted & abs( pairs$induced ) > med.ind

    cat( sprintf( "\n-- %s --  %d pairs refused by the corrected estimator, %d of them real\n",
                  tg, sum( !pairs$fitted ), sum( refused ) ) )

    for ( nc in fx.cluster.grid ) {

      id <- .fx11.cluster( u, n.clusters = nc, min.size = fx.min.cluster,
                           seed = fx.seed )

      out <- fx.lapply( seq_len( nrow( pairs ) ), function( i ) {
        .fx11.cluster.shear( u[ , pairs$source[ i ] ],
                             u[ , pairs$channel[ i ] ],
                             thr.mat[ , pairs$source[ i ] ],
                             id, min.cluster = fx.min.cluster )
      } )
      out <- do.call( rbind, out )

      err.mask   <- pairs$estimate - pairs$induced
      err.shear  <- out[ , "m.within" ] - pairs$induced
      combined   <- ifelse( pairs$fitted, pairs$estimate,
                            out[ , "m.within" ] )
      combined[ !is.finite( combined ) ] <- pairs$estimate[
        !is.finite( combined ) ]
      err.comb   <- combined - pairs$induced

      report <- function( take, label ) {
        if ( sum( take, na.rm = TRUE ) == 0L ) return( invisible( NULL ) )
        cat( sprintf( paste0( "  n %4d clusters, %-16s n %4d | med |err| ",
                              "mask %.5f  shear %.5f  combined %.5f\n" ),
                      nc, label, sum( take, na.rm = TRUE ),
                      stats::median( abs( err.mask[ take ] ),  na.rm = TRUE ),
                      stats::median( abs( err.shear[ take ] ), na.rm = TRUE ),
                      stats::median( abs( err.comb[ take ] ),  na.rm = TRUE ) ) )
      }

      report( rep( TRUE, nrow( pairs ) ), "all pairs" )
      report( pairs$fitted, "fitted" )
      report( refused, "refused, real" )

      fx.csv.z[[ length( fx.csv.z ) + 1L ]] <- data.frame(
        target = tg, n.clusters = nc,
        pairs, m.within = out[ , "m.within" ], m.floor = out[ , "m.floor" ],
        shear.clusters = out[ , "n.clusters" ],
        err.mask = err.mask, err.shear = err.shear, err.combined = err.comb,
        row.names = NULL, stringsAsFactors = FALSE )
    }

  } )

  fx.log.time( "Z", tg, t.elapsed[ "elapsed" ] )
  gc()
}

fx.write.csv( do.call( rbind, fx.csv.z ), "stage_Z_shear_rescue.csv" )
gc()

# ---------------------------------------------------------------------------
# STAGE W - the tail-excess detector on the corrected refusals (idea 1)
# ---------------------------------------------------------------------------
# Section 10 item 7, gate-rescue half only: the sign-dependent trust term was
# refuted in A9/A10 and is not re-tested. The symmetric statistic
# max( |hyper|, |hyper.high| ) detected damaging-unfitted pairs at AUC 0.75 /
# 0.87 against the mask-off estimator; re-scored here against the corrected
# estimator's refusals, which is the set a gate rescue would actually act on.

cat( "\n===================== STAGE W: tail-excess detector =====================\n" )

fx.csv.w <- list()

for ( tg in names( fx.prep ) ) {

  t.elapsed <- system.time( {

    p         <- fx.prep[[ tg ]]
    m.induced <- fx.stage.c[[ tg ]]$m.induced
    est.log   <- fx.stage.c[[ tg ]]$log

    u.full <- fx.stage.c[[ tg ]]$dec$abundance %*%
      solve( fx.stage.c[[ tg ]]$spillover )
    colnames( u.full ) <- p$panel

    set.seed( fx.seed )
    sub <- if ( nrow( u.full ) > fx.stat.events )
      sort( sample.int( nrow( u.full ), fx.stat.events ) ) else
        seq_len( nrow( u.full ) )

    u <- u.full[ sub, , drop = FALSE ]

    decu <- deconvolve.af.background( p$unstained, p$s.start.panel, p$af.basis,
                                      af.name = NULL )
    uu <- decu$abundance %*% solve( fx.stage.c[[ tg ]]$spillover )
    colnames( uu ) <- p$panel

    thr <- fx.threshold.margin *
      apply( uu, 2, stats::quantile, probs = fx.threshold.probs, names = FALSE )
    names( thr ) <- p$panel

    thr.mat <- get.spread.thresholds(
      unmixed = u, thresholds = thr,
      spillover.spread = p$spillover.spread,
      spread.kappa = fx.spread.kappa, verbose = FALSE )

    floor.low  <- apply( uu, 2, stats::quantile, probs = fx.hyper.q,
                         names = FALSE )
    floor.high <- apply( uu, 2, stats::quantile, probs = fx.upper.q,
                         names = FALSE )
    names( floor.low ) <- names( floor.high ) <- p$panel

    pairs <- .fx11.pairs( p$panel )

    out <- fx.lapply( seq_len( nrow( pairs ) ), function( i ) {

      src <- pairs$source[ i ]
      ch  <- pairs$channel[ i ]

      x <- u[ , src ]
      y <- u[ , ch ]

      positive <- x > thr.mat[ , src ]
      negative <- !positive

      if ( sum( positive ) < fx.min.events || sum( negative ) < fx.min.events )
        return( c( hyper = NA_real_, hyper.high = NA_real_ ) )

      cut.bright <- stats::quantile( x[ positive ],
                                     probs = 1 - fx.bright.frac,
                                     names = FALSE )
      bright <- positive & x >= cut.bright

      if ( sum( bright ) < fx.min.events )
        return( c( hyper = NA_real_, hyper.high = NA_real_ ) )

      c( hyper = mean( y[ bright ] < floor.low[ ch ] ) -
           mean( y[ negative ] < floor.low[ ch ] ),
         hyper.high = mean( y[ bright ] > floor.high[ ch ] ) -
           mean( y[ negative ] > floor.high[ ch ] ) )
    } )

    w.table <- cbind( pairs, do.call( rbind, out ) )

    ij <- cbind( match( w.table$source, p$panel ),
                 match( w.table$channel, p$panel ) )
    w.table$induced  <- m.induced[ ij ]
    w.table$estimate <- fx.stage.c[[ tg ]]$spillover[ ij ]

    key <- match( paste( w.table$source, w.table$channel ),
                  paste( est.log$source, est.log$channel ) )
    w.table$fitted <- !is.na( key ) & est.log$trust[ key ] > 0
    w.table$fitted[ is.na( w.table$fitted ) ] <- FALSE

    w.table$tail.excess <- pmax( abs( w.table$hyper ),
                                 abs( w.table$hyper.high ) )

    med.ind  <- stats::median( abs( w.table$induced ), na.rm = TRUE )
    damaging <- !w.table$fitted & abs( w.table$induced ) > med.ind

    ok <- is.finite( w.table$tail.excess )

    auc <- if ( sum( damaging & ok ) > 0 && sum( !damaging & ok ) > 0 )
      mean( outer( w.table$tail.excess[ damaging & ok ],
                   w.table$tail.excess[ !damaging & ok ], ">" ) ) else NA_real_

    cat( sprintf( paste0( "\n-- %s --  %d damaging-and-refused pairs; ",
                          "tail excess AUC %.3f\n" ),
                  tg, sum( damaging, na.rm = TRUE ), auc ) )
    cat( sprintf( "  median tail excess on them %.4f, elsewhere %.4f\n",
                  stats::median( w.table$tail.excess[ damaging ], na.rm = TRUE ),
                  stats::median( w.table$tail.excess[ !damaging ], na.rm = TRUE ) ) )

    fx.csv.w[[ tg ]] <- cbind( target = tg, w.table, damaging = damaging,
                               stringsAsFactors = FALSE )

  } )

  fx.log.time( "W", tg, t.elapsed[ "elapsed" ] )
  gc()
}

fx.write.csv( do.call( rbind, fx.csv.w ), "stage_W_tail_excess.csv" )
gc()

# ---------------------------------------------------------------------------
# Timing
# ---------------------------------------------------------------------------

cat( "\n===================== A11 timing =====================\n" )
print( fx.timing, digits = 3 )
cat( sprintf( "  total %.1f s\n", sum( fx.timing$elapsed ) ) )
fx.write.csv( fx.timing, "stage_timing_A11.csv" )

cat( "\nReport back these files from", fx.output.dir, ":\n" )
cat( "  stage_G_scatter_gate.csv\n" )
cat( "  stage_G_thresholds.csv\n" )
cat( "  stage_C_corrected_estimator.csv\n" )
cat( "  stage_C_cap_comparison.csv\n" )
cat( "  stage_C_bright_huber.csv\n" )
cat( "  stage_L_null_signature_bias.csv\n" )
cat( "  stage_L_holdout_validity.csv\n" )
cat( "  stage_R_gate_detail.csv\n" )
cat( "  stage_R_pass_table.csv\n" )
cat( "  stage_U_remaining_error.csv\n" )
cat( "  stage_Z_shear_rescue.csv\n" )
cat( "  stage_W_tail_excess.csv\n" )
cat( "  stage_timing_A11.csv\n" )

cat( "\n===================== STAGE F: spectral traces =====================\n" )
trace.dir <- "./figure_fix_my_unmix_A11"
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
