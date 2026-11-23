# dev/prototype_af_shrinkage.R
#
# Prototype of a noise-aware per-cell AF abundance, tested against
# sim.flow.data() ground truth on a real panel's spectra and AF library.
#
# In plain terms. A cell's AF abundance is read from the part of its spectrum
# that no dye in the panel can explain. On a cell that is bright for a dye,
# that part also carries the dye's photon noise, so the AF estimate becomes
# very noisy (on 7C, the AF spread on the brightest SBV515 cells is about 8
# times that of unstained cells, and a third of them are clamped at 0).
# Clamping a noisy estimate at 0 does not remove the noise; it piles cells up
# at zero and leaves their AF in the dye channels. The idea tested here is the
# one a careful analyst would use by eye: when a cell's own spectrum says
# little about its AF, assume it has about as much AF as similar cells, and
# trust its own spectrum more the more informative it is.
#
# The maths. For a cell with chosen AF variant j, the unclamped AF estimate is
# the Frisch-Waugh projection k_hat = r_j' y / D_j (r_j: the part of the AF
# spectrum outside the panel's span; D_j = r_j' r_j, floored as in the app).
# Its sampling variance is
#     s2 = sum_d r_jd^2 var_d / D_j^2,
# with var_d the detector noise at the cell's fitted signal (readout variance
# plus photon and spillover-sampling variance; here the simulation's own
# noise model, in Honeychrome the fitted noise model). The prior is the AF of
# similar cells: k ~ Normal(m, tau2) truncated at 0, with m a linear function
# of log SSC and tau2 the population AF variance with the estimation noise
# removed (both fitted on the sample's own cells; only cells with below-median
# s2 are used, as they carry the AF information). The posterior is then a
# normal truncated at 0 with
#     w      = tau2 / (tau2 + s2)            (weight on the cell's own estimate)
#     mu_p   = m + w (k_hat - m)
#     sd_p   = sqrt(w s2)
# and the reported AF is its mean,
#     k_post = mu_p + sd_p * phi(a) / (1 - Phi(a)),   a = -mu_p / sd_p,
# which is positive, smooth through zero, and close to k_hat on cells whose
# own spectrum determines their AF (s2 << tau2). The fluorophores follow as
# u - k v_j (u: the plain unmix without AF; v_j: the variant's leakage into
# the fluorophores), exactly as for any other AF value.
#
# Arms:
#   A  clamp:     max(k_hat, 0)                  (the app's per-cell AF today)
#   B  unclamped: k_hat
#   S1 shrink:    mu_p (normal prior, not truncated)
#   S2 shrink:    k_post (prior truncated at 0)
#
# Measures, against the simulation's truth:
#   * AF bias and RMSE by decile of the cell's true total dye abundance, and
#     the share of cells at exactly 0
#   * calibration: robust SD of (k_hat - k_true) / sqrt(s2) by decile (about
#     1 when the variance model is right; variant mismatch inflates it)
#   * background in the fluorophores: for each fluorophore, the mean (bias)
#     and SD of its estimate on cells where it is truly off, summarised over
#     the fluorophores AF leaks into most (|v| > 0.2) and over all
#   * the same background on the cells arm A clamps (the cells seen on plots)
#
# Deliberate v1 simplifications:
#   * the variance model is the simulation's own noise (an oracle stand-in for
#     the fitted noise model), with the fitted signal clamped at 0 as M2pg's
#     weights are
#   * the prior mean depends on log SSC only; the prior is fitted on the
#     sample being unmixed
#   * AF is assigned by assign.af.joint.cov.l2() (the app's per-cell scoring
#     without the negative envelope) and the variant is kept fixed
#   * sim.flow.data() subtracts its dark offset from every detector after
#     adding zero-mean readout noise, which leaves every detector biased by
#     -dark.offset; that offset is added back here so the AF estimate is not
#     biased by a constant baseline error
#
# Data: panel_spectra.csv and af_spectra.csv written by the Honeychrome
# scenario hc_harness/scenarios/af_m2pg_spread.py; the newest results folder
# is used, and the tables are written to its af_shrinkage subfolder. Check
# section 0, then run top to bottom.

# run from the AutoSpectral sandbox project (its folder is the working
# directory); the Honeychrome harness sits beside it
asp.dir <- "."
devtools::load_all( asp.dir )

# ---------------------------------------------------------------------------
# 0. Setup -- EDIT THIS SECTION
# ---------------------------------------------------------------------------

results.root <- file.path( "..", "hc_harness", "results", "af_m2pg_spread" )
spectra.dir  <- utils::tail( sort( list.dirs( path.expand( results.root ), recursive = FALSE ) ), 1L )
out.dir      <- file.path( spectra.dir, "af_shrinkage" )
cytometer   <- "aurora"

n.cells        <- 40000L
complexity     <- 0.15                 # share of the panel on per cell
af.scale.range <- c( 1e-4, 1.1e-3 )    # x expr.data.max: AF about 400-4,600
af.scatter     <- 0.3
seed           <- 11L

# ---------------------------------------------------------------------------
# 1. Panel, AF library, simulation
# ---------------------------------------------------------------------------

read.matrix <- function( file ) {
  m <- as.matrix( utils::read.csv( file, row.names = 1, check.names = FALSE ) )
  storage.mode( m ) <- "double"
  m
}
spectra    <- read.matrix( file.path( spectra.dir, "panel_spectra.csv" ) )
af.spectra <- read.matrix( file.path( spectra.dir, "af_spectra.csv" ) )
stopifnot( identical( colnames( spectra ), colnames( af.spectra ) ) )
dir.create( out.dir, showWarnings = FALSE, recursive = TRUE )

asp    <- get.autospectral.param( cytometer = cytometer, figures = FALSE )
preset <- .sim.detector.preset( asp$cytometer )

sim <- sim.flow.data(
  spectra, asp,
  n.cells          = n.cells,
  complexity       = complexity,
  af.spectra       = af.spectra,
  af.scale.range   = af.scale.range,
  af.scatter.slope = af.scatter,
  seed             = seed
)
raw    <- sim$raw + preset$dark.offset
truth  <- sim$truth$abundances[ , rownames( spectra ), drop = FALSE ]
k.true <- sim$truth$af.abundance
ssc    <- sim$truth$scatter[ , 2L ]

# ---------------------------------------------------------------------------
# 2. Per-cell AF: assignment, unclamped abundance, its variance
# ---------------------------------------------------------------------------

af.idx <- assign.af.joint.cov.l2( raw, spectra, af.spectra )
fit    <- unmix.af.fwl( raw, spectra, af.spectra, af.idx, denominator.floor = 0.01 )
k.hat  <- fit$af

unmixing.matrix <- solve.default( tcrossprod( spectra ), spectra )
v.library   <- unmixing.matrix %*% t( af.spectra )                  # fluors x variants
r.library   <- t( af.spectra ) - t( spectra ) %*% v.library         # detectors x variants
denominator <- colSums( r.library^2 )
denominator <- pmax( denominator, 0.01 * max( denominator ) )

# detector variance at each cell's fitted signal: readout, photon and
# spillover-sampling terms of the simulation's noise. The signal is taken
# from the plain unmix without AF and the sample's median AF, so a cell's
# variance does not depend on its own AF estimate (otherwise grouping cells
# by variance would also group them by their AF estimation error).
v.cell <- t( v.library )[ af.idx, , drop = FALSE ]
x.fit  <- pmax( fit$fluorophores + k.hat * v.cell, 0 )
mu.fit <- x.fit %*% spectra + max( stats::median( k.hat ), 0 ) * af.spectra[ af.idx, , drop = FALSE ]
p.spec <- pmax( spectra, 0 ) / rowSums( pmax( spectra, 0 ) )
cpu    <- preset$counts.per.unit
det.var <- preset$readout.sd^2 + pmax( mu.fit, 0 ) / cpu +
  ( x.fit %*% ( p.spec * ( 1 - p.spec ) ) ) / cpu
s2 <- rowSums( det.var * t( r.library^2 )[ af.idx, , drop = FALSE ] ) / denominator[ af.idx ]^2

# ---------------------------------------------------------------------------
# 3. Prior from the sample's own cells, posterior
# ---------------------------------------------------------------------------

informative <- s2 <= stats::median( s2 )
log.ssc     <- log( pmax( ssc, 1 ) )
prior.fit   <- stats::lm( k.hat ~ log.ssc, subset = informative )
m           <- stats::predict( prior.fit, newdata = data.frame( log.ssc = log.ssc ) )
resid.var   <- stats::var( k.hat[ informative ] - m[ informative ] )
tau2        <- max( resid.var - mean( s2[ informative ] ), 0.05 * resid.var )

w     <- tau2 / ( tau2 + s2 )
mu.p  <- m + w * ( k.hat - m )
sd.p  <- sqrt( w * s2 )
a     <- -mu.p / sd.p
mills <- exp( stats::dnorm( a, log = TRUE ) - stats::pnorm( a, lower.tail = FALSE, log.p = TRUE ) )

k.arm <- list(
  A  = pmax( k.hat, 0 ),
  B  = k.hat,
  S1 = mu.p,
  S2 = mu.p + sd.p * mills
)
fluor.arm <- lapply( k.arm, function( k ) fit$fluorophores + ( k.hat - k ) * v.cell )

true.tau2 <- stats::var( k.true - stats::predict( stats::lm( k.true ~ log.ssc ) ) )
cat( sprintf( "Prior: tau2 estimated %.0f (SD %.0f), true residual AF variance %.0f (SD %.0f)\n",
              tau2, sqrt( tau2 ), true.tau2, sqrt( true.tau2 ) ) )
cat( sprintf( "AF variant chosen correctly for %.1f%% of cells\n", 100 * mean( af.idx == sim$truth$af.row ) ) )

# ---------------------------------------------------------------------------
# 4. Measures
# ---------------------------------------------------------------------------

# cells are grouped by their true total dye abundance, which the noise cannot
# move: with photon noise the errors of the dye and AF estimates are
# correlated, so grouping by anything estimated (s2 included) would also
# group cells by their AF error
brightness <- rowSums( truth )
decile <- cut( brightness, stats::quantile( brightness, seq( 0, 1, 0.1 ) ), include.lowest = TRUE, labels = FALSE )
rsd    <- function( x ) 1.4826 * stats::median( abs( x - stats::median( x ) ) )

af.table <- do.call( rbind, lapply( seq_len( 10 ), function( d ) {
  sel <- decile == d
  row <- data.frame( decile = d, events = sum( sel ), dye.total = signif( stats::median( brightness[ sel ] ), 3 ),
                     sd.estimate = round( sqrt( stats::median( s2[ sel ] ) ) ),
                     k.true.median = round( stats::median( k.true[ sel ] ) ),
                     calibration = round( rsd( ( k.hat[ sel ] - k.true[ sel ] ) / sqrt( s2[ sel ] ) ), 2 ) )
  for ( arm in names( k.arm ) ) {
    err <- k.arm[[ arm ]][ sel ] - k.true[ sel ]
    row[[ paste0( arm, ".bias" ) ]] <- round( mean( err ) )
    row[[ paste0( arm, ".rmse" ) ]] <- round( sqrt( mean( err^2 ) ) )
  }
  row$A.share.zero <- round( mean( k.arm$A[ sel ] == 0 ), 3 )
  row
} ) )

v.mean   <- v.library[ , 1L ]
af.like  <- names( v.mean )[ abs( v.mean ) > 0.2 ]
clamped  <- k.hat <= 0

background <- function( sel.cells ) {
  do.call( rbind, lapply( colnames( truth ), function( f ) {
    off <- sel.cells & truth[ , f ] == 0
    row <- data.frame( fluor = f, v = round( v.mean[ f ], 3 ), events = sum( off ) )
    for ( arm in names( fluor.arm ) ) {
      x <- fluor.arm[[ arm ]][ off, f ]
      row[[ paste0( arm, ".bias" ) ]] <- round( mean( x ) )
      row[[ paste0( arm, ".sd" ) ]]   <- round( stats::sd( x ) )
    }
    row
  } ) )
}
bg.all     <- background( rep( TRUE, n.cells ) )
bg.clamped <- background( clamped )

summarise.bg <- function( bg, label ) {
  out <- do.call( rbind, lapply( list( all = bg, af.like = bg[ bg$fluor %in% af.like, ] ), function( b ) {
    data.frame( t( sapply( names( k.arm ), function( arm )
      c( median.abs.bias = stats::median( abs( b[[ paste0( arm, ".bias" ) ]] ) ),
         median.sd       = stats::median( b[[ paste0( arm, ".sd" ) ]] ) ) ) ) )
  } ) )
  out$fluorophores <- rep( c( "all", "AF-like" ), each = length( k.arm ) )
  out$arm   <- rep( names( k.arm ), 2 )
  out$cells <- label
  out
}
bg.summary <- rbind( summarise.bg( bg.all, "all cells" ), summarise.bg( bg.clamped, "cells clamped by A" ) )
rownames( bg.summary ) <- NULL

positive.rmse <- sapply( fluor.arm, function( x ) {
  on <- truth > 0
  sqrt( mean( ( x[ on ] - truth[ on ] )^2 ) )
} )

# ---------------------------------------------------------------------------
# 5. Output
# ---------------------------------------------------------------------------

cat( "\nAF abundance by decile of true total dye abundance\n" )
print( af.table, row.names = FALSE )
cat( "\nBackground in fluorophores on cells where they are off: median |bias| and median SD\n" )
print( bg.summary[ , c( "cells", "fluorophores", "arm", "median.abs.bias", "median.sd" ) ], row.names = FALSE )
cat( "\nRMSE of fluorophores on cells where they are on, by arm\n" )
print( round( positive.rmse ) )
cat( "\nAF-like fluorophores, all cells\n" )
print( bg.all[ bg.all$fluor %in% af.like, ], row.names = FALSE )

utils::write.csv( af.table,   file.path( out.dir, "af_by_brightness_decile.csv" ), row.names = FALSE )
utils::write.csv( bg.all,     file.path( out.dir, "background_all_cells.csv" ), row.names = FALSE )
utils::write.csv( bg.clamped, file.path( out.dir, "background_clamped_cells.csv" ), row.names = FALSE )
utils::write.csv( bg.summary, file.path( out.dir, "background_summary.csv" ), row.names = FALSE )
cat( "\nTables written to", out.dir, "\n" )
