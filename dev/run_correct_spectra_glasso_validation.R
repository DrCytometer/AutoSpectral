# run_correct_spectra_glasso_validation.R
#
# Recovery test for correct.spectra.glasso(): unlike
# run_fix_my_unmix_synthetic_null.R (which measures an estimator's own bias
# on data with zero injected library error), this script injects a KNOWN,
# quantified spectral error into a copy of a real reference spectra table,
# simulates a fully-stained sample from the TRUE (uncorrupted) spectra via
# sim.flow.data(), and checks whether correct.spectra.glasso() - started from
# the corrupted spectra - recovers both the injected spillover coefficient
# and the row shape it came from.
#
# Optionally runs fix.my.unmix() on the identical synthetic pair with the
# identical corrupted starting spectra, so the two estimators can be compared
# head to head on data where the ground truth is known exactly rather than
# only on real data where it is not.
#
# READING THE RECOVERY CHECK - the part of this script that was wrong before:
# corrupting fluorophore j's row with eps of fluorophore k's shape
# (spectra.wrong[j,] <- normalize(spectra[j,] + eps*spectra[k,])) does not
# show up as a positive coefficient of k on j (k's shape "leaking into" j's
# channel). Unmixing is a change of basis, and the two rows actually coupled
# by that change are the other way round: j's own negative population,
# regressed on k, reads close to -eps, because j's true photons need a
# little of k's channel subtracted to be represented correctly in the
# corrupted basis. So the pair check below reads spillover[target, source]
# (row = target, the corrupted row; column = source, the donor), expects it
# near -eps, and separately reports spillover[source, target] for contrast,
# which should stay near zero. Both correct.spectra.glasso() and
# fix.my.unmix() share this identification strategy and this same reading
# convention. See CONTEXT_information_sources.md section 4 (hypernegativity)
# for why a negative coefficient is the unambiguous one, and the "How the
# lasso row is estimated" section of correct_spectra_glasso.R for the change-
# of-basis argument in full.
#
# PREREQUISITES:
#   devtools::load_all() -- this script calls .sim.scatter() and
#     .sim.detector.preset(), non-exported helpers from
#     simulate_flow_data.R (same prerequisite as
#     run_fix_my_unmix_synthetic_null.R).
#   asp, flow.control -- already built for the target cytometer/panel.
#   gv.spectra -- a real reference spectra matrix (fluorophores x detectors,
#     L-infinity normalised), treated as ground truth. The AF row, if
#     present, is dropped: this test is deliberately AF-free so a failure to
#     recover the injected error cannot be attributed to unmodelled
#     background, the same caveat run_signature_correction_recovery_test.R
#     flags for its own real-data recovery test.
#
# Run top to bottom. Inspect gv.pair.recovery and gv.spectra.recovery at the
# end; if gv.compare.fix.my.unmix is TRUE, gv.pair.recovery and
# gv.spectra.recovery each carry both estimators' columns side by side.

# ---------------------------------------------------------------------------
# 0. Setup - EDIT THIS SECTION
# ---------------------------------------------------------------------------

asp.dir <- "/Users/oliverburton/Bioinformatics/AutoSpectral"
asp.rcpp.dir <- "/Users/oliverburton/Bioinformatics/AutoSpectralRcpp"
devtools::load_all(asp.rcpp.dir)
devtools::load_all(asp.dir)

asp <- get.autospectral.param( cytometer = "discover" )   # set to your cytometer

gv.spectra.file <- "A8 Cells_autospectral_spectra.csv"               # read.spectra() input
gv.spectra.dir  <- "./table_spectra"

gv.af.name <- "AF"

gv.spectra <- read.spectra( spectra.file = gv.spectra.file,
                            spectra.dir  = gv.spectra.dir,
                            remove.af    = TRUE )         # ground truth, panel only

flow.control <- reload.flow.control(
  control.dir      = "./Plate_Cells",
  control.def.file = "./Plate_Cells/fcs_control_file.csv",
  asp              = asp
)

# Injected corruptions: each row is one (target, source, eps) triple. Target
# j's starting spectrum is set to normalize( true_j + eps * true_k ), so the
# corruption's size in the small-error linear regime `fix.my.unmix()`'s own
# documentation assumes is exactly `eps`. Keep `eps` comfortably under
# `max.coefficient` (0.5 by default in both estimators) - this test is about
# whether a plausible library error is recovered, not about the boundary of
# what either estimator will accept.
gv.injections <- data.frame(
  target = c( "BUV805", "BV480" ),
  source = c( "BUV395", "BV510" ),
  eps    = c( 0.12,     0.08 ),
  stringsAsFactors = FALSE
)

gv.n.cells    <- 100000L
gv.complexity <- 0.1
gv.seed       <- 42L
gv.output.dir <- "./correct_spectra_glasso_validation"

gv.compare.fix.my.unmix <- TRUE

if ( !dir.exists( gv.output.dir ) ) dir.create( gv.output.dir, recursive = TRUE )

# Arguments forwarded to correct.spectra.glasso(). AF is off for this test
# (bg.mode = "none") to isolate whether the estimator recovers a panel-only
# injected error, per the caveat above.
gv.glasso.args <- list(
  bg.mode        = "none",
  large.gate     = TRUE,
  max.iter       = 5L,
  downsample     = FALSE,     # already sized via gv.n.cells
  update.spectra = TRUE,
  figures        = FALSE,
  save           = FALSE,
  verbose        = TRUE
)

# Matching arguments for fix.my.unmix(), where the two functions' parameter
# names differ (bg.mode/large.gate/downsample/update.spectra/figures/save
# are shared and reused as-is).
gv.fix.args <- c(
  gv.glasso.args,
  list( estimator = "truncated", min.bin.events = 50 )
)

# ---------------------------------------------------------------------------
# 1. Validation
# ---------------------------------------------------------------------------

if ( !setequal( colnames( gv.spectra ), flow.control$spectral.channel ) )
  stop( "gv.spectra column names must match flow.control$spectral.channel.",
        call. = FALSE )

if ( !all( c( gv.injections$target, gv.injections$source ) %in%
           rownames( gv.spectra ) ) )
  stop( "Every `target`/`source` in gv.injections must be a row of gv.spectra.",
        call. = FALSE )

fluorophores <- rownames( gv.spectra )

# ---------------------------------------------------------------------------
# 2. Inject the corruption
# ---------------------------------------------------------------------------

gv.spectra.wrong <- gv.spectra

for ( r in seq_len( nrow( gv.injections ) ) ) {
  
  j   <- gv.injections$target[ r ]
  k   <- gv.injections$source[ r ]
  eps <- gv.injections$eps[ r ]
  
  corrupted <- gv.spectra[ j, ] + eps * gv.spectra[ k, ]
  gv.spectra.wrong[ j, ] <- corrupted / max( corrupted )
}

# ---------------------------------------------------------------------------
# 3. Realistic scatter and the FCS writer, matching
#    run_fix_my_unmix_synthetic_null.R's own pattern
# ---------------------------------------------------------------------------

time.name <- asp$default.time.parameter

.gv.write.fcs <- function( mat, asp, file.name, output.dir ) {
  
  n.par <- ncol( mat )
  keys  <- list()
  
  for ( i in seq_len( n.par ) ) {
    keys[[ paste0( "$P", i, "N" ) ]] <- colnames( mat )[ i ]
    keys[[ paste0( "$P", i, "B" ) ]] <- "32"
    keys[[ paste0( "$P", i, "E" ) ]] <- "0,0"
    keys[[ paste0( "$P", i, "R" ) ]] <- as.character( ceiling(
      max( asp$expr.data.max, max( mat[ , i ], na.rm = TRUE ) + 1 ) ) )
  }
  
  keys[[ "$CYT" ]] <- asp$cytometer
  
  writeFCS( mat = mat, keys = keys, file.name = file.name,
            output.dir = output.dir )
  
  file.path( output.dir, file.name )
}

# A true unstained control (no AF, no fluorophore expression) is outside
# sim.flow.data()'s own generative model (it requires complexity > 0), so
# this reuses its detector-noise helper directly on an all-zero signal,
# exactly as run_fix_my_unmix_synthetic_null.R's .sn.sim.unstained() does
# with af.spectra = NULL.
.gv.sim.unstained <- function( n.cells, spectra, asp, seed ) {
  
  set.seed( seed )
  
  det.names   <- colnames( spectra )
  n.detectors <- ncol( spectra )
  
  detector.preset <- .sim.detector.preset( asp$cytometer )
  cpu             <- detector.preset$counts.per.unit
  
  signal <- matrix(
    stats::rpois( n.cells * n.detectors, lambda = 0 ),
    nrow = n.cells, ncol = n.detectors ) / cpu
  
  readout.noise <- matrix(
    stats::rnorm( n.cells * n.detectors, mean = 0,
                  sd = detector.preset$readout.sd ),
    nrow = n.cells, ncol = n.detectors )
  
  signal <- signal + readout.noise - detector.preset$dark.offset
  colnames( signal ) <- det.names
  
  signal
}

gv.fs.scatter <- .sim.scatter( n.cells = gv.n.cells, scatter.data = NULL,
                               fsc.mean.log = log( 200000 ), fsc.sd.log = 0.4,
                               ssc.mean.log = log( 80000 ),  ssc.sd.log = 0.5,
                               fsc.ssc.cor  = 0.6 )
gv.un.scatter <- .sim.scatter( n.cells = gv.n.cells, scatter.data = NULL,
                               fsc.mean.log = log( 200000 ), fsc.sd.log = 0.4,
                               ssc.mean.log = log( 80000 ),  ssc.sd.log = 0.5,
                               fsc.ssc.cor  = 0.6 )

# .sim.scatter() always names its output c("FSC","SSC") regardless of the
# instrument; rename to the real scatter channel names so the synthetic
# matrices match flow.control$scatter.and.channel.spectral downstream.
colnames( gv.fs.scatter ) <- flow.control$scatter.parameter
colnames( gv.un.scatter ) <- flow.control$scatter.parameter

# ---------------------------------------------------------------------------
# 4. Generate the synthetic pair from the TRUE spectra and write it out
# ---------------------------------------------------------------------------

gv.fs <- sim.flow.data(
  spectra            = gv.spectra,
  asp                = asp,
  n.cells            = gv.n.cells,
  complexity         = gv.complexity,
  af.spectra         = NULL,
  af.variation       = FALSE,
  spectral.variation = FALSE,
  scatter.data       = gv.fs.scatter[ , 1:2, drop = FALSE ],
  seed               = gv.seed
)

gv.un.raw <- .gv.sim.unstained(
  n.cells = gv.n.cells, spectra = gv.spectra, asp = asp, seed = gv.seed + 1L )

gv.fs.mat <- cbind(
  matrix( seq_len( gv.n.cells ), ncol = 1, dimnames = list( NULL, time.name ) ),
  gv.fs.scatter,
  gv.fs$raw[ , colnames( gv.spectra ), drop = FALSE ] )

gv.un.mat <- cbind(
  matrix( seq_len( gv.n.cells ), ncol = 1, dimnames = list( NULL, time.name ) ),
  gv.un.scatter,
  gv.un.raw[ , colnames( gv.spectra ), drop = FALSE ] )

if ( !setequal( colnames( gv.spectra ), flow.control$spectral.channel ) )
  stop( "gv.spectra column names must match flow.control$spectral.channel.",
        call. = FALSE )

if ( length( flow.control$scatter.parameter ) < 2L )
  stop( "flow.control$scatter.parameter needs at least two entries.",
        call. = FALSE )

gv.fs.path <- .gv.write.fcs( gv.fs.mat, asp, "synthetic_fully_stained.fcs",
                             gv.output.dir )
gv.un.path <- .gv.write.fcs( gv.un.mat, asp, "synthetic_unstained.fcs",
                             gv.output.dir )

# ---------------------------------------------------------------------------
# 5. Run correct.spectra.glasso() (and, optionally, fix.my.unmix()) starting
#    from the corrupted spectra
# ---------------------------------------------------------------------------

cat( "\n===================== correct.spectra.glasso() =====================\n" )

gv.glasso.fit <- do.call( correct.spectra.glasso, c(
  list(
    spectra              = gv.spectra.wrong,
    unstained.sample     = gv.un.path,
    fully.stained.sample = gv.fs.path,
    flow.control         = flow.control,
    asp                  = asp,
    af.name              = gv.af.name
  ),
  gv.glasso.args
) )

gv.fix.fit <- NULL

if ( gv.compare.fix.my.unmix ) {
  
  cat( "\n===================== fix.my.unmix() =====================\n" )
  
  gv.fix.fit <- do.call( fix.my.unmix, c(
    list(
      spectra              = gv.spectra.wrong,
      unstained.sample     = gv.un.path,
      fully.stained.sample = gv.fs.path,
      flow.control         = flow.control,
      asp                  = asp,
      af.name              = gv.af.name
    ),
    gv.fix.args
  ) )
}

# ---------------------------------------------------------------------------
# 6. Recovery comparison
# ---------------------------------------------------------------------------

.gv.angle <- function( a, b ) {
  cosang <- sum( a * b ) / ( sqrt( sum( a^2 ) ) * sqrt( sum( b^2 ) ) )
  180 / pi * acos( pmin( 1, pmax( -1, cosang ) ) )
}

# For each injected pair: the true coefficient is `eps`, and it belongs to
# spillover[target, source] (row = the corrupted row, column = the donor),
# expected close to -eps - see the header comment for why. spillover[source,
# target] is reported alongside it purely for contrast: this is the cell the
# old version of this script checked, and it should sit near zero, which is
# not a failure to recover anything, it is the correct answer to a different
# question than the one that matters here.
gv.pair.recovery <- gv.injections

gv.pair.recovery$recovered.glasso <- mapply(
  function( j, k ) gv.glasso.fit$spillover[ j, k ],
  gv.injections$target, gv.injections$source )
gv.pair.recovery$other.cell.glasso <- mapply(
  function( j, k ) gv.glasso.fit$spillover[ k, j ],
  gv.injections$target, gv.injections$source )

# Direct check on the mechanism the fix targets: does the corrupted target's
# active set actually contain its true donor now? This is the thing that
# determines whether extract.raw.signature() gets a chance to correct the
# row at all, independent of how the spillover matrix itself reads.
gv.pair.recovery$donor.in.active.set <- mapply(
  function( j, k ) k %in% gv.glasso.fit$active.set[[ j ]],
  gv.injections$target, gv.injections$source )

if ( !is.null( gv.fix.fit ) ) {
  
  gv.pair.recovery$recovered.fix <- mapply(
    function( j, k ) gv.fix.fit$spillover[ j, k ],
    gv.injections$target, gv.injections$source )
  gv.pair.recovery$other.cell.fix <- mapply(
    function( j, k ) gv.fix.fit$spillover[ k, j ],
    gv.injections$target, gv.injections$source )
  
  # fix.my.unmix() always hands extract.raw.signature() the whole panel
  # (active = fluorophores, ridge-penalised), never a lasso-selected
  # subset, so there is no active-set gap for it to have in the first
  # place - this column exists only so the two estimators' recovery numbers
  # sit side by side without a silent asymmetry in what was checked.
  gv.pair.recovery$donor.in.active.set.fix <- TRUE
}

# For every fluorophore's row, not only the corrupted ones: an estimator
# that recovers the injected pairs while quietly distorting everything else
# has not actually won anything. This check does not depend on which cell
# of the spillover matrix carries the signal, so it is unchanged from
# before and remains the primary readout of whether the fix actually works
# end to end.
gv.spectra.recovery <- data.frame(
  fluorophore   = fluorophores,
  deg.start     = vapply( fluorophores, function( j )
    .gv.angle( gv.spectra.wrong[ j, ], gv.spectra[ j, ] ), numeric( 1 ) ),
  deg.glasso    = vapply( fluorophores, function( j )
    .gv.angle( gv.glasso.fit$spectra[ j, ], gv.spectra[ j, ] ), numeric( 1 ) ),
  stringsAsFactors = FALSE
)

if ( !is.null( gv.fix.fit ) )
  gv.spectra.recovery$deg.fix <- vapply( fluorophores, function( j )
    .gv.angle( gv.fix.fit$spectra[ j, ], gv.spectra[ j, ] ), numeric( 1 ) )

utils::write.csv( gv.pair.recovery,
                  file.path( gv.output.dir, "gv_pair_recovery.csv" ),
                  row.names = FALSE )
utils::write.csv( gv.spectra.recovery,
                  file.path( gv.output.dir, "gv_spectra_recovery.csv" ),
                  row.names = FALSE )

cat( "\nInjected-pair recovery (want recovered.* close to -eps;",
     "other.cell.* should stay near zero; donor.in.active.set.* should be",
     "TRUE):\n" )
print( gv.pair.recovery, digits = 4 )

cat( "\nPer-fluorophore angle to ground truth, degrees",
     "(want deg.* well below deg.start for injected rows,",
     "and not meaningfully above deg.start elsewhere):\n" )
print( gv.spectra.recovery, digits = 4 )

cat( sprintf(
  "\nMean angle across the whole panel: start %.3f, glasso %.3f%s\n",
  mean( gv.spectra.recovery$deg.start ),
  mean( gv.spectra.recovery$deg.glasso ),
  if ( !is.null( gv.fix.fit ) )
    sprintf( ", fix.my.unmix %.3f", mean( gv.spectra.recovery$deg.fix ) )
  else "" ) )

cat( "\nReport back these files from", gv.output.dir, ":\n" )
cat( "  gv_pair_recovery.csv\n" )
cat( "  gv_spectra_recovery.csv\n" )

print( gv.glasso.fit$convergence.log )
print( gv.glasso.fit$lambda.log )