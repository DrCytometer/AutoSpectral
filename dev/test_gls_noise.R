# dev/test_gls_noise.R
#
# Validation harness for estimate.noise.model() and unmix.gls().
# Run interactively, section by section. Each section prints a PASS/FAIL
# line and the numbers behind it.

asp.dir <- "/Users/oliverburton/Bioinformatics/AutoSpectral"
asp.rcpp.dir <- "/Users/oliverburton/Bioinformatics/AutoSpectralRcpp"
devtools::load_all(asp.rcpp.dir)
devtools::load_all(asp.dir)

set.seed( 42 )

# =============================================================================
# 0. Setup
# =============================================================================

asp <- get.autospectral.param( cytometer = "Aurora" )

ref.lib <- read.csv(
  system.file( "extdata", "Aurora_spectral_reference_library.csv",
               package = "AutoSpectral" ),
  row.names = 1, check.names = FALSE
)

panel   <- c( "BUV395", "BUV496", "BUV661", "BUV737", "BV421", "BV510",
              "BV605", "BV650", "BV711", "BV785", "FITC", "PerCP",
              "PE", "PE-Dazzle594", "PE-Cy7", "APC", "AF700", "APC-Fire750" )
panel   <- intersect( panel, rownames( ref.lib ) )
spectra <- as.matrix( ref.lib[ panel, , drop = FALSE ] )
spectra <- spectra / apply( abs( spectra ), 1, max )

det.n   <- ncol( spectra )
fluor.n <- nrow( spectra )

cat( sprintf( "Panel: %d fluorophores, %d detectors\n", fluor.n, det.n ) )

# --- synthetic AF dictionary: smooth, red-shifted, correlated components -----
make.af.dictionary <- function( n.af, det.n, seed = 1 ) {
  set.seed( seed )
  centres <- seq( 0.25, 0.8, length.out = n.af ) * det.n
  af <- t( sapply( centres, function( cc )
    exp( -0.5 * ( ( seq_len( det.n ) - cc ) / ( det.n / 6 ) )^2 ) ) )
  af <- af + matrix( runif( n.af * det.n, 0, 0.08 ), nrow = n.af )
  af <- af / apply( af, 1, max )
  rownames( af ) <- paste0( "AF", seq_len( n.af ) )
  colnames( af ) <- colnames( spectra )
  af
}
af.spectra <- make.af.dictionary( 12L, det.n )

# --- synthetic variants: smooth low-rank perturbations ----------------------
make.variants <- function( spectra, n.var = 8L, scale = 0.03, seed = 2 ) {
  set.seed( seed )
  det.n <- ncol( spectra )
  smooth.bump <- function( centre, width )
    exp( -0.5 * ( ( seq_len( det.n ) - centre ) / width )^2 )
  out <- lapply( rownames( spectra ), function( fl ) {
    ref  <- spectra[ fl, ]
    peak <- which.max( ref )
    b1   <- smooth.bump( peak + 3, 4 )
    b2   <- smooth.bump( peak - 3, 4 )
    coef <- matrix( rnorm( n.var * 2, sd = scale ), ncol = 2 )
    v <- t( apply( coef, 1, function( cf ) ref + cf[ 1 ] * b1 + cf[ 2 ] * b2 ) )
    v <- pmax( v, 0 )
    v <- v / apply( v, 1, max )
    colnames( v ) <- colnames( spectra )
    v
  } )
  names( out ) <- rownames( spectra )
  out
}
variants   <- make.variants( spectra )
delta.list <- lapply( names( variants ), function( fl )
  sweep( variants[[ fl ]], 2L, spectra[ fl, ], "-" ) )
names( delta.list ) <- names( variants )

variant.basis <- build.variant.basis( delta.list, rank = 2L )
cat( sprintf( "Variant basis built for %d/%d fluorophores\n",
              length( variant.basis ), fluor.n ) )


# =============================================================================
# 1. Noise model recovery
# =============================================================================
# Simulate with the AF dictionary but WITHOUT spectral variation, so that
# the residual after projection is (almost) pure instrument noise. This is
# the clean case; section 1b repeats it with variation present.

sim.clean <- sim.flow.data(
  spectra            = spectra,
  asp                = asp,
  n.cells            = 30000L,
  af.spectra         = af.spectra,
  variants           = NULL,
  spectral.variation = FALSE,
  shot.noise         = TRUE,
  spillover.noise    = TRUE,
  detector.noise     = TRUE,
  seed               = 7
)

truth.kappa <- sim.clean$params$detector.preset$counts.per.unit
truth.sd    <- sim.clean$params$detector.preset$readout.sd

nm <- estimate.noise.model(
  raw.data   = sim.clean$raw,
  spectra    = spectra,
  af.spectra = af.spectra
)

kappa.ratio <- nm$counts.per.unit / truth.kappa
sd.ratio    <- sqrt( nm$read.var ) / truth.sd

cat( "\n--- 1. Noise model recovery ---\n" )
cat( sprintf( "counts.per.unit ratio (est/true): median %.3f, IQR %.3f-%.3f\n",
              median( kappa.ratio ),
              quantile( kappa.ratio, 0.25 ), quantile( kappa.ratio, 0.75 ) ) )
cat( sprintf( "read SD ratio (est/true):         median %.3f, IQR %.3f-%.3f\n",
              median( sd.ratio ),
              quantile( sd.ratio, 0.25 ), quantile( sd.ratio, 0.75 ) ) )
cat( sprintf( "median per-detector R2:           %.3f\n",
              median( nm$r.squared, na.rm = TRUE ) ) )
cat( sprintf( "PASS if both ratios in [0.7, 1.4] and R2 > 0.8: %s\n",
              if ( median( kappa.ratio ) > 0.7 && median( kappa.ratio ) < 1.4 &&
                   median( sd.ratio )    > 0.7 && median( sd.ratio )    < 1.4 &&
                   median( nm$r.squared, na.rm = TRUE ) > 0.8 ) "PASS" else "FAIL" ) )

# QC plot: mean-variance relationship
if ( interactive() ) {
  ft <- nm$fit.table
  plot( ft$mean, ft$variance, pch = 16, cex = 0.3, log = "xy",
        xlab = "fitted signal", ylab = "residual variance",
        main = "Mean-variance, all detectors" )
}


# =============================================================================
# 1b. Noise model under spectral variation (contamination check)
# =============================================================================

sim.var <- sim.flow.data(
  spectra            = spectra,
  asp                = asp,
  n.cells            = 30000L,
  af.spectra         = af.spectra,
  variants           = variants,
  spectral.variation = TRUE,
  seed               = 8
)

nm.var <- estimate.noise.model( sim.var$raw, spectra, af.spectra,
                                verbose = FALSE )

cat( "\n--- 1b. Noise model with spectral variation present ---\n" )
cat( sprintf( "kappa inflation vs clean fit: median %.2fx\n",
              median( nm$counts.per.unit / nm.var$counts.per.unit ) ) )
cat( "Expect kappa to DROP (variance up => slope up => kappa down).\n" )
cat( "This is the contamination the variant covariance term is meant to absorb.\n" )


# =============================================================================
# 2. Woodbury vs dense equivalence
# =============================================================================

sub <- sim.var$raw[ 1:200, , drop = FALSE ]
af.i <- sim.var$truth$af.row[ 1:200 ]

x.wb <- unmix.gls( sub, spectra, nm, variant.basis, af.spectra, af.i,
                   method = "woodbury", verbose = FALSE )
x.dn <- unmix.gls( sub, spectra, nm, variant.basis, af.spectra, af.i,
                   method = "dense",    verbose = FALSE )

rel.err <- max( abs( x.wb - x.dn ) ) / max( abs( x.dn ) )
cat( "\n--- 2. Woodbury vs dense ---\n" )
cat( sprintf( "max relative difference: %.3e  (%s)\n", rel.err,
              if ( rel.err < 1e-8 ) "PASS" else "FAIL" ) )


# =============================================================================
# 3. Accuracy and calibration
# =============================================================================
# Give OLS/WLS the SAME per-cell AF as GLS, so the comparison isolates the
# covariance model rather than re-measuring AF assignment.

x.ols <- matrix( NA_real_, n.eval, fluor.n + 1L )
x.wls <- x.ols
for ( g in unique( af.i.eval ) ) {
  idx <- which( af.i.eval == g )
  S.g <- rbind( spectra, AF = af.spectra[ g, ] )
  x.ols[ idx, ] <- unmix.ols.fast( y.eval[ idx, , drop = FALSE ], S.g )
  x.wls[ idx, ] <- unmix.wls.fast( y.eval[ idx, , drop = FALSE ], S.g, w.diag )
}
colnames( x.ols ) <- colnames( x.wls ) <- c( rownames( spectra ), "AF" )

# --- RMSE, stratified by whether the fluorophore is actually present -------
strat.rmse <- function( est, tru, fl, on ) {
  idx <- if ( on ) which( tru[ , fl ] > 0 ) else which( tru[ , fl ] == 0 )
  sqrt( mean( ( est[ idx, fl ] - tru[ idx, fl ] )^2 ) )
}

tab <- t( sapply( common, function( fl ) c(
  OLS.on  = strat.rmse( x.ols, truth, fl, TRUE  ),
  GLS.on  = strat.rmse( x.gls, truth, fl, TRUE  ),
  OLS.off = strat.rmse( x.ols, truth, fl, FALSE ),
  GLS.off = strat.rmse( x.gls, truth, fl, FALSE )
) ) )
print( round( tab, 1 ) )
cat( sprintf( "median GLS/OLS, ON cells:  %.3f\n",
              median( tab[ , "GLS.on" ]  / tab[ , "OLS.on" ]  ) ) )
cat( sprintf( "median GLS/OLS, OFF cells: %.3f\n",
              median( tab[ , "GLS.off" ] / tab[ , "OLS.off" ] ) ) )

# --- calibration: per-cell standardised residual ---------------------------
# This is the real CRLB test. z should be ~N(0,1) if Sigma is right.

z <- ( x.gls[ , common, drop = FALSE ] - truth[ , common, drop = FALSE ] ) /
  gls$se[ , common, drop = FALSE ]

z.summary <- t( sapply( common, function( fl ) {
  on  <- truth[ , fl ] > 0
  c( sd.on  = stats::sd( z[ on,  fl ] ),
     sd.off = stats::sd( z[ !on, fl ] ) )
} ) )
print( round( z.summary, 2 ) )
cat( sprintf( "\nmedian |z| SD, ON:  %.2f\nmedian |z| SD, OFF: %.2f\n",
              median( z.summary[ , "sd.on" ] ),
              median( z.summary[ , "sd.off" ] ) ) )
cat( "Want ~1.0 in both. >1 = Sigma too small (missing a variance source).\n" )
cat( "<1 = Sigma too large, or the fit is using information it shouldn't.\n" )

# =============================================================================
# 3. Accuracy vs the Cramer-Rao bound
# =============================================================================

n.eval <- 5000L
eval.i <- seq_len( n.eval )

y.eval    <- sim.var$raw[ eval.i, , drop = FALSE ]
af.i.eval <- sim.var$truth$af.row[ eval.i ]
truth     <- sim.var$truth$abundances[ eval.i, , drop = FALSE ]

w.diag <- 1 / pmax( abs( colMeans( y.eval ) ), 1e-6 )

S.af <- rbind( spectra, AF = colMeans( af.spectra ) )

x.ols <- unmix.ols.fast( y.eval, S.af )
x.wls <- unmix.wls.fast( y.eval, S.af, w.diag )
gls   <- unmix.gls( y.eval, spectra, nm, variant.basis,
                    af.spectra, af.i.eval,
                    n.iter = 2L, return.variance = TRUE, verbose = FALSE )
x.gls <- gls$unmixed

rmse <- function( est, tru, cols )
  sqrt( colMeans( ( est[ , cols, drop = FALSE ] -
                      tru[ , cols, drop = FALSE ] )^2 ) )

common <- intersect( colnames( x.gls ), colnames( truth ) )
common <- setdiff( common, "AF" )

r.ols <- rmse( x.ols, truth, common )
r.wls <- rmse( x.wls, truth, common )
r.gls <- rmse( x.gls, truth, common )

# empirical CRLB: the model's own predicted SE, evaluated at the fit
crlb <- colMeans( gls$se[ , common, drop = FALSE ] )

cat( "\n--- 3. Accuracy vs CRLB ---\n" )
print( round( data.frame(
  OLS        = r.ols,
  WLS        = r.wls,
  GLS        = r.gls,
  CRLB       = crlb,
  GLS.over.CRLB = r.gls / crlb,
  GLS.over.OLS  = r.gls / r.ols
), 3 ) )
cat( sprintf( "\nmedian GLS/OLS RMSE ratio: %.3f (want < 1)\n",
              median( r.gls / r.ols ) ) )
cat( sprintf( "median GLS/CRLB ratio:     %.3f (want ~1, ALARM if < 0.9)\n",
              median( r.gls / crlb ) ) )


# =============================================================================
# 4. Overfitting canary: negative-population spread
# =============================================================================
# For each fluorophore, take the cells where its true abundance is zero.
# Their unmixed spread should MATCH the model's predicted SE, not beat it.

canary <- sapply( common, function( fl ) {
  neg <- which( truth[ , fl ] == 0 )
  if ( length( neg ) < 200L ) return( c( NA, NA, NA ) )
  obs  <- mad( x.gls[ neg, fl ] )
  pred <- median( gls$se[ neg, fl ] )
  c( observed = obs, predicted = pred, ratio = obs / pred )
} )

cat( "\n--- 4. Negative-population spread (overfitting canary) ---\n" )
print( round( t( canary ), 3 ) )
cat( sprintf( "median observed/predicted: %.3f\n",
              median( canary[ "ratio", ], na.rm = TRUE ) ) )
cat( "Ratio << 1 means the fit is absorbing real signal. Ratio >> 1 means\n" )
cat( "the noise model is too optimistic (usually kappa too high).\n" )


# =============================================================================
# 5. Phantom fluorophore test
# =============================================================================
# Add a fluorophore to the unmixing matrix that is NOT in the data.
# Its abundance must stay at zero.

phantom.name <- setdiff( rownames( ref.lib ), panel )[ 1 ]
spectra.plus <- rbind( spectra, as.matrix( ref.lib[ phantom.name, , drop = FALSE ] ) )
spectra.plus <- spectra.plus / apply( abs( spectra.plus ), 1, max )

vb.plus <- variant.basis   # phantom has no variants, correctly omitted

x.phantom <- unmix.gls( y.eval, spectra.plus, nm, vb.plus,
                        af.spectra, af.i.eval, verbose = FALSE )

real.scale <- median( apply( truth[ , common, drop = FALSE ], 2,
                             function( z ) mad( z[ z > 0 ] ) ), na.rm = TRUE )

cat( "\n--- 5. Phantom fluorophore ---\n" )
cat( sprintf( "phantom: %s\n", phantom.name ) )
cat( sprintf( "median |phantom| = %.1f, 99th pct = %.1f, real-signal scale = %.1f\n",
              median( abs( x.phantom[ , phantom.name ] ) ),
              quantile( abs( x.phantom[ , phantom.name ] ), 0.99 ),
              real.scale ) )
cat( sprintf( "phantom 99th pct / real scale: %.3f (want << 0.1)\n",
              quantile( abs( x.phantom[ , phantom.name ] ), 0.99 ) / real.scale ) )


# =============================================================================
# 6. Held-out detector cross-validation (works on real data too)
# =============================================================================
# Hide k detectors, fit on the rest, predict the hidden ones.
# This is the model-selection criterion: use it to choose variant rank,
# AF dictionary size, and whether a term earns its parameters.

heldout.cv <- function( raw, spectra, nm, variant.basis, af.spectra, af.index,
                        n.folds = 4L, ... ) {
  det.names <- colnames( spectra )
  folds <- split( sample( det.names ), rep( seq_len( n.folds ),
                                            length.out = length( det.names ) ) )
  err <- sapply( folds, function( hidden ) {
    keep <- setdiff( det.names, hidden )
    nm.k <- list(
      read.var        = nm$read.var[ keep ],
      counts.per.unit = nm$counts.per.unit[ keep ],
      kappa.pooled    = nm$kappa.pooled
    )
    vb.k <- lapply( variant.basis, function( vb ) list(
      basis  = vb$basis[ , keep, drop = FALSE ],
      lambda = vb$lambda ) )
    x <- unmix.gls( raw[ , keep, drop = FALSE ], spectra[ , keep, drop = FALSE ],
                    nm.k, vb.k, af.spectra[ , keep, drop = FALSE ], af.index,
                    verbose = FALSE, ... )
    S.full <- rbind( spectra, AF = 0 )
    pred <- sapply( hidden, function( h ) {
      s.h <- c( spectra[ , h ], 0 )
      as.vector( x %*% c( spectra[ , h ], NA ) )  # placeholder, see note
    } )
    NULL
  } )
  err
}

# NOTE: the AF row differs per cell, so prediction of the hidden detectors
# must reconstruct it per cell. Simpler explicit version:

cv.heldout <- function( raw, spectra, nm, variant.basis, af.spectra, af.index,
                        hidden, ... ) {
  keep <- setdiff( colnames( spectra ), hidden )
  nm.k <- list( read.var        = nm$read.var[ keep ],
                counts.per.unit = nm$counts.per.unit[ keep ],
                kappa.pooled    = nm$kappa.pooled )
  vb.k <- lapply( variant.basis, function( vb )
    list( basis = vb$basis[ , keep, drop = FALSE ], lambda = vb$lambda ) )
  
  x <- unmix.gls( raw[ , keep, drop = FALSE ], spectra[ , keep, drop = FALSE ],
                  nm.k, vb.k, af.spectra[ , keep, drop = FALSE ], af.index,
                  verbose = FALSE, ... )
  
  S.hidden  <- spectra[ , hidden, drop = FALSE ]                # F x H
  af.hidden <- af.spectra[ af.index, hidden, drop = FALSE ]     # N x H
  
  pred <- x[ , rownames( spectra ), drop = FALSE ] %*% S.hidden +
    x[ , "AF" ] * af.hidden
  
  sqrt( mean( ( raw[ , hidden, drop = FALSE ] - pred )^2 ) )
}

set.seed( 11 )
folds <- split( sample( colnames( spectra ) ), rep( 1:4, length.out = det.n ) )

cv.with    <- mean( sapply( folds, function( h )
  cv.heldout( y.eval, spectra, nm, variant.basis, af.spectra, af.i.eval, h ) ) )
cv.without <- mean( sapply( folds, function( h )
  cv.heldout( y.eval, spectra, nm, list(), af.spectra, af.i.eval, h ) ) )

cat( "\n--- 6. Held-out detector CV ---\n" )
cat( sprintf( "RMSE with variant covariance:    %.2f\n", cv.with ) )
cat( sprintf( "RMSE without variant covariance: %.2f\n", cv.without ) )
cat( sprintf( "improvement: %.1f%% (%s)\n",
              100 * ( 1 - cv.with / cv.without ),
              if ( cv.with < cv.without ) "term earns its place" else
                "term is NOT helping" ) )


# =============================================================================
# 7. Timing
# =============================================================================

t.ols <- system.time( unmix.ols.fast( y.eval, S.af ) )[ "elapsed" ]
t.gls <- system.time( unmix.gls( y.eval, spectra, nm, variant.basis,
                                 af.spectra, af.i.eval, verbose = FALSE ) )[ "elapsed" ]

cat( "\n--- 7. Timing (R prototype, single thread) ---\n" )
cat( sprintf( "OLS: %.2f s, GLS: %.2f s for %d cells (%.0fx)\n",
              t.ols, t.gls, n.eval, t.gls / max( t.ols, 1e-3 ) ) )
cat( "The R loop is the bottleneck, not the algebra. Compare the Woodbury\n" )
cat( "rank against D before concluding anything about the C++ port:\n" )
cat( sprintf( "  detectors = %d, typical active set implies rank ~ %d\n",
              det.n, round( 3 * fluor.n / 3 ) ) )