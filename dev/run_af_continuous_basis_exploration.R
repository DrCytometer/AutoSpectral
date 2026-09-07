# dev/run_af_continuous_basis_exploration.R
#
# Decides whether the discrete AF library can be replaced by a continuous
# low-dimensional basis, and if so how many components it needs.
#
#   1. Out-of-span scree of the AF library
#   2. Quantisation error the discrete path currently accepts
#   3. Numerical equivalence of unmix.af.basis() against explicit joint OLS
#   4. Negative-AF fraction of the unconstrained continuous fit
#   5. Shortlist recall for the hybrid variant
#   6. Timing: discrete vs continuous vs hybrid
#
# Run interactively. Point `af.spectra` at a real library from get.af.spectra()
# if you have one -- the synthetic fallback below is smoother and lower-rank
# than real autofluorescence, so it will flatter the continuous approach.

asp.dir <- "/Users/oliverburton/Bioinformatics/AutoSpectral"
if ( requireNamespace( "devtools", quietly = TRUE ) )
  devtools::load_all( asp.dir )

set.seed( 42 )

pass <- function( label, ok, detail = "" )
  cat( sprintf( "[%s] %-52s %s\n", if ( ok ) "PASS" else "FAIL", label, detail ) )

# =============================================================================
# 0. Setup
# =============================================================================

asp <- get.autospectral.param( cytometer = "Aurora" )

ref.lib <- read.csv(
  system.file( "extdata", "Aurora_spectral_reference_library.csv",
               package = "AutoSpectral" ),
  row.names = 1, check.names = FALSE
)

panel <- c( "BUV395", "BUV496", "BUV661", "BUV737", "BUV805",
            "BV421", "BV510", "BV605", "BV650", "BV711", "BV785",
            "FITC", "PerCP", "PerCP-eFluor 710",
            "PE", "PE-Dazzle594", "PE-Cy5", "PE-Cy7",
            "APC", "AF700", "APC-Fire750" )
panel   <- intersect( panel, rownames( ref.lib ) )
spectra <- as.matrix( ref.lib[ panel, , drop = FALSE ] )
spectra <- spectra / apply( abs( spectra ), 1, max )

det.n   <- ncol( spectra )
fluor.n <- nrow( spectra )

# ---- SUBSTITUTE A REAL AF LIBRARY HERE IF YOU HAVE ONE ----------------------
# af.spectra <- get.af.spectra( ... )
#
# Synthetic fallback: n.latent underlying AF shapes, mixed to produce a library
# of n.af nodes plus node-level noise. n.latent sets the true out-of-span rank.
make.af.library <- function( n.af, n.latent, det.n, det.names,
                             node.noise = 0.02, seed = 1 ) {
  set.seed( seed )
  centres <- seq( 0.20, 0.85, length.out = n.latent ) * det.n
  latent  <- t( sapply( centres, function( cc )
    exp( -0.5 * ( ( seq_len( det.n ) - cc ) / ( det.n / 5 ) )^2 ) ) )
  mix <- matrix( runif( n.af * n.latent )^2, nrow = n.af )
  af  <- mix %*% latent + node.noise * matrix( runif( n.af * det.n ), nrow = n.af )
  af  <- pmax( af, 0 )
  af  <- af / apply( af, 1, max )
  rownames( af ) <- paste0( "AF", seq_len( n.af ) )
  colnames( af ) <- det.names
  af
}

n.af       <- 300L
af.spectra <- make.af.library( n.af, 6L, det.n, colnames( spectra ) )

sim <- sim.flow.data(
  spectra    = spectra,
  asp        = asp,
  n.cells    = 20000L,
  af.spectra = af.spectra,
  seed       = 7L
)
raw.data <- sim$raw[ , colnames( spectra ), drop = FALSE ]

cat( sprintf( "\nPanel %d x %d, AF library %d rows, %d cells\n\n",
              fluor.n, det.n, nrow( af.spectra ), nrow( raw.data ) ) )

U     <- solve.default( tcrossprod( spectra ), spectra )
Pperp <- diag( det.n ) - t( spectra ) %*% U

# =============================================================================
# 1. Out-of-span scree
# =============================================================================

cat( "--- 1. Out-of-span scree of the AF library ---\n" )

R.lib <- Pperp %*% t( af.spectra )          # D x nAF
sv    <- svd( R.lib )
cumul <- cumsum( sv$d^2 ) / sum( sv$d^2 )

# For contrast, the scree of the raw library (the wrong basis to use)
sv.raw    <- svd( t( af.spectra ) )
cumul.raw <- cumsum( sv.raw$d^2 ) / sum( sv.raw$d^2 )

cat( sprintf( "%6s %14s %14s\n", "q", "out-of-span", "raw library" ) )
for ( q in 1:10 )
  cat( sprintf( "%6d %14.5f %14.5f\n", q, cumul[ q ], cumul.raw[ q ] ) )

q.99  <- which( cumul >= 0.99  )[ 1 ]
q.999 <- which( cumul >= 0.999 )[ 1 ]
cat( sprintf( "  components for 99%%: %d    for 99.9%%: %d\n\n", q.99, q.999 ) )

# =============================================================================
# 2. Quantisation error of the discrete library
# =============================================================================

cat( "--- 2. Quantisation error currently accepted ---\n" )

# Distance from each library row's out-of-span part to the nearest OTHER row's,
# i.e. how far a cell can be forced to move when snapped to a library member.
Rn    <- sweep( R.lib, 2, pmax( sqrt( colSums( R.lib^2 ) ), 1e-12 ), "/" )
cosmat <- crossprod( Rn )
diag( cosmat ) <- -Inf
nearest.cos <- apply( cosmat, 2, max )

cat( sprintf( "  cosine to nearest neighbour (out-of-span): median %.4f, p05 %.4f\n",
              stats::median( nearest.cos ), stats::quantile( nearest.cos, 0.05 ) ) )
cat( sprintf( "  implied worst-case angular snap: median %.2f deg, p95 %.2f deg\n\n",
              stats::median( acos( pmin( nearest.cos, 1 ) ) * 180 / pi ),
              stats::quantile( acos( pmin( nearest.cos, 1 ) ) * 180 / pi, 0.95 ) ) )

# =============================================================================
# 3. unmix.af.basis() vs explicit joint OLS
# =============================================================================

cat( "--- 3. Continuous solve equivalence ---\n" )

af.basis <- get.af.basis( af.spectra, spectra, n.components = q.999 )
cat( sprintf( "  basis: %d components, %.5f of out-of-span variance\n",
              af.basis$n.components, af.basis$var.explained ) )

cont <- unmix.af.basis( raw.data, spectra, af.basis,
                        af.spectra = af.spectra, return.fitted.af = TRUE )

# reference: joint OLS on the stacked design [spectra ; t(basis)]
design <- rbind( spectra, t( af.basis$basis ) )
ref    <- unmix.ols.fast( raw.data, design )

e.f <- max( abs( cont$fluorophores - ref[ , seq_len( fluor.n ) ] ) ) /
  max( abs( ref[ , seq_len( fluor.n ) ] ) )
e.k <- max( abs( cont$k - ref[ , fluor.n + seq_len( af.basis$n.components ) ] ) ) /
  max( abs( ref[ , fluor.n + seq_len( af.basis$n.components ) ] ) )

pass( "fluorophores match joint OLS", e.f < 1e-8, sprintf( "rel.err = %.3e", e.f ) )
pass( "AF coefficients match joint OLS", e.k < 1e-8, sprintf( "rel.err = %.3e", e.k ) )

# the Gram matrix should be diagonal by construction
G <- t( af.basis$basis ) %*% Pperp %*% af.basis$basis
off <- max( abs( G - diag( diag( G ) ) ) ) / max( abs( diag( G ) ) )
pass( "B' Pperp B is diagonal", off < 1e-10, sprintf( "offdiag/diag = %.3e", off ) )
cat( "\n" )

# =============================================================================
# 4. Negative AF fraction
# =============================================================================

cat( "--- 4. Physicality of the unconstrained continuous fit ---\n" )

cat( sprintf( "  cells with meaningfully negative AF: %.2f%%\n",
              100 * mean( cont$negative ) ) )

# are they the low-AF cells?
af.rank <- rank( cont$af ) / length( cont$af )
cat( sprintf( "  median AF percentile of flagged cells: %.1f%% (of all cells: 50%%)\n\n",
              100 * stats::median( af.rank[ cont$negative ] ) ) )

# =============================================================================
# 5. Shortlist recall for the hybrid
# =============================================================================

cat( "--- 5. Hybrid shortlist recall ---\n" )

discrete.index <- assign.af.joint.cov( raw.data, spectra, af.spectra )

af.norm  <- af.spectra / sqrt( rowSums( af.spectra^2 ) )
fit.norm <- cont$fitted.af / pmax( sqrt( rowSums( cont$fitted.af^2 ) ), 1e-12 )
sim.mat  <- fit.norm %*% t( af.norm )                       # cells x nAF

sub <- sample( nrow( raw.data ), min( 4000L, nrow( raw.data ) ) )
cat( sprintf( "%6s %12s\n", "m", "recall" ) )
for ( m in c( 1L, 2L, 4L, 8L, 16L, 32L ) ) {
  hit <- vapply( sub, function( i ) {
    top <- order( sim.mat[ i, ], decreasing = TRUE )[ seq_len( m ) ]
    discrete.index[ i ] %in% top
  }, logical( 1 ) )
  cat( sprintf( "%6d %11.1f%%\n", m, 100 * mean( hit ) ) )
}
cat( "\n" )

# =============================================================================
# 6. Timing
# =============================================================================

cat( "--- 6. Timing ---\n" )

n.time <- 2e5L
y <- raw.data[ sample( nrow( raw.data ), n.time, replace = TRUE ), , drop = FALSE ]

t.discrete <- system.time( {
  idx <- assign.af.joint.cov( y, spectra, af.spectra )
  unmix.af.fwl( y, spectra, af.spectra, idx, return.fitted.af = TRUE )
} )[ "elapsed" ]

t.continuous <- system.time(
  unmix.af.basis( y, spectra, af.basis, return.fitted.af = TRUE )
)[ "elapsed" ]

cat( sprintf( "  %d cells, %d AF rows, %d basis components\n",
              n.time, nrow( af.spectra ), af.basis$n.components ) )
cat( sprintf( "  discrete   (assign + FWL solve): %.3fs\n", t.discrete ) )
cat( sprintf( "  continuous (basis solve)       : %.3fs   %.1fx\n",
              t.continuous, t.discrete / max( t.continuous, 1e-6 ) ) )
cat( "\nDone.\n" )
