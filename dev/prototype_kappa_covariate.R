# prototype_kappa_spillover_covariate.R
#
# Standalone prototype for the two-covariate noise-model regression
# ("Plan B" in CONTEXT_noise_model_from_controls.md):
#
#   Var(y_i,d) = read.var_d + mu_i,d / kappa_d
#                + x_i * p_f,d * (1 - p_f,d) / kappa_spillover_d
#
# where mu_i,d is event i's fitted mean at detector d from a single-row
# projection against its own control's reference spectrum, x_i is that
# event's fitted total abundance, and p_f,d = S_f[d] / sum_d'(S_f[d']) is
# the fraction of control f's total signal landing at detector d.
#
# Background: pooling single-stained-control residuals against a single
# covariate (mu alone) was tried in three variants and abandoned --
# CONTEXT_noise_model_from_controls.md, "What didn't work, and why".  The
# root cause: at a detector away from a control's own peak, the true
# variance is dominated by multinomial partitioning of that control's
# large, unrelated total photon flux (~ x_i * p_f,d * (1 - p_f,d)), not by
# the small local mean there. That term is what this script adds
# explicitly as a second covariate, instead of trying to mask it away.
#
# This is NOT a package function. It is meant to be run interactively,
# section by section, against 2-3 real single-stained controls whose
# reference spectra have visibly different shapes at the detector(s) you
# want to diagnose -- that shape difference across controls is what
# separates kappa_d from kappa_spillover_d; two controls with near-identical
# shapes at detector d leave that detector's split poorly conditioned.
#
# Validate against simulate.flow.data() before trusting this on real
# controls, and before any package change: is kappa_d recovered correctly
# in a simulation that includes the multinomial spillover term but not the
# single-covariate fit? Compare the two-covariate fit here against
# estimate.noise.model() run on the unstained control (Plan A) as a
# sanity check on the clean channels where both should roughly agree.
#
# Edit section 0, then run top to bottom.

library( AutoSpectral )

# ---------------------------------------------------------------------------
# 0. Setup -- EDIT THIS SECTION
# ---------------------------------------------------------------------------

asp <- get.autospectral.param( cytometer = "aurora" )   # set to your cytometer

spectra.file     <- "PLACEHOLDER_SPECTRA_FILE.csv"
spectra          <- read.spectra( spectra.file )
spectral.channel <- colnames( spectra )

# 2-3 single-stained control FCS files, named by fluorophore (must match
# rownames(spectra)).
control.file <- list(
  "PLACEHOLDER_FLUOR_1" = "PLACEHOLDER_CONTROL_1.fcs",
  "PLACEHOLDER_FLUOR_2" = "PLACEHOLDER_CONTROL_2.fcs",
  "PLACEHOLDER_FLUOR_3" = "PLACEHOLDER_CONTROL_3.fcs"
)
control.dir <- "PLACEHOLDER_CONTROL_DIR"

# Used only to set each control's own raw-channel positivity threshold (its
# peak detector), the same role raw.thresholds plays in get.fluor.variants().
unstained.file <- "PLACEHOLDER_UNSTAINED_CONTROL.fcs"

peak.channel <- c(                          # peak detector per control fluor
  "PLACEHOLDER_FLUOR_1" = "PLACEHOLDER_CHANNEL_1",
  "PLACEHOLDER_FLUOR_2" = "PLACEHOLDER_CHANNEL_2",
  "PLACEHOLDER_FLUOR_3" = "PLACEHOLDER_CHANNEL_3"
)

threshold.quantile     <- 0.995   # positivity cut, on the unstained control
max.events.per.control <- 20000L
n.bins.mu  <- 6L                  # quantile bins on the mu axis
n.bins.z   <- 4L                  # quantile bins on the z (partition) axis, within each mu bin
min.cell.n <- 30L                 # minimum events per 2D bin
seed       <- 42L

# ---------------------------------------------------------------------------
# 1. Read controls, single-row fit, compute the two covariates
# ---------------------------------------------------------------------------

unstained <- readFCS( unstained.file, columns = spectral.channel )
raw.thresholds <- apply( unstained, 2, function( col )
  stats::quantile( col, threshold.quantile ) )

fluors <- names( control.file )
stopifnot( all( fluors %in% rownames( spectra ) ) )

per.control <- lapply( fluors, function( fl ) {

  message( "Reading ", fl, " ..." )
  raw <- readFCS( file.path( control.dir, control.file[[ fl ]] ), columns = spectral.channel )

  # drop saturated events, matching get.fluor.variants()'s own filter
  keep <- rowSums( raw >= asp$expr.data.max ) == 0
  raw  <- raw[ keep, , drop = FALSE ]

  pos.idx <- which( raw[ , peak.channel[ fl ] ] > raw.thresholds[ peak.channel[ fl ] ] )
  raw <- raw[ pos.idx, , drop = FALSE ]

  if ( nrow( raw ) > max.events.per.control ) {
    set.seed( seed )
    raw <- raw[ sample( nrow( raw ), max.events.per.control ), , drop = FALSE ]
  }

  ref   <- spectra[ fl, , drop = FALSE ]
  x.hat <- as.numeric( unmix.ols( raw, ref ) )       # event-level fitted abundance
  mu    <- outer( x.hat, as.numeric( ref ) )         # events x detectors fitted mean
  resid <- raw - mu

  p.f <- as.numeric( ref ) / sum( ref )              # partition fraction per detector
  z   <- outer( x.hat, p.f * ( 1 - p.f ) )            # events x detectors partition term

  list(
    fluor = fl,
    n     = nrow( raw ),
    mu    = mu,
    z     = z,
    resid = resid
  )
} )
names( per.control ) <- fluors

message( sprintf( "Pooled events: %s",
                  paste( sprintf( "%s = %d", fluors, vapply( per.control, function( x ) x$n, integer( 1 ) ) ),
                         collapse = ", " ) ) )

pool.mu    <- do.call( rbind, lapply( per.control, function( x ) x$mu ) )
pool.z     <- do.call( rbind, lapply( per.control, function( x ) x$z  ) )
pool.resid <- do.call( rbind, lapply( per.control, function( x ) x$resid ) )
pool.fluor <- unlist( lapply( per.control, function( x ) rep( x$fluor, x$n ) ) )

colnames( pool.mu ) <- colnames( pool.z ) <- colnames( pool.resid ) <- spectral.channel

# ---------------------------------------------------------------------------
# 2. Per-detector two-covariate regression (2D quantile binning + IRLS)
# ---------------------------------------------------------------------------
# Mirrors .fit.noise.regression()'s binning/IRLS scheme (estimate_noise_model.R),
# extended to a second covariate. Binning (rather than a raw event-level
# squared-residual fit) keeps this comparable to the existing single-covariate
# diagnostics and keeps the chi-square-heavy tail of individual r^2 values
# from dominating the fit.

fit.two.covariate <- function( mu.d, z.d, resid.d, n.bins.mu, n.bins.z, min.cell.n ) {

  ok <- is.finite( mu.d ) & is.finite( z.d ) & is.finite( resid.d )
  mu.d <- mu.d[ ok ]; z.d <- z.d[ ok ]; resid.d <- resid.d[ ok ]
  if ( length( mu.d ) < min.cell.n * n.bins.mu * n.bins.z )
    return( NULL )

  mu.breaks <- unique( stats::quantile( mu.d, probs = seq( 0, 1, length.out = n.bins.mu + 1L ) ) )
  if ( length( mu.breaks ) < 3L ) return( NULL )
  mu.bin <- cut( mu.d, breaks = mu.breaks, include.lowest = TRUE, labels = FALSE )

  cell.rows <- list()
  for ( b in seq_len( length( mu.breaks ) - 1L ) ) {

    idx <- which( mu.bin == b )
    if ( length( idx ) < min.cell.n * 2L ) next

    z.breaks <- unique( stats::quantile( z.d[ idx ], probs = seq( 0, 1, length.out = n.bins.z + 1L ) ) )
    if ( length( z.breaks ) < 3L ) {
      if ( length( idx ) >= min.cell.n )
        cell.rows[[ length( cell.rows ) + 1L ]] <- data.frame(
          mean.mu = mean( mu.d[ idx ] ), mean.z = mean( z.d[ idx ] ),
          variance = stats::mad( resid.d[ idx ] )^2, n = length( idx ) )
      next
    }

    z.bin <- cut( z.d[ idx ], breaks = z.breaks, include.lowest = TRUE, labels = FALSE )
    for ( zb in seq_len( length( z.breaks ) - 1L ) ) {
      sub <- idx[ z.bin == zb ]
      if ( length( sub ) < min.cell.n ) next
      cell.rows[[ length( cell.rows ) + 1L ]] <- data.frame(
        mean.mu = mean( mu.d[ sub ] ), mean.z = mean( z.d[ sub ] ),
        variance = stats::mad( resid.d[ sub ] )^2, n = length( sub ) )
    }
  }

  if ( length( cell.rows ) < 6L ) return( NULL )   # need enough cells for a 3-parameter fit
  cell.df <- do.call( rbind, cell.rows )
  cell.df <- cell.df[ is.finite( cell.df$variance ) & cell.df$variance > 0, ]
  if ( nrow( cell.df ) < 6L ) return( NULL )

  w   <- cell.df$n / pmax( cell.df$variance, 1e-8 )^2
  fit <- stats::lm( variance ~ mean.mu + mean.z, data = cell.df, weights = w )
  for ( pass in seq_len( 2L ) ) {
    v.fit <- pmax( stats::fitted( fit ), 1e-8 )
    w     <- cell.df$n / v.fit^2
    fit   <- stats::lm( variance ~ mean.mu + mean.z, data = cell.df, weights = w )
  }

  cf <- stats::coef( fit )
  list(
    read.var        = max( cf[ "(Intercept)" ], 1e-8 ),
    kappa           = if ( is.finite( cf[ "mean.mu" ] ) && cf[ "mean.mu" ] > 0 ) 1 / cf[ "mean.mu" ] else NA_real_,
    kappa.spillover = if ( is.finite( cf[ "mean.z" ] )  && cf[ "mean.z" ]  > 0 ) 1 / cf[ "mean.z" ]  else NA_real_,
    r.squared       = summary( fit )$r.squared,
    cell.table      = cell.df
  )
}

results <- lapply( spectral.channel, function( d )
  fit.two.covariate( pool.mu[ , d ], pool.z[ , d ], pool.resid[ , d ],
                     n.bins.mu, n.bins.z, min.cell.n ) )
names( results ) <- spectral.channel

fit.ok <- !vapply( results, is.null, logical( 1 ) )
summary.df <- data.frame(
  detector        = spectral.channel[ fit.ok ],
  read.var        = vapply( results[ fit.ok ], function( r ) r$read.var, numeric( 1 ) ),
  kappa           = vapply( results[ fit.ok ], function( r ) r$kappa, numeric( 1 ) ),
  kappa.spillover = vapply( results[ fit.ok ], function( r ) r$kappa.spillover, numeric( 1 ) ),
  r.squared       = vapply( results[ fit.ok ], function( r ) r$r.squared, numeric( 1 ) )
)
print( summary.df )

cat( sprintf( "\nFit succeeded at %d / %d detectors.\n", sum( fit.ok ), length( spectral.channel ) ) )

# ---------------------------------------------------------------------------
# 3. Diagnostic: does adding z remove the dip-then-rise shape?
# ---------------------------------------------------------------------------
# Repeats the single-covariate fit (mu only, no z, no masking band -- the
# thing Attempts A-C showed was non-monotonic) on the same pooled data, for
# side-by-side comparison at a detector of interest.

fit.single.covariate <- function( mu.d, resid.d, n.bins, min.cell.n ) {
  ok <- is.finite( mu.d ) & is.finite( resid.d )
  mu.d <- mu.d[ ok ]; resid.d <- resid.d[ ok ]
  breaks <- unique( stats::quantile( mu.d, probs = seq( 0, 1, length.out = n.bins + 1L ) ) )
  bin    <- cut( mu.d, breaks = breaks, include.lowest = TRUE, labels = FALSE )
  bin.n  <- tabulate( bin, nbins = length( breaks ) - 1L )
  good   <- which( bin.n >= min.cell.n )
  data.frame(
    mean.mu  = vapply( good, function( b ) mean( mu.d[ bin == b ] ), numeric( 1 ) ),
    variance = vapply( good, function( b ) stats::mad( resid.d[ bin == b ] )^2, numeric( 1 ) ),
    n        = bin.n[ good ]
  )
}

detector.to.check <- "PLACEHOLDER_DETECTOR"   # e.g. "UV1-A" -- pick one flagged
# non-monotonic in the original
# single-covariate fit.table

single.cov <- fit.single.covariate(
  pool.mu[ , detector.to.check ], pool.resid[ , detector.to.check ],
  n.bins = n.bins.mu * n.bins.z, min.cell.n = min.cell.n )

cat( "\nSingle-covariate (mu only) bins at", detector.to.check, ":\n" )
print( single.cov )

cat( "\nTwo-covariate fit at", detector.to.check, ":\n" )
print( results[[ detector.to.check ]]$cell.table )

if ( interactive() ) {
  op <- par( mfrow = c( 1, 2 ) )
  plot( single.cov$mean.mu, single.cov$variance, pch = 16,
        xlab = "fitted mean (mu)", ylab = "residual variance",
        main = paste( detector.to.check, "-- mu only" ) )
  ct <- results[[ detector.to.check ]]$cell.table
  plot( ct$mean.z, ct$variance, pch = 16, cex = 0.3 + 2 * ( ct$mean.mu / max( ct$mean.mu ) ),
        xlab = "partition term (z)", ylab = "residual variance",
        main = paste( detector.to.check, "-- vs z, point size ~ mu" ) )
  par( op )
}

# ---------------------------------------------------------------------------
# PASS/fail heuristic for this prototype
# ---------------------------------------------------------------------------
# The two-covariate fit should show materially higher R^2 than the
# single-covariate fit at a detector flagged non-monotonic in the original
# attempt, and its per-cell variance, at roughly constant mean.z, should no
# longer reverse sign as mean.mu rises (check cell.table directly: group by
# similar mean.z and look for monotone increase in variance with mean.mu).
# If that holds up across 2-3 different control combinations, kappa/kappa.spillover
# from step 2 are the candidates for a package version of this fit; validate
# against simulate.flow.data() before wiring it into get.spectral.variants().
