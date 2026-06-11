# simulate_flow_data.R

#' @title Simulate Spectral Flow Cytometry Data
#'
#' @description
#' Generates synthetic spectral flow cytometry data for development and
#' validation of AutoSpectral pipelines. Cells are assigned random fluorophore
#' expression levels drawn from three discrete expression layers (high, medium,
#' low), mixed into detector space via provided spectra (and optionally spectral
#' variants), and corrupted by configurable noise sources: binomial spillover
#' sampling, Poisson shot noise, and detector readout noise. Autofluorescence
#' variation is supported via a provided `af.spectra` matrix.
#'
#' Ground-truth abundances, per-cell variant assignments, per-cell AF spectra,
#' and synthetic scatter are returned alongside the raw detector-space matrix
#' and a simple OLS unmixed result (no AF) for immediate inspection.
#'
#' @param spectra Numeric matrix of fluorophore spectral signatures, L-infinity
#'   normalised, fluorophores in rows and detectors in columns. Required.
#'   Defines both the panel and the detector layout.
#' @param asp The AutoSpectral parameter list from `get.autospectral.param`.
#'   Used to resolve `expr.data.max` and cytometer-specific detector noise
#'   presets.
#' @param n.cells Integer. Number of synthetic cells to generate. Default
#'   `10000`.
#' @param complexity Numeric in (0, 1]. Fraction of fluorophores that are "on"
#'   (non-zero) per cell, drawn without replacement. Default `1/3`.
#' @param layer.centroids Numeric vector of length 3, giving the high, medium,
#'   and low expression centroids as fractions of `expr.data.max`. Default
#'   `c(0.5, 0.05, 0.005)`.
#' @param layer.cv Numeric scalar. Coefficient of variation (on the log scale)
#'   for expression within each layer. Default `0.2`.
#' @param af.spectra Numeric matrix of autofluorescence spectra, L-infinity
#'   normalised, AF components in rows and detectors in columns. As returned by
#'   `get.af.spectra`. If `NULL` (default), no AF is added.
#' @param af.scale.range Numeric vector of length 2. Range of the uniform
#'   distribution from which per-cell AF abundance is drawn, as a fraction of
#'   `expr.data.max`. Default `c(0.001, 0.05)`.
#' @param af.scatter.slope Numeric scalar. Controls how much SSC drives AF
#'   abundance. AF abundance is multiplied by
#'   `1 + af.scatter.slope * (ssc - median(ssc)) / mad(ssc)`. Default `0.1`.
#' @param variants List of spectral variant matrices as returned by
#'   `get.spectral.variants` or `get.fluor.variants`. Each element is named by
#'   fluorophore and contains a matrix of variant spectra (variants in rows,
#'   detectors in columns). If `NULL` (default), the base `spectra` rows are
#'   used for all cells.
#' @param scatter.data Optional numeric matrix of real scatter data (cells in
#'   rows, at least two columns for FSC and SSC). When provided, rows are
#'   sampled with replacement to generate synthetic scatter, bypassing the
#'   log-normal model. Default `NULL`.
#' @param fsc.mean.log Numeric. Mean of the log-normal FSC distribution (on the
#'   natural-log scale). Default `log(200000)`.
#' @param fsc.sd.log Numeric. SD of the log-normal FSC distribution. Default
#'   `0.4`.
#' @param ssc.mean.log Numeric. Mean of the log-normal SSC distribution.
#'   Default `log(80000)`.
#' @param ssc.sd.log Numeric. SD of the log-normal SSC distribution. Default
#'   `0.5`.
#' @param fsc.ssc.cor Numeric in (-1, 1). Pearson correlation between log(FSC)
#'   and log(SSC) in the synthetic scatter population. Default `0.6`.
#' @param shot.noise Logical. Whether to apply Poisson shot noise. Default
#'   `TRUE`.
#' @param counts.per.unit Numeric scalar for `shot.noise`. Default `0.5`. Works
#'   inversely; higher numbers will decrease shot noise-driven spread.
#' @param spillover.noise Logical. Whether to apply binomial spillover sampling.
#'   Default `TRUE`.
#' @param detector.noise Logical. Whether to apply detector readout noise
#'   (Gaussian, cytometer-specific). Default `TRUE`.
#' @param af.variation Logical. Whether to include autofluorescence variation
#'   (requires `af.spectra`). Default `TRUE`.
#' @param spectral.variation Logical. Whether to sample per-cell spectral
#'   variants (requires `variants`). Default `TRUE`.
#' @param seed Integer random seed for reproducibility. Default `42`.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{`raw`}{Numeric matrix (cells x detectors) of simulated detector-space
#'     data, including all enabled noise sources.}
#'   \item{`unmixed.no.af`}{Numeric matrix (cells x fluorophores) of OLS-unmixed
#'     data without AF, for immediate inspection.}
#'   \item{`truth`}{Named list:
#'     \describe{
#'       \item{`abundances`}{Matrix (cells x fluorophores+AF), true abundance per
#'         cell. AF column appended as the last column when `af.spectra` is
#'         provided.}
#'       \item{`expression.layer`}{Character matrix (cells x fluorophores),
#'         `"high"`, `"mid"`, `"low"`, or `NA` (off).}
#'       \item{`variant.index`}{Integer matrix (cells x fluorophores), index of
#'         the variant spectrum used per cell per fluorophore. `NA` when off or
#'         no variants provided.}
#'       \item{`af.row`}{Integer vector (length cells), which row of `af.spectra`
#'         was assigned to each cell. `NA` when no AF.}
#'       \item{`af.abundance`}{Numeric vector (length cells), true AF abundance
#'         per cell in instrument units. `0` when no AF.}
#'       \item{`scatter`}{Numeric matrix (cells x 2, columns FSC and SSC),
#'         synthetic or resampled scatter values.}
#'     }
#'   }
#'   \item{`params`}{Named list of all simulation settings used, for
#'     reproducibility.}
#' }
#'
#' @export

sim.flow.data <- function(
    spectra,
    asp,
    n.cells              = 10000L,
    complexity           = 1/3,
    layer.centroids      = c( 0.5, 0.05, 0.005 ),
    layer.cv             = 0.2,
    af.spectra           = NULL,
    af.scale.range       = c( 0.001, 0.05 ),
    af.scatter.slope     = 0.1,
    variants             = NULL,
    scatter.data         = NULL,
    fsc.mean.log         = log( 200000 ),
    fsc.sd.log           = 0.4,
    ssc.mean.log         = log( 80000 ),
    ssc.sd.log           = 0.5,
    fsc.ssc.cor          = 0.6,
    shot.noise           = TRUE,
    counts.per.unit      = 0.5,
    spillover.noise      = TRUE,
    detector.noise       = TRUE,
    af.variation         = TRUE,
    spectral.variation   = TRUE,
    seed                 = 42L
) {

  set.seed( seed )

  # ---------------------------------------------------------------------------
  # Validation
  # ---------------------------------------------------------------------------

  stopifnot(
    is.matrix( spectra ) || is.data.frame( spectra ),
    is.numeric( n.cells ), length( n.cells ) == 1L, n.cells >= 1L,
    is.numeric( complexity ), length( complexity ) == 1L,
    complexity > 0, complexity <= 1,
    is.numeric( layer.centroids ), length( layer.centroids ) == 3L,
    all( layer.centroids > 0 ), all( layer.centroids <= 1 ),
    is.numeric( layer.cv ), length( layer.cv ) == 1L, layer.cv > 0
  )

  spectra     <- as.matrix( spectra )
  n.fluors    <- nrow( spectra )
  n.detectors <- ncol( spectra )
  fluor.names <- rownames( spectra )
  det.names   <- colnames( spectra )

  if ( is.null( fluor.names ) )
    fluor.names <- paste0( "F", seq_len( n.fluors ) )

  expr.max <- asp$expr.data.max
  if ( is.null( expr.max ) )
    stop( "asp$expr.data.max not found. Please use a valid asp from get.autospectral.param()." )

  # ---------------------------------------------------------------------------
  # Detector noise preset (bundled per cytometer)
  # ---------------------------------------------------------------------------
  # counts.per.unit : photon-equivalent counts per instrument unit (for shot
  #   noise and binomial sampling)
  # readout.sd      : Gaussian readout noise SD in instrument units
  # dark.offset     : mean dark current offset subtracted by instrument
  #   (centres near-zero signals around zero, producing negatives)

  detector.preset <- .sim.detector.preset( asp$cytometer )
  # allow user override of counts.per.unit for noise tuning
  if ( !is.null( counts.per.unit ) )
    detector.preset$counts.per.unit <- counts.per.unit

  # ---------------------------------------------------------------------------
  # Scatter
  # ---------------------------------------------------------------------------

  scatter <- .sim.scatter(
    n.cells      = n.cells,
    scatter.data = scatter.data,
    fsc.mean.log = fsc.mean.log,
    fsc.sd.log   = fsc.sd.log,
    ssc.mean.log = ssc.mean.log,
    ssc.sd.log   = ssc.sd.log,
    fsc.ssc.cor  = fsc.ssc.cor
  )

  # ---------------------------------------------------------------------------
  # Expression layer centroids in instrument units
  # ---------------------------------------------------------------------------

  layer.units <- layer.centroids * expr.max

  # log-normal parameters: mu and sigma on the natural-log scale
  # CV = sqrt( exp(sigma^2) - 1 )  =>  sigma = sqrt( log( CV^2 + 1 ) )
  layer.sigma <- sqrt( log( layer.cv^2 + 1 ) )
  layer.mu    <- log( layer.units ) - 0.5 * layer.sigma^2

  names( layer.mu )    <- c( "high", "mid", "low" )

  # ---------------------------------------------------------------------------
  # Per-cell fluorophore assignments: which are on, which layer, abundance
  # ---------------------------------------------------------------------------

  n.on.per.cell <- max( 1L, round( complexity * n.fluors ) )

  # abundance matrix: cells x fluorophores
  abundance     <- matrix( 0, nrow = n.cells, ncol = n.fluors,
                           dimnames = list( NULL, fluor.names ) )
  expr.layer    <- matrix( NA_character_, nrow = n.cells, ncol = n.fluors,
                           dimnames = list( NULL, fluor.names ) )
  variant.index <- matrix( NA_integer_,  nrow = n.cells, ncol = n.fluors,
                           dimnames = list( NULL, fluor.names ) )

  layer.names <- c( "high", "mid", "low" )

  for ( i in seq_len( n.cells ) ) {

    on.fluors <- sample( seq_len( n.fluors ), n.on.per.cell, replace = FALSE )

    for ( f in on.fluors ) {

      lyr <- sample( layer.names, 1L )
      expr.layer[ i, f ] <- lyr

      abundance[ i, f ] <- exp(
        stats::rnorm( 1L, mean = layer.mu[ lyr ], sd = layer.sigma )
      )

    }
  }

  # ---------------------------------------------------------------------------
  # Per-cell spectral mixing: build detector-space signal
  # ---------------------------------------------------------------------------
  # For each cell, mix fluorophore contributions using per-cell variant spectra
  # (or base spectra when variants are unavailable / disabled).

  # signal.clean: noiseless detector-space signal (cells x detectors)
  signal.clean <- matrix( 0, nrow = n.cells, ncol = n.detectors,
                          dimnames = list( NULL, det.names ) )

  use.variants <- spectral.variation && !is.null( variants )

  for ( f in seq_len( n.fluors ) ) {

    fname   <- fluor.names[ f ]
    f.abund <- abundance[ , f ]

    # cells where this fluorophore is on
    on.cells <- which( f.abund > 0 )
    if ( length( on.cells ) == 0L ) next

    if ( use.variants && fname %in% names( variants ) ) {

      var.mat <- as.matrix( variants[[ fname ]] )
      n.vars  <- nrow( var.mat )

      v.idx   <- sample( seq_len( n.vars ), length( on.cells ), replace = TRUE )
      variant.index[ on.cells, f ] <- v.idx

      for ( k in seq_along( on.cells ) ) {
        ci <- on.cells[ k ]
        signal.clean[ ci, ] <- signal.clean[ ci, ] +
          f.abund[ ci ] * var.mat[ v.idx[ k ], ]
      }

    } else {

      # base spectrum for all on-cells
      base.spec <- spectra[ f, ]
      signal.clean[ on.cells, ] <- signal.clean[ on.cells, ] +
        outer( f.abund[ on.cells ], base.spec )

    }
  }

  # ---------------------------------------------------------------------------
  # Autofluorescence
  # ---------------------------------------------------------------------------

  af.row       <- rep( NA_integer_, n.cells )
  af.abundance <- rep( 0,           n.cells )

  use.af <- af.variation && !is.null( af.spectra )

  if ( use.af ) {

    af.spectra <- as.matrix( af.spectra )
    n.af       <- nrow( af.spectra )

    # SSC z-score for AF-scatter coupling
    ssc.vals <- scatter[ , 2L ]
    ssc.med  <- stats::median( ssc.vals )
    ssc.mad  <- max( stats::mad( ssc.vals ), 1 )   # guard against zero MAD
    ssc.z    <- ( ssc.vals - ssc.med ) / ssc.mad

    # base AF abundance: uniform on af.scale.range * expr.max
    af.base <- stats::runif(
      n.cells,
      min = af.scale.range[ 1L ] * expr.max,
      max = af.scale.range[ 2L ] * expr.max
    )

    # scatter modulation: multiplicative, floored at 0
    af.abundance <- pmax( 0, af.base * ( 1 + af.scatter.slope * ssc.z ) )

    # assign one AF spectrum per cell (uniform draw)
    af.row <- sample( seq_len( n.af ), n.cells, replace = TRUE )

    for ( i in seq_len( n.cells ) ) {
      signal.clean[ i, ] <- signal.clean[ i, ] +
        af.abundance[ i ] * af.spectra[ af.row[ i ], ]
    }
  }

  # ---------------------------------------------------------------------------
  # Noise pipeline
  # ---------------------------------------------------------------------------

  signal <- signal.clean   # will be mutated through noise stages

  cpu <- detector.preset$counts.per.unit

  # Stage 1: Binomial spillover sampling
  # Each detector channel's contribution from each fluorophore is treated as a
  # binomial draw from the total photon count for that fluorophore.
  # Here we approximate at the per-detector level: we convert the mixed signal
  # to counts, binomially sample (n = counts, p = 1 with variance correction),
  # and convert back. The variance introduced is Binomial(n, p_ij) where p_ij
  # is the spectral coefficient -- this adds variance proportional to signal.

  if ( spillover.noise ) {

    # Binomial spillover sampling operates per-fluorophore before mixing.
    # For fluorophore f, each cell emits N_f ~ Poisson( abundance_f * cpu )
    # total photons. Those photons are multinomially distributed across
    # detectors with probabilities given by the spectral row. The variance
    # added per detector j is N_f * p_fj * (1 - p_fj).
    # We use the Normal approximation: add N( 0, sqrt( N_f * p_fj * (1-p_fj) ) )
    # to the already-mixed signal (equivalent to adding the residual variance
    # from the multinomial after the mean has already been placed by mixing).

    spillover.noise.mat <- matrix( 0, nrow = n.cells, ncol = n.detectors )

    for ( f in seq_len( n.fluors ) ) {

      f.abund    <- abundance[ , f ]
      on.cells   <- which( f.abund > 0 )
      if ( length( on.cells ) == 0L ) next

      # per-cell photon counts for this fluorophore
      n.photons <- f.abund[ on.cells ] * detector.preset$counts.per.unit

      # spectral probabilities for this fluorophore (use base spectrum;
      # variant spectra are L-inf normalised so we normalise to sum = 1
      # here for valid multinomial probabilities)
      spec.row <- spectra[ f, ]
      spec.row <- pmax( spec.row, 0 )
      spec.sum <- sum( spec.row )
      if ( spec.sum <= 0 ) next
      p.vec <- spec.row / spec.sum

      # Normal approximation to Binomial variance per detector
      for ( j in seq_len( n.detectors ) ) {
        p_j   <- p.vec[ j ]
        binom.sd <- sqrt( n.photons * p_j * ( 1 - p_j ) )
        spillover.noise.mat[ on.cells, j ] <- spillover.noise.mat[ on.cells, j ] +
          stats::rnorm( length( on.cells ), mean = 0, sd = binom.sd )
      }
    }

    # convert noise from photon counts back to instrument units and add
    signal <- signal + spillover.noise.mat / detector.preset$counts.per.unit

  }

  # Stage 2: Poisson shot noise
  # Applied to the integrated photon count per detector per cell.

  if ( shot.noise ) {

    signal.counts <- pmax( 0, signal * cpu )

    # rpois requires integer lambda; use rounded counts as lambda
    signal.counts.noisy <- matrix(
      stats::rpois( n.cells * n.detectors, lambda = as.vector( signal.counts ) ),
      nrow = n.cells,
      ncol = n.detectors
    )

    signal <- signal.counts.noisy / cpu

  }

  # Stage 3: Detector readout noise + dark current subtraction
  # Gaussian readout noise is added, then the dark offset is subtracted,
  # replicating the instrument's background correction and naturally
  # producing negative values for near-zero signals.

  if ( detector.noise ) {

    readout.noise <- matrix(
      stats::rnorm( n.cells * n.detectors,
                    mean = 0,
                    sd   = detector.preset$readout.sd ),
      nrow = n.cells,
      ncol = n.detectors
    )

    signal <- signal + readout.noise - detector.preset$dark.offset

  }

  colnames( signal ) <- det.names

  # ---------------------------------------------------------------------------
  # OLS unmix (no AF) for immediate inspection
  # ---------------------------------------------------------------------------

  unmixed.no.af <- unmix.ols.fast( signal, spectra )

  # ---------------------------------------------------------------------------
  # Ground truth: abundances with AF appended
  # ---------------------------------------------------------------------------

  if ( use.af ) {
    truth.abundances <- cbind( abundance, AF = af.abundance )
  } else {
    truth.abundances <- abundance
  }

  # ---------------------------------------------------------------------------
  # Return
  # ---------------------------------------------------------------------------

  list(
    raw           = signal,
    unmixed.no.af = unmixed.no.af,
    truth = list(
      abundances       = truth.abundances,
      expression.layer = expr.layer,
      variant.index    = variant.index,
      af.row           = af.row,
      af.abundance     = af.abundance,
      scatter          = scatter
    ),
    params = list(
      n.cells            = n.cells,
      complexity         = complexity,
      layer.centroids    = layer.centroids,
      layer.cv           = layer.cv,
      af.scale.range     = af.scale.range,
      af.scatter.slope   = af.scatter.slope,
      shot.noise         = shot.noise,
      spillover.noise    = spillover.noise,
      detector.noise     = detector.noise,
      af.variation       = af.variation && !is.null( af.spectra ),
      spectral.variation = spectral.variation && !is.null( variants ),
      cytometer          = asp$cytometer,
      expr.data.max      = expr.max,
      counts.per.unit    = detector.preset$counts.per.unit,
      detector.preset    = detector.preset,
      seed               = seed
    )
  )
}


# ---------------------------------------------------------------------------
# Private helper: scatter generation
# ---------------------------------------------------------------------------

.sim.scatter <- function(
    n.cells,
    scatter.data,
    fsc.mean.log,
    fsc.sd.log,
    ssc.mean.log,
    ssc.sd.log,
    fsc.ssc.cor
) {

  if ( !is.null( scatter.data ) ) {

    scatter.data <- as.matrix( scatter.data )
    if ( ncol( scatter.data ) < 2L )
      stop( "scatter.data must have at least two columns (FSC and SSC)." )

    idx <- sample( seq_len( nrow( scatter.data ) ), n.cells, replace = TRUE )
    out <- scatter.data[ idx, 1:2, drop = FALSE ]
    colnames( out ) <- c( "FSC", "SSC" )
    return( out )

  }

  # Correlated bivariate log-normal via Cholesky
  cov.mat <- matrix(
    c( fsc.sd.log^2,
       fsc.ssc.cor * fsc.sd.log * ssc.sd.log,
       fsc.ssc.cor * fsc.sd.log * ssc.sd.log,
       ssc.sd.log^2 ),
    nrow = 2
  )

  chol.mat <- chol( cov.mat )
  z        <- matrix( stats::rnorm( n.cells * 2L ), nrow = n.cells, ncol = 2L )
  log.vals <- z %*% chol.mat

  fsc <- exp( log.vals[ , 1L ] + fsc.mean.log )
  ssc <- exp( log.vals[ , 2L ] + ssc.mean.log )

  matrix( c( fsc, ssc ), nrow = n.cells, ncol = 2L,
          dimnames = list( NULL, c( "FSC", "SSC" ) ) )
}


# ---------------------------------------------------------------------------
# Private helper: cytometer detector noise presets
# ---------------------------------------------------------------------------
#
# counts.per.unit : photon-equivalent counts per instrument intensity unit.
#   Used to convert intensity to photon counts for shot noise and binomial
#   spillover sampling. Higher = finer quantisation relative to signal.
#
# readout.sd      : SD of additive Gaussian readout noise in instrument units.
#   Represents amplification variance (PMT excess noise factor, APD readout,
#   SiPM dark count contribution) after the photon-detection stage.
#
# dark.offset     : mean dark current subtracted by the instrument's
#   background correction. Near-zero signals centre around this value, so
#   subtracting it produces the negative values seen in real data.
#   Expressed in instrument units.

.sim.detector.preset <- function( cytometer ) {

  presets <- list(

    # Cytek Aurora / Northern Lights (SiPM)
    # SiPMs have low excess noise but non-trivial dark count rates.
    # expr.data.max ~ 4.2M; moderate counts.per.unit.
    "Aurora" = list(
      counts.per.unit = 5e-3,
      readout.sd      = 80,
      dark.offset     = 50
    ),

    # Sony ID7000 (PMT)
    # PMT dynode cascade variance; expr.data.max ~ 1M.
    "ID7000" = list(
      counts.per.unit = 8e-3,
      readout.sd      = 60,
      dark.offset     = 40
    ),

    # BD FACSDiscover S8 (PMT + APD mixed array)
    # Larger dynamic range (expr.data.max ~ 24M); higher readout noise
    # consistent with observed deeper negatives.
    "FACSDiscover S8" = list(
      counts.per.unit = 2e-3,
      readout.sd      = 200,
      dark.offset     = 120
    ),

    # BD FACSDiscover A8 (PMT + APD mixed array)
    # Same detector architecture as S8.
    "FACSDiscover A8" = list(
      counts.per.unit = 2e-3,
      readout.sd      = 200,
      dark.offset     = 120
    ),

    # Beckman Coulter CytoFLEX Opteon (APD)
    # APD excess noise lower than PMT; expr.data.max ~ 1M.
    "Opteon" = list(
      counts.per.unit = 8e-3,
      readout.sd      = 50,
      dark.offset     = 30
    ),

    # Luminex Mosaic (PMT)
    "Mosaic" = list(
      counts.per.unit = 6e-3,
      readout.sd      = 70,
      dark.offset     = 45
    ),

    # Becton Dickinson Xenith (PMT)
    "Xenith" = list(
      counts.per.unit = 6e-3,
      readout.sd      = 70,
      dark.offset     = 45
    ),

    # BD FACSymphony A5 SE (PMT)
    "Symphony" = list(
      counts.per.unit = 7e-3,
      readout.sd      = 65,
      dark.offset     = 40
    )
  )

  preset <- presets[[ cytometer ]]

  if ( is.null( preset ) ) {
    warning(
      "No detector preset found for cytometer '", cytometer,
      "'. Using generic PMT defaults."
    )
    preset <- list(
      counts.per.unit = 6e-3,
      readout.sd      = 80,
      dark.offset     = 50
    )
  }

  preset
}
