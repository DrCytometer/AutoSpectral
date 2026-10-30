# residual_spread_model.R
#
# Stand-alone replica of the Residual Model spillover-spread estimate made by
# get.fluor.variants() / get.spectral.variants() (Cai et al. 2026), for ONE
# single-stained control. It follows the package code step by step:
#
#  1. remove saturated events; keep events above the 99.5th percentile of the
#     unstained control in the raw peak channel (top 2 * n.cells by peak);
#  2. subtract the scatter-matched universal negative (k nearest unstained
#     neighbours in FSC/SSC);
#  3. low-rank fit with unmix.ols(): AF principal components (SVD of the
#     unstained control) + the fluorophore's own spectrum (fluorophore alone
#     if it is nearly collinear with an AF component), giving the control's
#     own abundance b.h and the residual of that fit;
#  4-5. project the residual through the FULL panel and fit the Residual
#     Model line Var(unmixed_k | b.h) = intercept_k + slope_k * b.h for every
#     channel k with the package's own estimator (internal
#     .residual.model.line(), the code get.fluor.variants() uses; default
#     "binned.mad": Huber line of binned squared MADs on the bins' mean b.h).
#     The slope is SS(fluorophore, k) in the units used throughout AutoSpectral.
#
# Differences from the package: (a) it works on events already gated by the
# caller (the package reads the control ungated), and (b) it does not need a
# control file or the whole panel of controls. Steps 1-3 mirror
# get.fluor.variants(); the estimator in steps 4-5 is the package's.

#' @param stained.data,unstained.data Matrices of GATED raw events containing
#'   the spectral channels (`colnames(spectra)`) and `scatter.channels`.
#' @param spectra Full-panel spectral matrix (fluorophores x detectors) used
#'   for unmixing (the correct, cell-derived matrix).
#' @param fluorophore Row of `spectra` for the stained control.
#' @param peak.channel Peak detector of the fluorophore.
#' @param asp AutoSpectral parameter list.
#' @param scatter.channels FSC/SSC channel names (default
#'   `asp$default.scatter.parameter`).
#' @param spread.estimator,spread.n.bins,spread.min.bin.events,huber.k,huber.max.iter
#'   As in AutoSpectral get.fluor.variants() (default "binned.mad").
#' @return list(spread = named numeric slope per unmixed channel (NA for the
#'   fluorophore itself; negative slopes set to 0, as in
#'   get.spectral.variants()), intercept, n.events, af.collinear, af.cosine,
#'   estimator).
calculate.residual.spread <- function(
    stained.data,
    unstained.data,
    spectra,
    fluorophore,
    peak.channel,
    asp,
    scatter.channels = asp$default.scatter.parameter,
    raw.threshold.quantile = 0.995,
    n.cells = 10000L,
    k.neighbors = 3L,
    n.af.pcs = 4L,
    af.collinear.threshold = 0.95,
    spread.estimator = c( "binned.mad", "huber.squared" ),
    spread.n.bins = 20L,
    spread.min.bin.events = 20L,
    huber.k = 1.345,
    huber.max.iter = 100L,
    min.events = 50L
) {

  spread.estimator <- match.arg( spread.estimator )
  spec.ch <- colnames( spectra )

  # saturation filter
  drop.saturated <- function( m ) m[ rowSums( m[ , spec.ch, drop = FALSE ] >= asp$expr.data.max ) == 0, , drop = FALSE ]
  pos <- drop.saturated( stained.data )
  neg <- drop.saturated( unstained.data )

  pos.spectral <- pos[ , spec.ch, drop = FALSE ]
  pos.scatter  <- pos[ , scatter.channels, drop = FALSE ]
  neg.spectral <- neg[ , spec.ch, drop = FALSE ]
  neg.scatter  <- neg[ , scatter.channels, drop = FALSE ]

  # positive events: above the 99.5th percentile of the raw unstained peak channel
  raw.threshold <- stats::quantile( neg.spectral[ , peak.channel ], raw.threshold.quantile )
  pos.idx <- which( pos.spectral[ , peak.channel ] > raw.threshold )

  if ( length( pos.idx ) > n.cells * 2 ) {
    top <- order( pos.spectral[ pos.idx, peak.channel ], decreasing = TRUE )[ seq_len( n.cells * 2 ) ]
    pos.idx <- pos.idx[ top ]
  }
  if ( length( pos.idx ) < min.events )
    stop( "Only ", length( pos.idx ), " events above the raw positivity threshold; ",
          "cannot estimate spillover spread.", call. = FALSE )

  # scatter-matched universal-negative subtraction
  knn.idx <- FNN::knnx.index(
    data  = neg.scatter,
    query = pos.scatter[ pos.idx, , drop = FALSE ],
    k     = k.neighbors
  )
  bg <- matrix( 0, nrow( knn.idx ), length( spec.ch ), dimnames = list( NULL, spec.ch ) )
  for ( ki in seq_len( k.neighbors ) )
    bg <- bg + neg.spectral[ knn.idx[ , ki ], , drop = FALSE ]
  bg <- bg / k.neighbors
  pos.corrected <- pos.spectral[ pos.idx, , drop = FALSE ] - bg

  # AF principal directions from the unstained control
  af.data <- neg.spectral
  if ( nrow( af.data ) > asp$gate.downsample.n.cells ) {
    set.seed( asp$bird.seed )
    af.data <- af.data[ sample( nrow( af.data ), asp$gate.downsample.n.cells ), , drop = FALSE ]
  }
  sv <- svd( af.data, nu = 0, nv = n.af.pcs )
  af.pcs <- t( sv$v )
  colnames( af.pcs ) <- spec.ch
  rownames( af.pcs ) <- paste0( "AFPC", seq_len( nrow( af.pcs ) ) )

  orig <- spectra[ fluorophore, , drop = FALSE ]
  af.cosine <- max( abs( cosine.similarity( rbind( orig, af.pcs ) )[ 1, -1 ] ) )
  af.collinear <- af.cosine >= af.collinear.threshold

  # low-rank fit: AF PCs + fluorophore (or fluorophore alone if collinear)
  design <- if ( af.collinear ) orig else rbind( af.pcs, orig )
  fit <- unmix.ols( pos.corrected, design )
  b.h <- as.numeric( fit[ , ncol( fit ) ] )
  resid.lr <- pos.corrected - fit %*% design

  # Residual Model line per channel: the package's estimator
  residual.model.line <- utils::getFromNamespace( ".residual.model.line", "AutoSpectral" )
  coef <- residual.model.line(
    b.h, resid.lr, spectra,
    spread.estimator = spread.estimator,
    n.bins = spread.n.bins, min.bin.events = spread.min.bin.events,
    huber.k = huber.k, huber.max.iter = huber.max.iter
  )

  intercept <- coef[ 1, ]
  spread <- coef[ 2, ]
  names( intercept ) <- names( spread ) <- rownames( spectra )

  spread[ !is.na( spread ) & spread < 0 ] <- 0
  intercept[ !is.na( intercept ) & intercept < 0 ] <- 0
  spread[ fluorophore ] <- NA_real_
  intercept[ fluorophore ] <- NA_real_

  list(
    spread       = spread,
    intercept    = intercept,
    n.events     = length( b.h ),
    af.collinear = af.collinear,
    af.cosine    = af.cosine,
    estimator    = spread.estimator
  )
}
