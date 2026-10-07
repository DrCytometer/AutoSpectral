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
#  4. project the residual through the FULL panel with unmix.ols();
#  5. for each unmixed channel k, regress proj_k^2 on b.h with a Huber-robust
#     line: Var(unmixed_k | b.h) = intercept_k + slope_k * b.h. The slope is
#     SS(fluorophore, k) in the units used throughout AutoSpectral.
#
# Differences from the package: (a) it works on events already gated by the
# caller (the package reads the control ungated), and (b) it does not need a
# control file or the whole panel of controls.

# Huber IRLS line fit, used only if the package's .fix.huber.slope() (needs
# AutoSpectralRcpp) is unavailable. Returns c(intercept, slope).
.huber.line.fallback <- function( x, y, k = 1.345, max.iter = 100L, tol = 1e-4 ) {
  X <- cbind( 1, x )
  beta <- stats::lm.fit( X, y )$coefficients
  for ( i in seq_len( max.iter ) ) {
    r <- y - drop( X %*% beta )
    s <- stats::mad( r )
    if ( !is.finite( s ) || s == 0 ) break
    w <- pmin( 1, k * s / abs( r ) )
    new.beta <- stats::lm.wfit( X, y, w )$coefficients
    converged <- max( abs( new.beta - beta ) / ( abs( beta ) + 1e-12 ) ) < tol
    beta <- new.beta
    if ( converged ) break
  }
  as.numeric( beta )
}

.huber.line <- function( x, y, k, max.iter ) {
  pkg.fn <- tryCatch(
    utils::getFromNamespace( ".fix.huber.slope", "AutoSpectral" ),
    error = function( e ) NULL
  )
  if ( !is.null( pkg.fn ) ) {
    res <- tryCatch( pkg.fn( x, y, k = k, max.iter = max.iter ), error = function( e ) NULL )
    if ( !is.null( res ) ) return( as.numeric( res ) )
  }
  .huber.line.fallback( x, y, k = k, max.iter = max.iter )
}

#' @param stained.data,unstained.data Matrices of GATED raw events containing
#'   the spectral channels (`colnames(spectra)`) and `scatter.channels`.
#' @param spectra Full-panel spectral matrix (fluorophores x detectors) used
#'   for unmixing (the correct, cell-derived matrix).
#' @param fluorophore Row of `spectra` for the stained control.
#' @param peak.channel Peak detector of the fluorophore.
#' @param asp AutoSpectral parameter list.
#' @param scatter.channels FSC/SSC channel names (default
#'   `asp$default.scatter.parameter`).
#' @return list(spread = named numeric slope per unmixed channel (NA for the
#'   fluorophore itself; negative slopes set to 0), intercept, n.events,
#'   af.collinear, af.cosine).
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
    huber.k = 1.345,
    huber.max.iter = 100L,
    min.events = 50L
) {

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

  # project the residual through the full panel and regress squared value on b.h
  proj <- unmix.ols( resid.lr, spectra )
  coef <- apply( proj, 2, function( p ) .huber.line( b.h, p^2, huber.k, huber.max.iter ) )

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
    af.cosine    = af.cosine
  )
}
