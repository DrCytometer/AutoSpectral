# noise_model_groups.R

# Building blocks of the groups `estimate.noise.model()` fits: event splits,
# projection designs and bead negatives.

# Sorted training and held-out indices; the held-out set capped at max.test.
.holdout.split <- function( n, train.frac, seed, max.test = NULL ) {
  if ( !( train.frac > 0 && train.frac < 1 ) )
    stop( "The training fraction must be between 0 and 1.", call. = FALSE )
  set.seed( seed )
  order.n <- sample.int( n )
  n.train <- round( n * train.frac )
  train <- order.n[ seq_len( n.train ) ]
  test  <- order.n[ setdiff( seq_len( n ), seq_len( n.train ) ) ]
  if ( !is.null( max.test ) && length( test ) > max.test ) test <- test[ seq_len( max.test ) ]
  list( train = sort( train ), test = sort( test ) )
}

# Sorted indices of at most cap of n events.
.subsample.idx <- function( n, cap, seed ) {
  if ( n <= cap ) return( seq_len( n ) )
  set.seed( seed )
  sort( sample.int( n, cap ) )
}

# Unsaturated events above threshold in the peak detector.
.positive.rows <- function( raw, peak, threshold, ceiling ) {
  unsaturated <- rowSums( raw >= 0.999 * ceiling ) == 0
  which( unsaturated & raw[ , peak ] > threshold )
}

# Leading right singular vectors of the uncentred events, unit length and
# signed so that each sums positive: the first follows the mean AF spectrum,
# the rest its main variation.
.af.basis.rows <- function( raw, n.rows = 4L, max.events = 20000L, seed = 0L ) {
  raw <- as.matrix( raw )
  raw <- raw[ .subsample.idx( nrow( raw ), max.events, seed ), , drop = FALSE ]
  n.keep <- max( 1L, min( n.rows, nrow( raw ), ncol( raw ) ) )
  rows <- t( svd( raw, nu = 0, nv = n.keep )$v )
  sgn  <- sign( rowSums( rows ) )
  sgn[ sgn == 0 ] <- 1
  rows <- rows * sgn
  colnames( rows ) <- colnames( raw )
  rows
}

# The spectrum restricted to each laser's detectors, one row per laser where
# that part and the rest each reach split.min of the spectrum's peak.
.laser.rows <- function( spectrum, lasers, split.min = 0.01 ) {
  spectrum <- as.numeric( spectrum )
  las  <- as.character( lasers )
  peak <- max( abs( spectrum ) )
  rows <- list()
  for ( laser in unique( las ) ) {
    block <- ifelse( las == laser, spectrum, 0 )
    if ( max( abs( block ) ) >= split.min * peak && max( abs( spectrum - block ) ) >= split.min * peak )
      rows[[ length( rows ) + 1L ]] <- block
  }
  if ( length( rows ) ) do.call( rbind, rows ) else NULL
}

# Projection design of one control: its reference spectrum (row 1, the
# binning variable), its per-laser rows when lasers are given, the mean
# variant delta when it exceeds mean.delta.min of the peak, the variant basis
# (an entry of build.variant.basis()) and the AF rows.
.control.design <- function( spectrum, variant.basis = NULL, af.rows = NULL, mean.delta.min = 1e-3,
                             lasers = NULL, split.min = 0.01 ) {
  rows <- list( as.numeric( spectrum ) )
  if ( !is.null( lasers ) ) rows <- c( rows, list( .laser.rows( spectrum, lasers, split.min ) ) )
  if ( !is.null( variant.basis ) ) {
    if ( max( abs( variant.basis$mean.delta ) ) > mean.delta.min )
      rows <- c( rows, list( as.numeric( variant.basis$mean.delta ) ) )
    rows <- c( rows, list( variant.basis$basis ) )
  }
  if ( !is.null( af.rows ) ) rows <- c( rows, list( af.rows ) )
  design <- do.call( rbind, Filter( Negate( is.null ), rows ) )
  dimnames( design ) <- NULL
  design
}

# Otsu split of values.
.otsu.threshold <- function( values, n.bins = 256L ) {
  rng   <- range( values )
  edges <- seq( rng[ 1 ], rng[ 2 ], length.out = n.bins + 1L )
  idx   <- pmin( findInterval( values, edges ), n.bins )
  hist  <- tabulate( idx, nbins = n.bins )
  centres <- 0.5 * ( edges[ -1 ] + edges[ -( n.bins + 1L ) ] )
  w0 <- cumsum( hist )
  w1 <- w0[ n.bins ] - w0
  m0 <- cumsum( hist * centres )
  between <- w0 * w1 * ( m0 / w0 - ( m0[ n.bins ] - m0 ) / w1 )^2
  between[ !is.finite( between ) ] <- -1
  centres[ which.max( between ) ]
}

# Dark level of bead events per detector: the median, over bead controls
# whose dye emits at most hard.frac of its peak there, of that control's
# median signal in the detector. `medians` is a list of
# list(median = (D), fluor = row of spectra). NA where no control is that dark.
.bead.reference.levels <- function( medians, spectra, hard.frac = 0.005 ) {
  d <- ncol( spectra )
  vapply( seq_len( d ), function( j ) {
    v <- unlist( lapply( medians, function( m ) if ( spectra[ m$fluor, j ] <= hard.frac ) m$median[ j ] ) )
    if ( length( v ) ) stats::median( v ) else NA_real_
  }, numeric( 1 ) )
}

# Unstained beads inside a bead single-stain control: an Otsu split of log10
# of the peak detector, accepted when the lower population's median is at
# least `factor` times below the upper one's, it has min.events events, and
# its median is within `tolerance` times the dark level other bead controls
# show in this detector (ignored when NA). Returns list(mask or NULL, reason).
.internal.bead.negative <- function( events, peak, reference.level, factor = 10, tolerance = 3,
                                     min.events = 200L ) {
  p <- log10( pmax( events[ , peak ], 1 ) )
  if ( length( p ) < 2L * min.events || diff( range( p ) ) < log10( factor ) )
    return( list( mask = NULL, reason = "no bimodal peak channel" ) )
  cut <- .otsu.threshold( p )
  low <- p < cut
  if ( sum( low ) < min.events || sum( !low ) < min.events )
    return( list( mask = NULL, reason = "too few events on one side of the split" ) )
  lo.med <- stats::median( events[ low, peak ] )
  hi.med <- stats::median( events[ !low, peak ] )
  if ( hi.med < factor * max( lo.med, 1 ) )
    return( list( mask = NULL, reason = sprintf( "populations only %.1fx apart", hi.med / max( lo.med, 1 ) ) ) )
  if ( is.finite( reference.level ) && lo.med > tolerance * max( reference.level, 1 ) )
    return( list( mask = NULL, reason = sprintf(
      "lower population (%.0f) above the bead dark level (%.0f)", lo.med, reference.level ) ) )
  low <- low & events[ , peak ] < hi.med / factor
  if ( sum( low ) < min.events )
    return( list( mask = NULL, reason = "too few events below the positives / factor" ) )
  list( mask = low, reason = sprintf( "%d events, %.0fx below the positives", sum( low ), hi.med / max( lo.med, 1 ) ) )
}
