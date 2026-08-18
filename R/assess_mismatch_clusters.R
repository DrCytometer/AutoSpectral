# assess_mismatch_clusters.R

#' @title Cluster-Based Permutation Test for Spectral Mismatch Regions
#'
#' @description
#' Identifies contiguous regions of the detector array where the mismatch
#' between two particle types' spectra (e.g. beads vs. cells) is larger than
#' expected by chance, while accounting for two structural features of the
#' data: (1) adjacent detectors are correlated because a fluorophore's
#' emission is a smooth curve, so per-detector tests cannot be treated as
#' independent, and (2) a fluorophore can only show mismatch at a detector
#' where it has appreciable signal in the first place, so raw per-detector
#' averages are confounded by how many/which fluorophores peak there.
#'
#' Both issues are handled together: mismatch is only evaluated at
#' fluorophore/detector pairs where the fluorophore has on-peak signal
#' (`peak.threshold`), a per-detector statistic is computed only from those
#' on-peak observations, and significance is assessed with a sign-flip
#' permutation test using cluster mass (summed |statistic| over contiguous
#' detectors) as the test statistic. Because permutations flip the sign of
#' whole fluorophore traces rather than shuffling individual detector
#' values, the detector-to-detector correlation structure already present
#' in the data is preserved in the null distribution, so no independence
#' assumption across detectors is required.
#'
#' @param dist.mat Numeric matrix of signed per-detector differences
#'   (fluorophore x detector), e.g. \code{comparison[[b]]$Distance} from
#'   \code{run.bead.cell.comparison()} or the output of \code{bead.cell.dist()}.
#'   Column order is assumed to be true spectral/detector order (as in the
#'   package's per-cytometer reference libraries); contiguity for clustering
#'   is defined by that order, not by column name.
#' @param ref.spectra Numeric matrix (fluorophore x detector) giving the
#'   reference particle type's spectral intensity, used only to build the
#'   on-peak mask. Must share row and column names with \code{dist.mat}
#'   (extra rows/columns are dropped; a warning is issued if the shared set
#'   is smaller than either input).
#' @param peak.threshold Numeric scalar in `[0, 1]`. A fluorophore/detector
#'   pair is included in the test if \code{ref.spectra[f, d] >=
#'   peak.threshold * max(ref.spectra[f, ])}. Defaults to \code{0.1} (10\%
#'   of that fluorophore's peak height).
#' @param min.n Integer scalar. Minimum number of on-peak fluorophores
#'   required at a detector for it to receive a test statistic; detectors
#'   with fewer contributing fluorophores are set to \code{NA} and cannot
#'   seed or extend a cluster. Defaults to \code{3L}.
#' @param cluster.threshold Numeric scalar. Primary (per-detector)
#'   threshold on \code{abs(t-statistic)} used to form candidate clusters,
#'   in the usual cluster-based-permutation sense: this sets what counts as
#'   "locally unusual" before cluster mass is compared across permutations.
#'   Defaults to \code{1.96} (approx. two-sided alpha = 0.05 for a single
#'   detector, before cluster-level correction).
#' @param n.perm Integer scalar. Number of sign-flip permutations used to
#'   build the null distribution of the maximum cluster mass. Defaults to
#'   \code{2000L}.
#' @param bird.seed Integer scalar or \code{NULL}. Seed for the permutation
#'   draws, for reproducibility. Defaults to \code{NULL} (no seeding).
#'
#' @details
#' At each detector \code{d}, let \code{x} be the on-peak entries of
#' \code{dist.mat[, d]} (those \code{f} with \code{mask[f, d] == TRUE}). The
#' per-detector statistic is the one-sample t-statistic
#' \code{mean(x) / (sd(x) / sqrt(length(x)))}, testing the null hypothesis
#' that the reference and comparison particle types produce the same
#' spectrum at that detector, restricted to fluorophores that actually
#' fluoresce there. Detectors are then thresholded at
#' \code{cluster.threshold} and contiguous runs of exceedance are collapsed
#' into clusters, each summarized by its cluster mass (sum of \code{abs(t)}
#' across its member detectors).
#'
#' The null distribution is built by, for each of \code{n.perm} iterations,
#' independently flipping the sign of every fluorophore's entire on-peak
#' row, recomputing the per-detector statistic and cluster masses under
#' that permutation, and retaining only the single largest cluster mass
#' observed. This max-statistic approach controls the family-wise error
#' rate across the whole detector array without assuming detectors are
#' independent.
#'
#' Each observed cluster's p-value is the proportion of permutation draws
#' whose maximum cluster mass met or exceeded that cluster's mass.
#'
#' @return A list with components:
#' \describe{
#'   \item{stat}{Named numeric vector (length \code{ncol(dist.mat)}) of the
#'     observed per-detector t-statistic, \code{NA} where \code{min.n} was
#'     not met.}
#'   \item{mask}{Logical matrix (fluorophore x detector), the on-peak mask
#'     used.}
#'   \item{n.detector}{Named integer vector of the number of on-peak
#'     fluorophores contributing to each detector's statistic.}
#'   \item{clusters}{Data frame, one row per observed cluster, with columns
#'     \code{start}, \code{end} (detector names), \code{n.detectors},
#'     \code{mass}, and \code{p.value}. Empty (0-row) if no detector met
#'     \code{cluster.threshold}.}
#'   \item{null.max.mass}{Numeric vector of length \code{n.perm}, the
#'     permutation null distribution of the maximum cluster mass.}
#' }
#'
#' @export

assess.mismatch.clusters <- function(
    dist.mat,
    ref.spectra,
    peak.threshold    = 0.1,
    min.n             = 3L,
    cluster.threshold = 1.96,
    n.perm            = 2000L,
    bird.seed         = NULL
) {

  if ( !is.matrix( dist.mat ) || !is.numeric( dist.mat ) )
    stop( "`dist.mat` must be a numeric matrix.", call. = FALSE )

  if ( !is.matrix( ref.spectra ) || !is.numeric( ref.spectra ) )
    stop( "`ref.spectra` must be a numeric matrix.", call. = FALSE )

  common.fluor <- intersect( rownames( dist.mat ), rownames( ref.spectra ) )
  common.det   <- intersect( colnames( dist.mat ), colnames( ref.spectra ) )

  if ( length( common.fluor ) == 0 || length( common.det ) == 0 )
    stop(
      "No shared fluorophores/detectors between `dist.mat` and ",
      "`ref.spectra` -- check rownames/colnames.", call. = FALSE
    )

  if ( length( common.fluor ) < nrow( dist.mat ) ||
       length( common.det )   < ncol( dist.mat ) )
    warning(
      "Restricting to ", length( common.fluor ), " shared fluorophore(s) and ",
      length( common.det ), " shared detector(s).", call. = FALSE
    )

  # preserve dist.mat's detector (spectral) order -- contiguity is defined
  # by column order, so this must be the real spectral order, not
  # alphabetical or match-derived order
  common.det <- colnames( dist.mat )[ colnames( dist.mat ) %in% common.det ]

  dist.mat    <- dist.mat[ common.fluor, common.det, drop = FALSE ]
  ref.spectra <- ref.spectra[ common.fluor, common.det, drop = FALSE ]

  n.fluor <- nrow( dist.mat )

  # -- on-peak mask -----------------------------------------------------------
  row.max <- apply( ref.spectra, 1, max, na.rm = TRUE )
  mask <- ref.spectra >= ( peak.threshold * row.max )
  mask[ is.na( mask ) ] <- FALSE

  masked.delta <- dist.mat
  masked.delta[ !mask ] <- NA

  # -- per-detector statistic --------------------------------------------------
  .detector.stat <- function( x ) {
    x <- x[ !is.na( x ) ]
    n <- length( x )
    if ( n < min.n ) return( c( stat = NA_real_, n = n ) )
    s <- stats::sd( x )
    if ( !is.finite( s ) || s == 0 ) return( c( stat = NA_real_, n = n ) )
    c( stat = mean( x ) / ( s / sqrt( n ) ), n = n )
  }

  obs <- apply( masked.delta, 2, .detector.stat )
  obs.stat <- obs[ "stat", ]
  obs.n    <- obs[ "n", ]
  names( obs.stat ) <- names( obs.n ) <- common.det

  # -- cluster identification --------------------------------------------------
  .find.clusters <- function( stat.vec, threshold ) {
    exceed <- which( !is.na( stat.vec ) & abs( stat.vec ) >= threshold )
    if ( length( exceed ) == 0 )
      return( list() )
    breaks <- c( 0, which( diff( exceed ) != 1 ), length( exceed ) )
    lapply( seq_len( length( breaks ) - 1 ), function( i ) {
      exceed[ ( breaks[ i ] + 1 ):breaks[ i + 1 ] ]
    } )
  }

  .cluster.mass <- function( idx, stat.vec ) sum( abs( stat.vec[ idx ] ) )

  obs.clusters <- .find.clusters( obs.stat, cluster.threshold )

  # -- permutation null of max cluster mass ------------------------------------
  if ( !is.null( bird.seed ) ) set.seed( bird.seed )

  null.max.mass <- vapply( seq_len( n.perm ), function( i ) {

    flip <- sample( c( -1, 1 ), n.fluor, replace = TRUE )
    perm.delta <- masked.delta * flip

    perm.stat <- apply( perm.delta, 2, function( x ) {
      x <- x[ !is.na( x ) ]
      n <- length( x )
      if ( n < min.n ) return( NA_real_ )
      s <- stats::sd( x )
      if ( !is.finite( s ) || s == 0 ) return( NA_real_ )
      mean( x ) / ( s / sqrt( n ) )
    } )

    perm.clusters <- .find.clusters( perm.stat, cluster.threshold )
    if ( length( perm.clusters ) == 0 ) return( 0 )
    max( vapply( perm.clusters, .cluster.mass, numeric( 1 ), stat.vec = perm.stat ) )

  }, numeric( 1 ) )

  # -- assemble cluster results -------------------------------------------------
  if ( length( obs.clusters ) == 0 ) {
    clusters.df <- data.frame(
      start = character(0), end = character(0),
      n.detectors = integer(0), mass = numeric(0), p.value = numeric(0),
      stringsAsFactors = FALSE
    )
  } else {
    clusters.df <- do.call( rbind, lapply( obs.clusters, function( idx ) {
      mass <- .cluster.mass( idx, obs.stat )
      p    <- ( 1 + sum( null.max.mass >= mass ) ) / ( n.perm + 1 )
      data.frame(
        start       = common.det[ idx[ 1 ] ],
        end         = common.det[ idx[ length( idx ) ] ],
        n.detectors = length( idx ),
        mass        = mass,
        p.value     = p,
        stringsAsFactors = FALSE
      )
    } ) )
    rownames( clusters.df ) <- NULL
  }

  list(
    stat          = obs.stat,
    mask          = mask,
    n.detector    = obs.n,
    clusters      = clusters.df,
    null.max.mass = null.max.mass
  )
}
