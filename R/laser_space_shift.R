# laser_space_shift.R
#
# Instrument-independent description of how a fluorophore's spectrum changes
# between two particle types (e.g. cells and compensation beads). Detectors
# differ between cytometers but the excitation lasers are shared, so each
# laser's bank of detectors is treated as one axis ("laser space").
# Moved into the package from the bead-vs-cell manuscript figure functions
# (inst/scripts/bead_cell/figure_functions); the spider-chart plotting stays
# there.

#' @title Cell-to-Bead Spectral Shift in Laser Space
#'
#' @description
#' Describes the change between two spectra of the same fluorophore (for
#' example, measured on cells and on compensation beads) per excitation laser,
#' so that cytometers with different detectors can be compared.
#'
#' 1. Both spectra are restricted to the detectors on `lasers` and normalised
#'    to unit Euclidean (L2) norm over those detectors.
#' 2. The bank norm of laser \eqn{L} is \eqn{n_L = \|s_L\|_2}, the L2 norm of
#'    the spectrum over the detectors of that laser. Because the spectrum has
#'    unit norm, \eqn{\sum_L n_L^2 = 1}, so the bank norms form a unit vector
#'    in laser space.
#' 3. The shift per laser is either the centred log2 ratio (`metric = "clr"`),
#'    \eqn{\log_2(n_{L,beads} / n_{L,cells})} minus its mean over lasers, or
#'    the difference \eqn{n_{L,beads} - n_{L,cells}} (`metric = "delta"`).
#'    A gain applied to a whole laser bank multiplies the cell and bead bank
#'    norms alike and cancels in the ratio, and any common scaling cancels in
#'    the centring, so "clr" shifts can be compared between cytometers. The
#'    centred log ratio is the standard transform for compositional data
#'    (Aitchison, 1986).
#'
#' The spectral angle between the two unit spectra is decomposed into the
#' angle between the two laser-space vectors (`angle.between`; redistribution
#' of signal between lasers, at most `angle.total` by the Cauchy-Schwarz
#' inequality) and, per laser, the angle between the two bank sub-spectra
#' (`angle.within`; change of shape within a laser's emission profile).
#'
#' @param cells,beads Named numeric vectors or one-row matrices over the same
#'   detectors: the reference spectrum (`cells`) and the test spectrum
#'   (`beads`).
#' @param detector.lasers Detector-to-laser assignment: a named character
#'   vector (names = detectors) or a data frame with columns `Detector` and
#'   `Laser`. If `NULL` (default), the lasers are looked up in
#'   `cytometer_database.csv` for the cytometer in `asp`.
#' @param lasers Character vector of the lasers to use (the axes), in display
#'   order. Every laser must have at least one detector.
#' @param metric `"clr"` (centred log2 ratio, default) or `"delta"`.
#' @param asp AutoSpectral parameter list, used only when `detector.lasers` is
#'   `NULL`.
#'
#' @return A list with `shift` (named by laser), `cell.norms`, `bead.norms`,
#'   `angle.total`, `angle.between` (degrees) and `angle.within` (degrees,
#'   named by laser).
#'
#' @references Aitchison J (1986). The Statistical Analysis of Compositional
#'   Data. Chapman and Hall, London.
#'
#' @seealso [laser.shift.concordance()], [assess.mismatch.angle()]
#'
#' @export
laser.space.shift <- function(
    cells,
    beads,
    detector.lasers = NULL,
    lasers,
    metric = c( "clr", "delta" ),
    asp = NULL
) {

  metric <- match.arg( metric )
  c.v <- .laser.shift.as.vector( cells )
  b.v <- .laser.shift.as.vector( beads )
  if ( !identical( names( c.v ), names( b.v ) ) )
    stop( "Cell and bead spectra must have identical detector names.", call. = FALSE )

  if ( is.null( detector.lasers ) ) {
    detector.lasers <- .detector.lasers( names( c.v ), asp = asp )
    if ( is.null( detector.lasers ) )
      stop( "Could not assign lasers to the detectors; supply `detector.lasers` ",
            "or an `asp` for this cytometer.", call. = FALSE )
  }
  if ( is.data.frame( detector.lasers ) )
    detector.lasers <- stats::setNames( as.character( detector.lasers$Laser ),
                                        as.character( detector.lasers$Detector ) )

  las <- unname( detector.lasers[ names( c.v ) ] )
  if ( anyNA( las ) )
    warning( sum( is.na( las ) ), " detector(s) without a laser assignment were ignored: ",
             paste( names( c.v )[ is.na( las ) ], collapse = ", " ), call. = FALSE )
  missing.lasers <- setdiff( lasers, las )
  if ( length( missing.lasers ) > 0 )
    stop( "Laser(s) not present on this instrument: ", paste( missing.lasers, collapse = ", " ),
          call. = FALSE )

  keep <- !is.na( las ) & las %in% lasers
  c.v <- c.v[ keep ]; b.v <- b.v[ keep ]; las <- las[ keep ]
  c.v <- c.v / sqrt( sum( c.v^2 ) )
  b.v <- b.v / sqrt( sum( b.v^2 ) )

  bank.norm <- function( v ) vapply( lasers, function( L ) sqrt( sum( v[ las == L ]^2 ) ), numeric( 1 ) )
  n.c <- bank.norm( c.v )
  n.b <- bank.norm( b.v )

  if ( any( c( n.c, n.b ) < 1e-6 ) && metric == "clr" )
    warning( "A laser bank has (near-)zero signal; its log ratio is unstable.", call. = FALSE )

  shift <- if ( metric == "clr" ) {
    lr <- log2( n.b / n.c )
    lr - mean( lr )
  } else {
    n.b - n.c
  }

  within <- vapply( lasers, function( L ) {
    i <- las == L
    .laser.shift.acos.deg( sum( c.v[ i ] * b.v[ i ] ) / ( n.c[[ L ]] * n.b[[ L ]] ) )
  }, numeric( 1 ) )

  list(
    shift         = shift,
    cell.norms    = n.c,
    bead.norms    = n.b,
    angle.total   = .laser.shift.acos.deg( sum( c.v * b.v ) ),
    angle.between = .laser.shift.acos.deg( sum( n.c * n.b ) ),
    angle.within  = within
  )
}

#' @title Directional Concordance of Laser-Space Shifts Across Cytometers
#'
#' @description
#' Measures whether the laser-space shift of a fluorophore (from
#' [laser.space.shift()]) points in the same direction on different
#' cytometers. The statistic is the mean pairwise cosine similarity of the
#' cytometers' shift vectors. Its null distribution is obtained by permuting
#' the laser labels independently within each cytometer. P = (1 + number of
#' permutations with a mean cosine at least the observed value) /
#' (1 + `n.perm`), so the smallest attainable P is 1 / (1 + `n.perm`).
#'
#' @param shift.matrix Numeric matrix, cytometers (rows, named) by lasers
#'   (columns).
#' @param n.perm Number of random permutations. Default `10000`.
#' @param seed Random seed, set with `set.seed()` before permuting.
#'   Default `42`.
#'
#' @return A list with `pairwise` (cosine matrix), `mean.cosine`, `p.value`,
#'   `n.perm`, `leave.one.out` (cosine of each cytometer's shift with the mean
#'   shift of the others) and `null` (the permutation distribution).
#'
#' @seealso [laser.space.shift()]
#'
#' @export
laser.shift.concordance <- function( shift.matrix, n.perm = 10000, seed = 42 ) {

  m <- as.matrix( shift.matrix )
  if ( nrow( m ) < 2 ) stop( "Need at least two instruments.", call. = FALSE )

  cos.matrix <- function( x ) { u <- x / sqrt( rowSums( x^2 ) ); tcrossprod( u ) }
  mean.pairwise <- function( x ) { cm <- cos.matrix( x ); mean( cm[ upper.tri( cm ) ] ) }

  observed <- mean.pairwise( m )

  set.seed( seed )
  null <- vapply( seq_len( n.perm ), function( i ) {
    mean.pairwise( t( apply( m, 1, sample ) ) )
  }, numeric( 1 ) )
  p.value <- ( 1 + sum( null >= observed ) ) / ( 1 + n.perm )

  loo <- vapply( seq_len( nrow( m ) ), function( i )
    cosine.similarity( rbind( m[ i, ], colMeans( m[ -i, , drop = FALSE ] ) ) )[ 1, 2 ], numeric( 1 ) )
  names( loo ) <- rownames( m )

  pairwise <- cos.matrix( m )
  dimnames( pairwise ) <- list( rownames( m ), rownames( m ) )

  list(
    pairwise      = pairwise,
    mean.cosine   = observed,
    p.value       = p.value,
    n.perm        = n.perm,
    leave.one.out = loo,
    null          = null
  )
}

# first row of a matrix/data frame as a named vector; vectors unchanged
.laser.shift.as.vector <- function( x ) {
  if ( is.matrix( x ) || is.data.frame( x ) ) {
    v <- as.numeric( x[ 1, ] ); names( v ) <- colnames( x ); v
  } else x
}

.laser.shift.acos.deg <- function( x ) acos( pmin( 1, pmax( -1, x ) ) ) * 180 / pi
