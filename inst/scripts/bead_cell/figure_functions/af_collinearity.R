# af_collinearity.R
#
# Is the spectral error left by subtracting the wrong cell type's background
# collinear with autofluorescence (AF)?
#
# Background. For a spectrum extracted as the slope of every detector on the
# peak detector over positive + negative events, replacing the correct
# negative (median background B_r) by another (B_w) changes the spectrum by
# approximately
#     delta = S_w - S_r  ~  (B_r - B_w) / (P_peak - B_peak)  -  c * S
# i.e. a scaled difference between the two backgrounds, minus a small
# multiple of the spectrum itself (from the change in the peak-detector
# background). The test is therefore whether the difference between the AF
# of two cell types lies along the main directions of AF variation in the
# unstained cells, rather than in some other spectral direction.
#
# Statistics for one delta (all detectors):
#   cos.af      cosine similarity (signed) between delta and the median AF
#               spectrum of the unstained cells;
#   r2.af       fraction of the squared norm of delta that lies in the
#               subspace spanned by the first k AF principal components
#               (uncentred SVD of the unstained cells):  ||P delta||^2/||delta||^2;
#   null        r2 of every spectrum in the cytometer's reference library
#               (except the stained fluorophore) in the same subspace, i.e.
#               how much of an arbitrary dye spectrum AF can absorb;
#   p.value     (1 + #{library r2 >= r2.af}) / (1 + N);
#   rank.af     rank of AF among AF + the panel fluorophores by |cosine| with
#               delta (1 = AF is the most collinear);
#   most.collinear  panel fluorophore (not AF, not the stained dye) with the
#               largest |cosine| with delta.
# Bootstrap percentile intervals for cos.af and r2.af come from
# bootstrap.delta() (fig3_extraction.R).

#' AF reference: median spectrum and uncentred principal components
#'
#' @param unstained Matrix of unstained events (spectral columns used).
#' @param spectral.channels Detector names.
#' @param k Number of components.
#' @param max.events Downsample to this many events for the SVD.
#' @return list(median = named vector (max = 1), basis = k x detectors matrix
#'   with orthonormal rows, var.explained = fraction of the total squared norm
#'   per component).
af.reference <- function( unstained, spectral.channels, k = 3L, max.events = 2e4, seed = 42 ) {
  x <- unstained[ , spectral.channels, drop = FALSE ]
  med <- apply( x, 2, stats::median )
  med <- med / max( med )
  if ( nrow( x ) > max.events ) {
    set.seed( seed )
    x <- x[ sample( nrow( x ), max.events ), , drop = FALSE ]
  }
  sv <- svd( x, nu = 0, nv = k )
  basis <- t( sv$v )
  # orient each component so that its largest-magnitude element is positive
  basis <- t( apply( basis, 1, function( v ) v * sign( v[ which.max( abs( v ) ) ] ) ) )
  dimnames( basis ) <- list( paste0( "AF.PC", seq_len( k ) ), spectral.channels )
  list(
    median = med,
    basis = basis,
    var.explained = ( sv$d^2 / sum( sv$d^2 ) )[ seq_len( k ) ]
  )
}

.cos <- function( a, b ) sum( a * b ) / sqrt( sum( a^2 ) * sum( b^2 ) )

.r2.subspace <- function( v, basis ) {
  # basis rows orthonormal
  proj <- crossprod( basis, basis %*% v )
  sum( proj^2 ) / sum( v^2 )
}

#' Collinearity statistics for one delta
#'
#' @param delta Named numeric vector (detectors).
#' @param af List from af.reference().
#' @param library.spectra Matrix (dyes x detectors), the reference library.
#' @param panel.spectra Matrix (dyes x detectors), the unmixing panel WITHOUT AF.
#' @param stained Name of the stained fluorophore (excluded from the null and
#'   from the most-collinear search).
#' @return list of scalars + the null vector and the per-dye cosines.
af.collinearity <- function( delta, af, library.spectra, panel.spectra, stained ) {

  ch <- names( delta )
  basis <- af$basis[ , ch, drop = FALSE ]

  cos.af <- .cos( delta, af$median[ ch ] )
  r2.af  <- .r2.subspace( delta, basis )

  lib <- library.spectra[ rownames( library.spectra ) != stained, ch, drop = FALSE ]
  null.r2 <- apply( lib, 1, .r2.subspace, basis = basis )
  p.value <- ( 1 + sum( null.r2 >= r2.af ) ) / ( 1 + length( null.r2 ) )

  pan <- panel.spectra[ rownames( panel.spectra ) != stained, ch, drop = FALSE ]
  cos.panel <- apply( pan, 1, .cos, b = delta )
  rank.af <- 1 + sum( abs( cos.panel ) > abs( cos.af ) )
  mc <- names( cos.panel )[ which.max( abs( cos.panel ) ) ]

  list(
    cos.af = cos.af,
    r2.af = r2.af,
    p.value = p.value,
    null.r2 = null.r2,
    null.median = stats::median( null.r2 ),
    rank.af = rank.af,
    n.compared = length( cos.panel ) + 1,
    most.collinear = mc,
    cos.most.collinear = unname( cos.panel[ mc ] ),
    cos.panel = cos.panel
  )
}

#' Bootstrap intervals for cos.af and r2.af
af.collinearity.ci <- function( delta.boot, af, level = 0.95 ) {
  ch <- colnames( delta.boot )
  basis <- af$basis[ , ch, drop = FALSE ]
  cos.b <- apply( delta.boot, 1, .cos, b = af$median[ ch ] )
  r2.b  <- apply( delta.boot, 1, .r2.subspace, basis = basis )
  a <- ( 1 - level ) / 2
  c(
    cos.af.lo = unname( stats::quantile( cos.b, a ) ),
    cos.af.hi = unname( stats::quantile( cos.b, 1 - a ) ),
    r2.af.lo  = unname( stats::quantile( r2.b, a ) ),
    r2.af.hi  = unname( stats::quantile( r2.b, 1 - a ) ),
    n.boot    = nrow( delta.boot )
  )
}
