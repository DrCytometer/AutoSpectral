# unmix_af_fwl.r

#' @title Unmix With Per-Cell Autofluorescence By Frisch-Waugh
#'
#' @description
#' Solves the joint least-squares problem `y = t(spectra) %*% f + k * af` for
#' every cell at once, where each cell may use a different autofluorescence
#' spectrum. The Frisch-Waugh-Lovell decomposition splits the joint solve into
#' two precomputable library projections plus one pass over the data, so no
#' per-cell or per-group design matrix is ever formed. Results are identical to
#' solving each AF group separately with `unmix.ols.fast()` on the stacked
#' spectra, to floating-point precision.
#'
#' @param raw.data Expression data from raw FCS files. Cells in rows and
#' detectors in columns. Columns must match the columns in `spectra`.
#' @param spectra Spectral signatures of fluorophores, with fluorophores in
#' rows and detectors in columns.
#' @param af.spectra Spectral signatures of autofluorescences, with variants in
#' rows and detectors in columns. Prepare using `get.af.spectra`.
#' @param af.index Integer vector, one entry per cell, giving the row of
#' `af.spectra` assigned to that cell.
#' @param unmixed.no.af Optional numeric matrix (cells x fluorophores) holding
#' the AF-free unmixing `raw.data %*% t(U)`. Supply it when the caller has
#' already computed it to avoid repeating the largest matrix product. Default
#' `NULL`, in which case it is computed here.
#' @param return.fitted.af Logical, default `FALSE`. Whether to also return the
#' fitted autofluorescence in detector space, `k * af.spectra[af.index, ]`.
#' @param denominator.floor Numeric, default `0`. When positive, floors each
#' AF candidate's out-of-span self-dot at this fraction of the largest. An AF
#' variant lying almost inside the fluorophore span has a vanishing out-of-span
#' direction and an unidentifiable abundance; flooring caps the amplification.
#' `0` reproduces the unregularised least-squares solution exactly.
#' @param chunk.size Integer, default `1e5`. Number of cells processed per
#' block, bounding the size of the intermediate detector-wide gather.
#'
#' @return A list with elements `fluorophores` (cells x fluorophores), `af`
#' (numeric vector of per-cell AF abundance) and, when requested,
#' `fitted.af` (cells x detectors).
#'
#' @export

unmix.af.fwl <- function(
    raw.data,
    spectra,
    af.spectra,
    af.index,
    unmixed.no.af     = NULL,
    return.fitted.af  = FALSE,
    denominator.floor = 0,
    chunk.size        = 1e5L
) {

  cell.n <- nrow( raw.data )

  if ( length( af.index ) != cell.n )
    stop( "`af.index` must have one entry per cell.", call. = FALSE )
  if ( any( af.index < 1L ) || any( af.index > nrow( af.spectra ) ) )
    stop( "`af.index` contains rows outside `af.spectra`.", call. = FALSE )

  # pseudoinverse of the fluorophore panel
  unmixing.matrix <- solve.default( tcrossprod( spectra ), spectra )

  # how much each AF variant looks like each fluorophore
  v.library <- unmixing.matrix %*% t( af.spectra )

  # the part of each AF variant that no combination of fluorophores can explain
  r.library <- t( af.spectra ) - ( t( spectra ) %*% v.library )

  denominator <- colSums( r.library^2 )
  if ( denominator.floor > 0 )
    denominator <- pmax( denominator,
                         denominator.floor * max( denominator, 1e-10 ) )
  denominator[ denominator <= 0 ] <- 1e-10

  # AF-free unmixing; the caller may already hold this
  if ( is.null( unmixed.no.af ) )
    unmixed.no.af <- raw.data %*% t( unmixing.matrix )

  r.library.t <- t( r.library )
  v.library.t <- t( v.library )

  af <- numeric( cell.n )

  starts <- seq.int( 1L, cell.n, by = as.integer( chunk.size ) )
  for ( start in starts ) {
    end <- min( start + as.integer( chunk.size ) - 1L, cell.n )
    idx <- start:end
    j   <- af.index[ idx ]

    af[ idx ] <- rowSums( raw.data[ idx, , drop = FALSE ] *
                            r.library.t[ j, , drop = FALSE ] ) / denominator[ j ]
  }

  fluorophores <- unmixed.no.af -
    af * v.library.t[ af.index, , drop = FALSE ]
  colnames( fluorophores ) <- rownames( spectra )

  out <- list( fluorophores = fluorophores, af = af )

  if ( return.fitted.af )
    out$fitted.af <- af * af.spectra[ af.index, , drop = FALSE ]

  out
}
