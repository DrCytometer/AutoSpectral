#' @title Unmix WLS Per Cell
#'
#' @description
#' Faster solver for per-cell optimization workflow. Performs spectral unmixing
#' using weighted least squares.
#'
#' @param cell.raw Expression data from raw FCS files. A single row of cells with
#' detectors in columns. Columns should be fluorescent data only and must
#' match the columns in spectra.
#' @param spectra Spectral signatures of fluorophores, normalized between 0
#' and 1, with fluorophores in rows and detectors in columns.
#' @param weights The weighting values for the detectors. No checks are performed
#' in this function; `unmix.wls` should be used for standard cases.
#'
#' @return Unmixed data (one row), fluorophores in columns.
#'
#' @export

unmix.wls.per.cell <- function(
    cell.raw,
    spectra,
    weights
) {
  # apply weights directly
  Sw <- spectra * rep( weights, each = nrow( spectra ) )

  # use crossprod to compute (S W S')
  XtX <- tcrossprod( Sw, spectra )

  # compute (S W y)
  Xty <- Sw %*% t( cell.raw )

  # solve the system Ax = b directly
  return( t( solve( XtX, Xty ) ) )
}
