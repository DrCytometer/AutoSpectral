# assign_af_residuals.r

#' @title Assign AF Spectrum By Residual Alignment
#'
#' @description
#' Aligns the residuals with the autofluorescence spectral variants to determine
#' which best fits each cell (event). A fast approximation for brute force
#' sequential unmixing method in early versions of AutoSpectral. Provides similar
#' results to residual minimization, and when combined with per-cell optimization,
#' works well. Substantially faster.
#'
#' @param raw.data Expression data from raw FCS files. Cells in rows and
#' detectors in columns. Columns should be fluorescent data only and must
#' match the columns in spectra.
#' @param spectra Spectral signatures of fluorophores, normalized between 0
#' and 1, with fluorophores in rows and detectors in columns.
#' @param af.spectra Spectral signatures of autofluorescences, normalized
#' between 0 and 1, with fluorophores in rows and detectors in columns.
#'
#' @return Row indices for best-fitting AF spectra (from `af.spectra`)
#'
#' @export

assign.af.residuals <- function(
    raw.data,
    spectra,
    af.spectra
) {

  # unmix without autofluorescence
  unmixed <- unmix.ols( raw.data, spectra )

  # calculate initial residual
  residuals <- raw.data - ( unmixed %*% spectra )

  # which AF spectrum best aligns with the residual for each cell?
  scores <- residuals %*% t( af.spectra )

  # maximum score corresponds to the minimum L2 error
  return( max.col( scores ) )
}
