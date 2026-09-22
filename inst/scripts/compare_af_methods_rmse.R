# compare_af_methods_rmse.R

#' @title Compare Autofluorescence Extraction Methods By Reconstruction RMSE
#'
#' @description
#' Unmixes a raw (unstained) sample five ways - no autofluorescence
#' correction, a single population-mean autofluorescence spectrum, and three
#' per-cell autofluorescence assignment methods (`assign.af.residuals()`,
#' `assign.af.fluorophores()`, and the joint covariance-weighted assignment,
#' `assign.af.joint.cov()` or `assign.af.joint.cov.l2()`) - and reports the
#' reconstruction RMSE (`compute.unmix.rmse()`) for each, together with the
#' percent reduction in RMSE relative to the uncorrected (No AF) baseline.
#'
#' For the three per-cell methods, the reconstruction used for RMSE includes
#' the fitted per-cell autofluorescence contribution (`unmix.af.fwl()`'s
#' `fitted.af`), not just the fluorophore term, since the autofluorescence
#' spectrum is being fit as part of the model.
#'
#' @param raw.data Expression data from a raw (unstained) FCS file. Cells in
#' rows and detectors in columns. Columns must match the columns in
#' `spectra` and `af.spectra`.
#' @param spectra Spectral signatures of fluorophores, normalized between 0
#' and 1, with fluorophores in rows and detectors in columns.
#' @param af.spectra Spectral signatures of autofluorescences, normalized
#' between 0 and 1, with AF variants in rows and detectors in columns.
#' Row 1 is used as the single population-mean spectrum for the `SingleAF`
#' method. Prepare using `get.af.spectra()`.
#' @param af.assign.method Character, one of `l1` (default, uses
#' `assign.af.joint.cov()`) or `l2` (uses `assign.af.joint.cov.l2()`),
#' selecting the joint-assignment solver used for the `PerCellAF_Joint`
#' method.
#'
#' @seealso
#' * [compute.unmix.rmse()]
#' * [assign.af.residuals()]
#' * [assign.af.fluorophores()]
#' * [assign.af.joint.cov()]
#' * [assign.af.joint.cov.l2()]
#' * [unmix.af.fwl()]
#'
#' @return A data frame with one row per method (`NoAF`, `SingleAF`,
#' `PerCellAF_Residuals`, `PerCellAF_Fluorophores`, `PerCellAF_Joint`) and
#' columns `Method`, `RMSE`, and `PercentReduction` (percent reduction in
#' RMSE relative to `NoAF`; `0` for the `NoAF` row itself).
#'
#' @export

compare.af.methods.rmse <- function(
    raw.data,
    spectra,
    af.spectra,
    af.assign.method = c( "l1", "l2" )
) {

  af.assign.method <- match.arg( af.assign.method )

  if ( "AF" %in% rownames( spectra ) )
    spectra <- spectra[ rownames( spectra ) != "AF", , drop = FALSE ]

  unmixing.matrix <- solve.default( tcrossprod( spectra ), spectra )
  unmixed.no.af <- raw.data %*% t( unmixing.matrix )

  # ---- No AF (baseline) ----------------------------------------------------

  fitted.no.af <- unmixed.no.af %*% spectra
  rmse.no.af <- compute.unmix.rmse( raw.data, fitted.no.af )

  # ---- Single AF (one population-mean spectrum, plain OLS) -----------------

  single.af.spectrum <- af.spectra[ 1, , drop = FALSE ]
  rownames( single.af.spectrum ) <- "AF"
  single.spectra <- rbind( spectra, single.af.spectrum )

  unmixed.single.af <- unmix.ols( raw.data, single.spectra )
  fitted.single.af <- unmixed.single.af %*% single.spectra
  rmse.single.af <- compute.unmix.rmse( raw.data, fitted.single.af )

  # ---- Per-cell AF, residual alignment --------------------------------------

  af.index.residuals <- assign.af.residuals( raw.data, spectra, af.spectra )
  fit.residuals <- unmix.af.fwl(
    raw.data = raw.data, spectra = spectra, af.spectra = af.spectra,
    af.index = af.index.residuals, unmixed.no.af = unmixed.no.af,
    return.fitted.af = TRUE
  )
  fitted.residuals <- ( fit.residuals$fluorophores %*% spectra ) + fit.residuals$fitted.af
  rmse.residuals <- compute.unmix.rmse( raw.data, fitted.residuals )

  # ---- Per-cell AF, fluorophore projection ----------------------------------

  af.index.fluorophores <- assign.af.fluorophores( raw.data, spectra, af.spectra )
  fit.fluorophores <- unmix.af.fwl(
    raw.data = raw.data, spectra = spectra, af.spectra = af.spectra,
    af.index = af.index.fluorophores, unmixed.no.af = unmixed.no.af,
    return.fitted.af = TRUE
  )
  fitted.fluorophores <- ( fit.fluorophores$fluorophores %*% spectra ) + fit.fluorophores$fitted.af
  rmse.fluorophores <- compute.unmix.rmse( raw.data, fitted.fluorophores )

  # ---- Per-cell AF, joint covariance-weighted assignment --------------------

  af.index.joint <- if ( af.assign.method == "l1" ) {
    assign.af.joint.cov( raw.data, spectra, af.spectra )
  } else {
    assign.af.joint.cov.l2( raw.data, spectra, af.spectra )
  }
  fit.joint <- unmix.af.fwl(
    raw.data = raw.data, spectra = spectra, af.spectra = af.spectra,
    af.index = af.index.joint, unmixed.no.af = unmixed.no.af,
    return.fitted.af = TRUE
  )
  fitted.joint <- ( fit.joint$fluorophores %*% spectra ) + fit.joint$fitted.af
  rmse.joint <- compute.unmix.rmse( raw.data, fitted.joint )

  # ---- assemble ---------------------------------------------------------

  result <- data.frame(
    Method = c(
      "NoAF", "SingleAF", "PerCellAF_Residuals",
      "PerCellAF_Fluorophores", "PerCellAF_Joint"
    ),
    RMSE = c(
      rmse.no.af, rmse.single.af, rmse.residuals, rmse.fluorophores, rmse.joint
    ),
    stringsAsFactors = FALSE
  )

  result$PercentReduction <- 100 * ( rmse.no.af - result$RMSE ) / rmse.no.af

  result
}
