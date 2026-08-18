#' @title Denoised Per-Detector Variability (MAD)
#'
#' @description
#' As `assess.variability()`, but returns the full per-detector robust SD
#' (MAD) profile across a fluorophore's variant spectra, rather than a
#' single summary value. Detectors where the reference spectrum falls
#' below `peak.height` are zeroed out before computing MAD, so that
#' off-peak noise does not contribute to the profile.
#'
#' @param variant.list Named list of spectral variant matrices (one per
#'   fluorophore), as returned in the `variants` element of
#'   `get.spectral.variants()`. Row 1 of each matrix is treated as the
#'   reference spectrum; rows 2+ are treated as variants.
#' @param peak.height Numeric in `[0, 1]`, default `0.1`. Detectors where
#'   the reference spectrum falls below this fraction of that fluorophore's
#'   own peak value (`peak.height * max(reference.spectrum)`) are set to
#'   zero before MAD is computed. Relative to the peak rather than an
#'   absolute cutoff, so behavior is unchanged whether `variant.list` is
#'   L-infinity or L2 normalized.
#'
#' @return A numeric matrix with one row per fluorophore and one column
#'   per detector, giving the denoised per-detector MAD of the variant
#'   spectra.
#'
#' @export

assess.variability.mad <- function(
    variant.list,
    peak.height = 0.1
) {

  fluorophores <- names(variant.list)

  variability <- lapply(
    fluorophores,
    function(fluor) {

      tryCatch(
        {
          # pull best spectrum from first row in the matrix
          fluor.variants <- variant.list[[fluor]]
          reference.spectrum <- fluor.variants[1,]
          variant.spectra <- fluor.variants[2:nrow(fluor.variants),]

          # eliminate off-peak variation (noise) -- relative to this
          # fluorophore's own peak so the threshold is agnostic to whether
          # the spectrum is L-infinity or L2 normalized
          peak.idx <- reference.spectrum > ( peak.height * max( reference.spectrum, na.rm = TRUE ) )
          no.noise.variant.spectra <- variant.spectra
          no.noise.variant.spectra[,!peak.idx] <- 0
          apply(no.noise.variant.spectra, 2, stats::mad)
        },
        error = function(e) {
          message(
            sprintf(
              "[ERROR] Fluorophore %s failed: %s",
              fluor,
              conditionMessage(e)
            )
          )
          return(NULL)
        }
      )

    })

  names(variability) <- fluorophores

  variability <- do.call(rbind, variability)


  return(variability)
}
