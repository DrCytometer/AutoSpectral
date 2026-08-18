#' @title Assess Within-Fluorophore Spectral Variability
#'
#' @description
#' Summarises how much a fluorophore's spectral variants (rows 2+ of its
#' variant matrix, as produced by `get.spectral.variants()`) disagree with
#' its reference spectrum (row 1) and with each other, returning several
#' complementary summary metrics per fluorophore.
#'
#' @param variant.list Named list of spectral variant matrices (one per
#'   fluorophore), as returned in the `variants` element of
#'   `get.spectral.variants()`. Row 1 of each matrix is treated as the
#'   reference spectrum; rows 2+ are treated as variants.
#' @param peak.height Numeric in `[0, 1]`, default `0.1`. Detectors where
#'   the reference spectrum falls below this fraction of that fluorophore's
#'   own peak value (`peak.height * max(reference.spectrum)`) are excluded
#'   from the denoised variability metric (`Denoised_rSD`). Relative to the
#'   peak rather than an absolute cutoff, so behavior is unchanged whether
#'   `variant.list` is L-infinity or L2 normalized.
#'
#' @return A numeric matrix with one row per fluorophore and four columns:
#' \describe{
#'   \item{`rSD`}{Sum, across detectors, of the robust SD (MAD) of the
#'     variant spectra.}
#'   \item{`Scaled`}{`rSD` weighted by the reference spectrum, summed
#'     across detectors, to downweight low-signal channels.}
#'   \item{`Sim_rSD`}{Robust SD of the pairwise cosine similarities among
#'     all variants (including the reference).}
#'   \item{`Denoised_rSD`}{As `rSD`, but restricted to detectors where the
#'     reference spectrum exceeds `peak.height`.}
#' }
#' Returns `NULL` (with a message) if no fluorophore has usable
#' variability data.
#'
#' @export

assess.variability <- function(
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

          # Ensure we have a matrix with at least 2 rows
          if (is.null(fluor.variants) || nrow(fluor.variants) < 2) {
            stop("Insufficient data (less than 2 rows)")
          }

          reference.spectrum <- fluor.variants[1,]
          variant.spectra <- fluor.variants[2:nrow(fluor.variants),]

          # get robust SD of variant spectra across channels
          rsd <- apply(variant.spectra, 2, stats::mad)

          # multiply by best spectrum to downscale background noise
          variation <- rsd * reference.spectrum

          # eliminate off-peak variation (noise) -- relative to this
          # fluorophore's own peak so the threshold is agnostic to whether
          # the spectrum is L-infinity or L2 normalized
          peak.idx <- reference.spectrum > ( peak.height * max( reference.spectrum, na.rm = TRUE ) )

          # Handle cases where no peaks are found
          if (!any(peak.idx)) {
            no.noise.rsd <- rep(0, length(rsd))
          } else {
            no.noise.variant.spectra <- variant.spectra
            no.noise.variant.spectra[, !peak.idx] <- 0
            no.noise.rsd <- apply(no.noise.variant.spectra, 2, stats::mad)
          }

          # get rSD of distribution of cosine similarity values
          sim.matrix <- cosine.similarity(fluor.variants)
          sim.values <- sim.matrix[lower.tri(sim.matrix)]
          sim.rsd <- if(length(sim.values) > 0) stats::mad(sim.values) else 0

          return(c(
            sum(rsd, na.rm = TRUE),
            sum(abs(variation), na.rm = TRUE),
            sim.rsd,
            sum(no.noise.rsd, na.rm = TRUE)
          ))

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

  # Remove NULLs before rbinding
  variability <- variability[!sapply(variability, is.null)]

  if (length(variability) == 0) {
    message("No valid variability data calculated.")
    return(NULL)
  }

  variability <- do.call(rbind, variability)

  colnames(variability) <- c("rSD", "Scaled", "Sim_rSD", "Denoised_rSD")

  return(variability)
}
