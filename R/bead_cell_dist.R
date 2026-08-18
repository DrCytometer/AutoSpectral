#' @title Per-Detector Spectral Distance Between Particle Types
#'
#' @description
#' For each fluorophore present in `cell.variants`, computes the raw
#' per-detector difference (cell reference spectrum minus bead reference
#' spectrum) between two particle types' extracted reference spectra.
#'
#' `cell.variants` and `bead.variants` are expected to come from
#' independent `get.spectral.variants()` calls (e.g. one per
#' particle-type folder) and are not guaranteed to share the same
#' detector set or column order, so the comparison is aligned explicitly
#' by detector name rather than position. A warning is issued if the two
#' particle types were extracted with different detector sets.
#'
#' @param cell.variants Named list of spectral variant matrices (one per
#'   fluorophore), as returned in the `variants` element of
#'   `get.spectral.variants()`, for the reference particle type
#'   (typically Cells). Row 1 of each matrix is treated as the reference
#'   spectrum.
#' @param bead.variants Named list of spectral variant matrices in the
#'   same format as `cell.variants`, for the particle type being compared
#'   against the reference (e.g. a bead type).
#'
#' @return A numeric matrix of per-detector differences (cell minus
#'   bead), one row per fluorophore common to both lists, columns are the
#'   shared detector channels.
#'
#' @export

bead.cell.dist <- function(
    cell.variants,
    bead.variants
) {

  fluorophores <- intersect(names(cell.variants), names(bead.variants))

  if (length(fluorophores) == 0) {
    stop(
      "No fluorophore names in common between `cell.variants` and ",
      "`bead.variants` -- check that both particle types resolved the ",
      "same fluorophore set.",
      call. = FALSE
    )
  }

  # find the first fluorophore with a properly dimnamed matrix on both
  # sides to use as the reference channel set; report anything malformed
  # along the way rather than assuming fluorophores[1] is representative
  malformed <- character(0)
  ref.fluor <- NA_character_

  for (fl in fluorophores) {
    cell.mat <- cell.variants[[fl]]
    bead.mat <- bead.variants[[fl]]
    cell.ok <- !is.null(cell.mat) && !is.null(colnames(cell.mat))
    bead.ok <- !is.null(bead.mat) && !is.null(colnames(bead.mat))
    if (cell.ok && bead.ok) {
      ref.fluor <- fl
      break
    }
    if (!cell.ok) malformed <- c(malformed, paste0(fl, " (cell)"))
    if (!bead.ok) malformed <- c(malformed, paste0(fl, " (bead)"))
  }

  if (is.na(ref.fluor)) {
    stop(
      "None of the shared fluorophore(s) have a properly dimnamed variant ",
      "matrix in both `cell.variants` and `bead.variants` -- check that ",
      "both particle types were processed with the same asp$cytometer, ",
      "and that cached extraction .rds files are not stale. Malformed ",
      "entries: ", paste(malformed, collapse = ", "), ".",
      call. = FALSE
    )
  }

  cell.channels <- colnames(cell.variants[[ref.fluor]])
  bead.channels <- colnames(bead.variants[[ref.fluor]])
  common.channels <- intersect(cell.channels, bead.channels)

  if (length(common.channels) == 0) {
    stop(
      "No shared detector channels between cell and bead spectra -- check ",
      "that both particle types were processed with the same asp$cytometer.",
      call. = FALSE
    )
  }

  if (!setequal(cell.channels, bead.channels)) {
    warning(
      sprintf(
        paste(
          "Cell and bead detector sets differ (%d shared of %d cell / %d",
          "bead channels). Comparing on the shared channel(s) only.",
          "Missing from bead: %s. Missing from cell: %s."
        ),
        length(common.channels), length(cell.channels), length(bead.channels),
        paste(setdiff(cell.channels, bead.channels), collapse = ", "),
        paste(setdiff(bead.channels, cell.channels), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (length(malformed) > 0)
    warning(
      "Fluorophore(s) with a missing/undimnamed variant matrix will fail ",
      "individually below: ", paste(malformed, collapse = ", "), ".",
      call. = FALSE
    )

  mismatch <- lapply(
    fluorophores,
    function(fluor) {

      tryCatch(
        {
          # pull best spectrum from first row in each matrix, aligned by
          # detector name rather than position
          cell.spectrum <- cell.variants[[fluor]][1, common.channels]
          bead.spectrum <- bead.variants[[fluor]][1, common.channels]

          cell.spectrum - bead.spectrum
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

  names(mismatch) <- fluorophores

  mismatch <- do.call(rbind, mismatch)

  return(mismatch)
}
