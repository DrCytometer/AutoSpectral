#' @title Assess Spectral Mismatch Between Reference and Test Spectra
#'
#' @description
#' For each fluorophore common to both `reference.variants` and
#' `test.variants`, computes the cosine similarity between the reference
#' spectrum (row 1 of each fluorophore's variant matrix) and the
#' corresponding test spectrum. Typically used to compare a spectral
#' reference library built from one particle type (e.g. Cells) against one
#' built from another (e.g. a bead type), by passing the reference particle
#' type as `reference.variants`.
#'
#' `reference.variants` and `test.variants` are expected to come from
#' independent `get.spectral.variants()` calls and are not guaranteed to
#' share the same detector set or column order, so the comparison is
#' aligned explicitly by detector name rather than position. A warning is
#' issued if the two detector sets differ.
#'
#' @param reference.variants Named list of spectral variant matrices (one
#'   per fluorophore), as returned in the `variants` element of
#'   `get.spectral.variants()`. Row 1 of each matrix is treated as the
#'   reference spectrum.
#' @param test.variants Named list of spectral variant matrices in the same
#'   format as `reference.variants`, to be compared against it.
#'
#' @return A one-column numeric matrix of cosine similarity values, one row
#'   per fluorophore common to both lists.
#'
#' @export

assess.mismatch <- function(
    reference.variants,
    test.variants
) {
  fluorophores <- intersect(names(reference.variants), names(test.variants))
  
  ref.channels <- colnames(reference.variants[[fluorophores[1]]])
  test.channels <- colnames(test.variants[[fluorophores[1]]])
  common.channels <- intersect(ref.channels, test.channels)
  
  if (length(common.channels) == 0) {
    stop(
      "No shared detector channels between reference and test spectra.",
      call. = FALSE
    )
  }
  
  if (!setequal(ref.channels, test.channels)) {
    warning(
      sprintf(
        paste(
          "Reference and test detector sets differ (%d shared of %d",
          "reference / %d test channels). Comparing on the shared",
          "channel(s) only. Missing from test: %s. Missing from reference: %s."
        ),
        length(common.channels), length(ref.channels), length(test.channels),
        paste(setdiff(ref.channels, test.channels), collapse = ", "),
        paste(setdiff(test.channels, ref.channels), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  mismatches <- lapply(
    fluorophores,
    function(fluor) {
      tryCatch({
        # find corresponding best spectrum in each matrix, aligned by
        # detector name rather than position
        reference <- reference.variants[[fluor]][1, common.channels]
        test <- test.variants[[fluor]][1, common.channels]
        
        # assess similarity
        sim <- cosine.similarity(rbind(reference, test))
        sim[lower.tri(sim)]
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
      })
    })
  
  names(mismatches) <- fluorophores
  
  return(do.call(rbind, mismatches))
}