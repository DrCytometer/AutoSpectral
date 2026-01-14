# sample_fcs_file.r

#' @title Sample FCS File
#'
#' @description
#' This function samples events from an FCS file based on specified parameters
#' and downsampling criteria.
#'
#' @importFrom flowCore read.FCS exprs
#'
#' @param file.name A character string specifying the name of the FCS file.
#' @param control.dir A character string specifying the directory containing
#' the control FCS file.
#' @param downsample.n A numeric value indicating the number of events to
#' downsample.
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#'
#' @return A matrix containing the sampled events from the FCS file.

sample.fcs.file <- function( file.name, control.dir, downsample.n, asp ) {

  ff <- suppressWarnings(
    flowCore::read.FCS(
      file.path( control.dir, file.name ),
      transformation = FALSE,
      truncate_max_range = FALSE,
      emptyValue = FALSE
    )
  )

  ff <- flowCore::exprs( ff )[ , asp$default.scatter.parameter ]

  event.n <- nrow( ff )

  if ( event.n <= asp$min.cell.warning.n ) {
    stop(
      paste0(
        "Fewer than0",
        asp$min.cell.warning.n,
        " events in ",
        file.name
      ),
      call. = FALSE
    )
  }

  downsample.n <- if ( event.n < downsample.n ) event.n else downsample.n

  fcs.idx <- sample( 1:event.n, downsample.n )

  ff <- ff[ fcs.idx, ]

  return( ff )
}
