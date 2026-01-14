# get_ungated_flow_expression_data.r

#' @title Get Ungated Flow Expression Data
#'
#' @description
#' Retrieves flow cytometry expression data for specified samples, without
#' gating, removing out-of-range events.
#'
#' @importFrom flowCore read.FCS exprs
#'
#' @param samp The sample identifier.
#' @param file.name A vector of file names for the samples.
#' @param control.dir The directory containing the control files.
#' @param scatter.and.spectral.channel A vector of scatter and spectral channels.
#' @param spectral.channel A vector of spectral channels.
#' @param set.resolution The resolution limit for the spectral channels.
#'
#' @return A matrix with the flow expression data.

get.ungated.flow.expression.data <- function( samp, file.name, control.dir,
                                              scatter.and.spectral.channel,
                                              spectral.channel, set.resolution ) {

  flow.file <- file.name[ samp ]

  fcs.data <- suppressWarnings(
    flowCore::read.FCS(
      file.path( control.dir, flow.file ),
      transformation = NULL,
      truncate_max_range = FALSE,
      emptyValue = FALSE
    )
  )

  # read exprs for spectral channels only
  expr.data <- flowCore::exprs( fcs.data )[ , scatter.and.spectral.channel ]

  rm( fcs.data )

  # remove any out-of-range events
  below.resolution.limit <- apply(
    expr.data[ , spectral.channel ], 1, function( flow.event ) {
      all( flow.event < set.resolution ) }
    )

  expr.data <- expr.data[ below.resolution.limit, ]

  return( expr.data )

}

