# get_gated_flow_expression_data.r

#' @title Get Gated Flow Expression Data
#'
#' @description
#' Retrieves gated flow cytometry expression data for specified
#' samples, removing out-of-range events and applying gating boundaries.
#'
#' @importFrom flowCore read.FCS exprs
#' @importFrom sp point.in.polygon
#'
#' @param samp The sample identifier.
#' @param file.name A vector of file names for the samples.
#' @param control.dir The directory containing the control files.
#' @param scatter.and.spectral.channel A vector of scatter and spectral channels.
#' @param spectral.channel A vector of spectral channels.
#' @param set.resolution The resolution limit for the spectral channels.
#' @param flow.gate A list of flow gates for the samples.
#' @param gate.list A list of gating boundaries.
#' @param scatter.param A vector of scatter parameters.
#' @param scatter.and.channel.label A label for scatter and channel.
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#'
#' @return A matrix with the gated expression data.

get.gated.flow.expression.data <- function( samp, file.name, control.dir,
                                            scatter.and.spectral.channel,
                                            spectral.channel, set.resolution,
                                            flow.gate, gate.list, scatter.param,
                                            scatter.and.channel.label, asp
                                            ) {

  flow.file <- file.name[ samp ]

  fcs.data <- suppressWarnings(
    flowCore::read.FCS(
      file.path( control.dir, flow.file ),
      transformation = NULL,
      truncate_max_range = FALSE,
      emptyValue = FALSE
      )
    )

  # read exprs for scatter and spectral channels only
  expr.data <- flowCore::exprs( fcs.data )[ , scatter.and.spectral.channel ]

  rm( fcs.data )

  # remove any out-of-range events
  below.resolution.limit <- apply(
    expr.data[ , spectral.channel ], 1, function( flow.event ) {
    all( flow.event < set.resolution ) }
    )

  expr.data <- expr.data[ below.resolution.limit, ]

  # gate
  gate.idx <- flow.gate[[ samp ]]

  gate.population.boundary <- gate.list[[ gate.idx ]]

  gate.data <- expr.data[ , scatter.param ]

  gate.population.pip <- sp::point.in.polygon(
    gate.data[ , 1 ], gate.data[ , 2 ],
    gate.population.boundary$x, gate.population.boundary$y
    )

  gate.population.idx <- which( gate.population.pip != 0 )

  # plot gate applied to sample
  if ( ! is.null( asp$figure.gate.dir ) ) {
    message( paste0( "\033[34m", "Plotting gate for: ", samp, "\033[0m" ) )
    suppressWarnings(
      gate.sample.plot(
        samp,
        gate.data,
        scatter.param,
        gate.population.boundary,
        scatter.and.channel.label,
        "cells",
        asp
      )
    )
  }

  return( expr.data[ gate.population.idx, ] )

}



