# downsample_control.r


#' @title Downsample Control Data
#' @description
#' Downsamples control data by selecting a specified number of
#' positive and negative events based on peak channel values.
#'
#' @param clean.expr.data A list containing cleaned expression data for each
#' sample.
#' @param samp The sample identifier.
#' @param peak.channels A vector of peak channels for the samples.
#' @param negative.n Number of negative events to select, default `500`.
#' @param positive.n Number of positive events to select, default `1000`.
#' @param verbose Logical. Default is `TRUE`.
#'
#' @return A matrix with the selected expression data.

downsample.control <- function( clean.expr.data, samp, peak.channels,
                                negative.n = 500, positive.n = 1000,
                                verbose = TRUE ) {

  if ( verbose ) message( paste0( "\033[34m", "Downsampling ", samp, "\033[0m" ) )

  # should just pass control's data
  pos.control.expr <- clean.expr.data[[ samp ]]
  # ditto
  peak.channel <- peak.channels[ samp ]

  pos.peak.channel <- pos.control.expr[ , peak.channel ]

  pos.selected <- sort( pos.peak.channel, decreasing = TRUE )[ 1:positive.n ]
  neg.selected <- sort( pos.peak.channel, decreasing = FALSE )[ 1:negative.n ]

  # recover full data
  selected.expr <- pos.control.expr[ names( c( pos.selected, neg.selected ) ), ]

  return( selected.expr )

}
