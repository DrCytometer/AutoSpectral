# trim_extreme_events.r

#' @title Trim Extreme Events
#'
#' @description
#' This function trims extreme events from expression data based on a specified
#' peak channel and trim factor. Not recommended for most panels.
#'
#' @param expr.data A  matrix containing the expression data.
#' @param peak.channel A character string specifying the peak channel to be
#' used for trimming.
#' @param trim.factor A numeric value indicating the proportion of extreme
#' events to trim from both ends of the peak channel distribution.
#'
#' @return A matrix with the extreme events trimmed.

trim.extreme.events <- function( expr.data, peak.channel, trim.factor ){

  # get expr.data for peak channel
  peak.channel.expr <- expr.data[ , peak.channel ]

  peak.channel.expr.n <- length( peak.channel.expr )

  expr.trim.n <- round( peak.channel.expr.n * trim.factor ) + 1

  peak.channel.expr.low <- sort( peak.channel.expr )[ expr.trim.n ]
  peak.channel.expr.high <- sort(
    peak.channel.expr, decreasing = TRUE )[ expr.trim.n ]

  expr.trim.idx <- which (
    peak.channel.expr > peak.channel.expr.low &
      peak.channel.expr < peak.channel.expr.high
    )

  return( expr.data[ expr.trim.idx, ] )

}
