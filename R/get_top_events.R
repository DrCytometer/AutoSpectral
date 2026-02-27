# get_top_events.r

#' @title Get Top (Brightest) Events
#'
#' @description
#' Retrieves the brightest `n` events in the specified FCS file in the specified
#' `channel`.
#'
#' @param fcs.path File path and name for the FCS file to be read.
#' @param channel Channel (detector) name for the selection of the brightest events.
#' @param n Number of brightest events to retrieve, default is `2000`.
#' @param scatter.param The names of the scatter channels to define which data
#' are returned.
#'
#' @return Returns a matrix of scatter (FSC, SSC) data
#'
#' @export
#'

get.top.events <- function( fcs.path, channel, n = 2000, scatter.param ) {
  mat <- readFCS( fcs.path, return.keywords = FALSE )

  # Check if channel exists
  if ( !channel %in% colnames( mat ) ) return( NULL )
  if ( !all( scatter.param %in% colnames( mat ) ) ) return( NULL )

  # Order by the specific channel descending
  ord <- order( mat[ , channel ], decreasing = TRUE )

  # Return top N rows
  return( mat[ utils::head( ord, n ), scatter.param ] )
}
