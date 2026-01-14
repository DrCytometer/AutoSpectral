# run_trim_events.r

#' @title Run Trim Events
#'
#' @description
#' This function trims extreme events from multiple samples based on specified
#' peak channels and trim factors.
#'
#' @param trim.sample.data A list containing the expression data for each
#' sample.
#' @param trim.sample A character vector specifying the names of the samples.
#' @param trim.peak.channels A list containing the peak channels for each
#' sample.
#' @param trim.factor A numeric value indicating the proportion of extreme
#' events to trim from both ends of the peak channel distribution.
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#'
#' @return A list containing the trimmed expression data for each sample.

run.trim.events <- function( trim.sample.data, trim.sample,
                            trim.peak.channels, trim.factor, asp ){

  trimmed.expr <- lapply( names( trim.sample.data ), function( sample.name ){
      trim.extreme.events( trim.sample.data[[ sample.name ]],
                           trim.peak.channels[[ sample.name ]], trim.factor )
    } )

  names( trimmed.expr ) <- names( trim.sample.data )

  # issue warning if fewer than 500 events for any sample
  trimmed.expr.n <- sapply( trimmed.expr, nrow )

  low.sample.n <- which( trimmed.expr.n < asp$min.cell.warning.n )

  if ( any( trimmed.expr.n < asp$min.cell.warning.n ) ) {
    warning(
      paste0(
        "\033[31m",
        "Warning! Fewer than ",
        asp$min.cell.warning.n,
        " gated events in ",
        names( low.sample.n ),
        "\033[0m",
        "\n"
      )
    )

  }

  low.sample.n <- which( trimmed.expr.n < asp$min.cell.stop.n )

  # for any samples that have been trimmed to fewer than 50 events,
  # return untrimmed data
  if ( any( trimmed.expr.n < asp$min.cell.stop.n ) ){
    for ( low.n in names( low.sample.n ) ) {
      trimmed.expr[[ low.n ]] <- trim.sample.data[[ low.n ]]
    }
  }

  rm( trim.sample.data )

  if ( any( trimmed.expr.n < asp$min.cell.stop.n ) ) {
    warning(
      paste0(
        "\033[31m",
        "Warning! Fewer than ",
        asp$min.cell.stop.n,
        " gated events in ",
        names( low.sample.n ),
        "\033[0m",
        "\n"
      )
    )
  }

  return( trimmed.expr )
}
