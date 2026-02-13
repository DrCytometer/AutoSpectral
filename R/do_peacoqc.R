# do_peacoqc.r

#' @title Perform Quality Control on Flow Cytometry Data using PeacoQC
#' @description
#' This function performs quality control on flow cytometry data using the
#' PeacoQC method. It transforms the data, removes margins, and identifies
#' good cells. The function also optionally plots the results and saves
#' the cleaned data.
#'
#' @importFrom flowCore transform flowFrame exprs keyword transformList keyword<-
#' @importFrom PeacoQC PeacoQC PlotPeacoQC RemoveMargins
#'
#' @param dirty.expr A matrix containing the raw expression data, pre-cleaning.
#' @param sample.name The name of the sample.

#' @param spectral.channel The spectral channels to be used.
#' @param biexp.transform The biexponential transformation function.
#' @param transform.inv The inverse transformation function.
#' @param output.dir The directory to save output files.
#' @param time.param The time channel parameter.
#' @param all.channels A vector of all channels to be included in the final
#' cleaned data.
#' @param method The PeacoQC method to use. Inherited from `PeacoQC`. Options
#' are `all`, `MAD` or `IT`. Default is `MAD`.
#' @param figures Logical. Controls whether PeacoQC plots are created. Default
#' is `TRUE`.
#' @param verbose Logical, default is `TRUE`. Set to `FALSE` to suppress messages.
#'
#' @return A matrix with the cleaned expression data.

do.peacoQC <- function(
    dirty.expr,
    sample.name,
    spectral.channel,
    biexp.transform,
    transform.inv,
    output.dir,
    time.param,
    all.channels,
    method = "MAD",
    figures = TRUE,
    verbose = TRUE
  ) {

  if ( !dir.exists( output.dir ) & figures )
    dir.create( output.dir )

  if ( verbose )
    message( paste0( "\033[34m", "Performing time-based cleaning of ", sample.name,
                     " using PeacoQC", "\033[0m" ) )

  transform.list <- flowCore::transformList( spectral.channel, biexp.transform )

  dirty.ff <- flowCore::flowFrame( dirty.expr )

  flowCore::keyword( dirty.ff )$FILENAME <- sample.name

  dirty.ff <- flowCore::transform( dirty.ff, transform.list )

  dirty.ff <- PeacoQC::RemoveMargins( dirty.ff, spectral.channel )

  peacoQC.result <- suppressWarnings(
    PeacoQC::PeacoQC(
      ff = dirty.ff,
      channels = spectral.channel,
      determine_good_cells = method,
      plot = FALSE, save_fcs = FALSE,
      output_directory = output.dir,
      report = FALSE, time_channel_parameter = time.param
    ) )

  if ( figures )
    PeacoQC::PlotPeacoQC(
      dirty.ff, spectral.channel, output.dir,
      display_peaks = peacoQC.result
      )

  transform.list <- flowCore::transformList( spectral.channel, transform.inv )

  peacoQC.result$FinalFF <- flowCore::transform(
    peacoQC.result$FinalFF, transform.list )

  clean.expr <- flowCore::exprs( peacoQC.result$FinalFF )[ , all.channels ]

  return( clean.expr )

}
