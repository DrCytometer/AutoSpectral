# calculate_weights.r

#' @title Calculate Weights
#'
#' @description
#' This function calculates weights for generating a weighted least-squares (WLS)
#' unmixing matrix. It does so using mean expression levels in the provided FCS
#' file. The provided FCS file should thus be representative of the data to which
#' the unmixing matrix will be applied. One option would be to concatenate samples
#' of all FCS files to be unmixed and use that as input to `fcs.file`.
#'
#' @importFrom flowCore read.FCS exprs
#'
#' @param fcs.file A character string specifying the path to the FCS file.
#' @param spectral.channels Character string specifying the names of the spectral
#' detectors in the FCS file. Can be obtained from `flow.control$spectral.channel`
#' or `colnames(spectra)`.
#' @param save Logical, if `TRUE`, save the weights to a CSV file in `output.dir`.
#' @param output.dir Character string specifying the directory to save the CSV
#' file, if generated.
#' @param filename Character string specifying the filename for the CSV file.
#' Default is `weights.csv`.
#'
#' @return A named numeric vector with weighting for each detector.
#'
#' @export


calculate.weights <- function( fcs.file, spectral.channels,
                               save = FALSE,
                               output.dir = "./table_spectra",
                               filename = "weights.csv" ) {

  if ( save & !dir.exists( output.dir ) )
    dir.create( output.dir )

  # import fcs, without warnings for fcs 3.2
  spectral.exprs <- suppressWarnings(
    flowCore::exprs(
      flowCore::read.FCS(
        fcs.file,
        transformation = FALSE,
        truncate_max_range = FALSE,
        emptyValue = FALSE
        )
    )[ , spectral.channels ]
  )

  weights <- pmax( abs( colMeans( spectral.exprs ) ), 1e-6 )
  weights <- 1 / weights

  names( weights ) <- spectral.channels

  if ( save )
    write.csv( weights, file = file.path( output.dir, filename ) )

  return( weights )
}
