# get_autospectral_param_cytostellar.r

#' @title Get Autospectral Parameters for CytoStellar Cytometer
#'
#' @description
#' Returns parameters for running a calculation of unmixing with
#' AutoSpectral for the CytoStellar, without creating any figures or tables.
#'
#' The CytoStellar allows users to acquire spectral parameters as `-A`
#' (area), `-H` (height), or both. `AutoSpectral` prefers `-A` when both are
#' present; the specific suffix actually present in a given FCS file is
#' resolved at read time by `.resolve.cytostellar.suffix()`, not fixed here.
#' Detector names in `cytometer_database.csv` and `fluorophore_database.csv`
#' are therefore stored WITHOUT a suffix for this cytometer.
#'
#' @param autosp.param A list of initial AutoSpectral parameters.
#'
#' @return A list of AutoSpectral parameters specific to the CytoStellar
#' cytometer.
#'
#' @export

get.autospectral.param.cytostellar <- function( autosp.param )
{
  # add cytometer-specific parameters
  autosp.param$cytometer <- "CytoStellar"

  autosp.param$scatter.data.min.x <- 0
  autosp.param$scatter.data.max.x <- 4000000
  autosp.param$scatter.data.min.y <- 0
  autosp.param$scatter.data.max.y <- 4000000

  autosp.param$expr.data.min <- -4000
  autosp.param$expr.data.max <- 4000000

  autosp.param$default.scatter.parameter <- c( "FSC-A", "BSSC-A" )

  autosp.param$default.transformation.param <- list(
    length = 256,
    max.range = 4000000,
    pos = 5.6,
    neg = 0,
    width = -1000
  )

  # NOTE: deliberately does NOT include "-H"/"-W" here, unlike the Discover
  # family. For CytoStellar those suffixes may be the chosen acquisition
  # parameter, not junk to discard; suffix resolution happens downstream in
  # `.derive.spectral.channels()` / `check.channels()`.
  autosp.param$non.spectral.channel <- c( "Time", "SSC", "FSC" )

  # bare name; suffix resolved and appended at read time
  autosp.param$af.channel <- "V7"

  # preferred suffix when both -A and -H are present
  autosp.param$default.suffix <- "-A"

  autosp.param$data.step <- 5e5

  autosp.param$large.gate.scaling.x <- 1.5
  autosp.param$large.gate.scaling.y <- 4.5

  autosp.param$ribbon.breaks <- c( -1e3, 0, 1e3, 1e4, 1e5, 1e6 )

  return( autosp.param )
}
