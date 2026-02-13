# get_autospectral_param_xenith.r

#' @title Get Autospectral Parameters for Xenith Cytometer
#'
#' @description
#' Returns parameters for running a calculation of unmixing with
#' AutoSpectral, without creating any figures or tables.
#'
#' @param autosp.param A list of initial AutoSpectral parameters.
#'
#' @return A list of AutoSpectral parameters specific to the Xenith cytometer.
#'
#' @export

get.autospectral.param.xenith <- function( autosp.param )
{
  # add cytometer-specific parameters
  autosp.param$cytometer <- "Xenith"

  autosp.param$scatter.data.min.x <- 0

  autosp.param$scatter.data.max.x <- 100000

  autosp.param$scatter.data.min.y <- 0

  autosp.param$scatter.data.max.y <- 100000

  autosp.param$expr.data.min <- -20

  autosp.param$expr.data.max <- 100000

  autosp.param$default.scatter.parameter <- c( "FSC51-A", "SSC52-A" )

  autosp.param$default.time.parameter <- "Time"

  autosp.param$default.transformation.param <- list(
          length = 256,
          max.range = 100000,
          pos = 5.00,
          neg = 0,
          width = -20
        )

  autosp.param$non.spectral.channel <- c(
    "Time", "SSC", "FSC", "-H", "-W", "Comp", "Event", "Gate", "Sort"
  )

  autosp.param$af.channel <- "FL13-A"

  autosp.param$data.step <- 1e4

  autosp.param$large.gate.scaling.x <- 1.5
  autosp.param$large.gate.scaling.y <- 4.5

  autosp.param$ribbon.breaks <- c( -1e3, 0, 1e3, 3e3, 1e4, 3e4 )

  # lower the similarity threshold since we have fewer detectors
  autosp.param$sim.threshold <- 0.975

  return( autosp.param )

}

