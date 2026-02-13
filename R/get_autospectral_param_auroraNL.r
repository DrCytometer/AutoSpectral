# get_autospectral_param_auroraNL.r

#' @title Get Autospectral Parameters for the Aurora Northern Lights Cytometer
#'
#' @description
#' Returns parameters for running a calculation of unmixing with
#' AutoSpectral, without creating any figures or tables.
#'
#' @param autosp.param A list of initial AutoSpectral parameters.
#'
#' @return A list of AutoSpectral parameters specific to the Aurora cytometer.
#'
#' @export

get.autospectral.param.auroraNL <- function( autosp.param )
{
  # add cytometer-specific parameters
  autosp.param$cytometer <- "Aurora"
  autosp.param$cytometer.version <- "NL"

  autosp.param$scatter.data.min.x <- 0

  autosp.param$scatter.data.max.x <- 4194304

  autosp.param$scatter.data.min.y <- 0

  autosp.param$scatter.data.max.y <- 4194304

  autosp.param$expr.data.min <- -111

  autosp.param$expr.data.max <- 4194304

  autosp.param$default.scatter.parameter <- c( "FSC-A", "SSC-A" )

  autosp.param$default.time.parameter <- "Time"

  autosp.param$default.transformation.param <- list(
          length = 256,
          max.range = 4194304,
          pos = 5.62,
          neg = 0,
          width = -1000
        )

  autosp.param$non.spectral.channel <- c( "Time", "SSC", "FSC" )

  autosp.param$af.channel <- "V7-A"

  autosp.param$data.step <- 5e5

  autosp.param$large.gate.scaling.x <- 1.5
  autosp.param$large.gate.scaling.y <- 4.5

  autosp.param$ribbon.breaks <- c( -1e3, 0, 1e3, 1e4, 1e5, 1e6 )

  # lower the similarity threshold since we have fewer detectors
  autosp.param$sim.threshold <- 0.975

  return( autosp.param )

}

