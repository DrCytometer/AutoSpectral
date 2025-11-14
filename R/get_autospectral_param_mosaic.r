# get_autospectral_param_mosaic.r

#' @title Get Autospectral Parameters for Mosaic Cytometer
#'
#' @description
#' Returns parameters for running a calculation of unmixing with
#' AutoSpectral, without creating any figures or tables.
#'
#' @param autosp.param A list of initial AutoSpectral parameters.
#'
#' @return A list of AutoSpectral parameters specific to the Mosaic cytometer.
#'
#' @export

get.autospectral.param.mosaic <- function( autosp.param )
{
  # add cytometer-specific parameters
  autosp.param$cytometer <- "Mosaic"

  autosp.param$scatter.data.min.x <- 0

  autosp.param$scatter.data.max.x <- 4000000

  autosp.param$scatter.data.min.y <- 0

  autosp.param$scatter.data.max.y <- 4000000

  autosp.param$expr.data.min <- -111

  autosp.param$expr.data.max <- 16777215

  autosp.param$default.scatter.parameter <- c( "FSC-A", "BSSC-A" )

  autosp.param$default.time.parameter <- "Time"

  autosp.param$default.transformation.param <- list(
          length = 256,
          max.range = 16777215,
          pos = 6.22,
          neg = 0,
          width = -500
        )

  autosp.param$non.spectral.channel <- c( "Time", "SSC", "FSC", "-H" )

  autosp.param$af.channel <- "V8-A"

  autosp.param$data.step <- 5e5

  autosp.param$large.gate.scaling.x <- 1.5
  autosp.param$large.gate.scaling.y <- 4.5

  autosp.param$ribbon.breaks <- c( -1e3, 0, 1e3, 1e4, 1e5, 1e6 )

  autosp.param

}

