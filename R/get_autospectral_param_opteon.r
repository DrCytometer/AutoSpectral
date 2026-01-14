# get_autospectral_param_opteon.r

#' @title Get Autospectral Parameters for Opteon Cytometer
#'
#' @description
#' Returns parameters for running a calculation of unmixing with
#' AutoSpectral for the Opteon, without creating any figures or tables.
#'
#' @param autosp.param A list of initial AutoSpectral parameters.
#'
#' @return A list of AutoSpectral parameters specific to the Opteon cytometer.
#'
#' @export

get.autospectral.param.opteon <- function( autosp.param )
{
  # add cytometer-specific parameters
  autosp.param$cytometer <- "Opteon"

  autosp.param$scatter.data.min.x <- 0

  autosp.param$scatter.data.max.x <- 13841000

  autosp.param$scatter.data.min.y <- 0

  autosp.param$scatter.data.max.y <- 8687233

  autosp.param$expr.data.min <- -111

  autosp.param$expr.data.max <- 16777216

  autosp.param$default.scatter.parameter <- c( "FSC-A", "VSSC-A" )

  autosp.param$default.time.parameter <- "Time"

  autosp.param$default.transformation.param <- list(
          length = 256,
          max.range = 1e7,
          pos = 6.68,
          neg = 0,
          width = -1000
        )

  autosp.param$non.spectral.channel <- c( "Time", "SSC", "FSC", "-H", "Width" )

  autosp.param$af.channel <- "UV508-A"

  autosp.param$data.step <- 3e6

  # spectral parameters

  autosp.param$ribbon.breaks <- c( -1e3, 0, 1e3, 1e4, 1e5, 1e6, 1e7 )

  autosp.param

}

