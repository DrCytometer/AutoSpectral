# get_autospectral_param_id7000.r

#' @title Get Autospectral Parameters for ID7000 Cytometer
#'
#' @description
#' Returns parameters for running a calculation of unmixing with
#' AutoSpectral for the ID7000, without creating any figures or tables.
#'
#' @param autosp.param A list of initial AutoSpectral parameters.
#'
#' @return A list of AutoSpectral parameters specific to the ID7000 cytometer.
#'
#' @export

get.autospectral.param.id7000 <- function( autosp.param )
{
  # add cytometer-specific parameters
  autosp.param$cytometer <- "ID7000"

  autosp.param$scatter.data.min.x <- 0

  autosp.param$scatter.data.max.x <- 1048576

  autosp.param$scatter.data.min.y <- 0

  autosp.param$scatter.data.max.y <- 1048576

  autosp.param$expr.data.min <- -111

  autosp.param$expr.data.max <- 1000000

  autosp.param$default.scatter.parameter <- c( "FSC-A", "SSC-A" )

  autosp.param$default.time.parameter <- "TIME"

  autosp.param$default.transformation.param <- list(
    length = 256,
    max.range = 1000000,
    pos = 5,
    neg = 0,
    width = -250
  )

  autosp.param$non.spectral.channel <- c( "TIME", "SSC", "FSC" )

  autosp.param$af.channel <- "405CH7-A"

  autosp.param$data.step <- 1e5

  autosp.param$large.gate.scaling.x <- 3
  autosp.param$large.gate.scaling.y <- 6

  autosp.param$ribbon.breaks <- c( -1e3, 0, 1e3, 1e4, 1e5, 1e6 )

  return( autosp.param )

}
