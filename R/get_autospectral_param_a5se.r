# get_autospectral_param_a5se.r

#' @title Get AutoSpectral Parameters for the FACSymphony A5SE Cytometer
#'
#' @description
#' Returns parameters for running a calculation of unmixing with
#' AutoSpectral, without creating any figures or tables.
#'
#' @param autosp.param A list of initial AutoSpectral parameters.
#'
#' @return A list of AutoSpectral parameters specific to the A5SE cytometer.
#'
#' @export

get.autospectral.param.a5se <- function( autosp.param )
{
  # add cytometer-specific parameters
  autosp.param$cytometer <- "Symphony"

  autosp.param$scatter.data.min.x <- 0

  autosp.param$scatter.data.max.x <- 262144

  autosp.param$scatter.data.min.y <- 0

  autosp.param$scatter.data.max.y <- 262144

  autosp.param$expr.data.min <- -111

  # This should be 262144, but the A5SE appears to record data well above this
  # so we will extend the window to better capture positives in the controls
  autosp.param$expr.data.max <- 2621440

  autosp.param$default.scatter.parameter <- c( "FSC-A", "SSC-A" )

  autosp.param$default.time.parameter <- "Time"

  # for now, keep the default transformation settings
  # with the expected upper limit
  autosp.param$default.transformation.param <- list(
          length = 256,
          max.range = 262144,
          pos = 4.42,
          neg = 0,
          width = -100
        )

  autosp.param$non.spectral.channel <- c( "Time", "SSC", "FSC" )

  autosp.param$af.channel <- "UV515-A"

  autosp.param$data.step <- 5e4

  autosp.param$large.gate.scaling.x <- 2
  autosp.param$large.gate.scaling.y <- 6

  autosp.param$ribbon.breaks <- c( -1e2, 0, 1e2, 1e3, 1e4, 1e5 )

  # lower the similarity threshold since we have fewer detectors
  autosp.param$sim.threshold <- 0.975

  return( autosp.param )

}

