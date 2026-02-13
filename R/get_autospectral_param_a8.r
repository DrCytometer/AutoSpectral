# get_autospectral_param_a8.r

#' @title Get Autospectral Parameters for DiscoverA8 Cytometer
#'
#' @description
#' Returns parameters for running a calculation of unmixing with
#' AutoSpectral, without creating any figures or tables.
#'
#' @param autosp.param A list of initial AutoSpectral parameters.
#'
#' @return A list of AutoSpectral parameters specific to the DiscoverA8
#' cytometer.
#'
#' @export

get.autospectral.param.a8 <- function( autosp.param )
{
  # add cytometer-specific parameters
  autosp.param$cytometer <- "FACSDiscover A8"

  autosp.param$scatter.data.min.x <- 0

  autosp.param$scatter.data.max.x <- 24140237

  autosp.param$scatter.data.min.y <- 0

  autosp.param$scatter.data.max.y <- 16488107

  autosp.param$expr.data.min <- -111

  autosp.param$expr.data.max <- 24140237

  autosp.param$default.scatter.parameter <- c( "FSC-A", "SSC (Violet)-A" )

  autosp.param$default.transformation.param <- list(
          length = 256,
          max.range = 2147483648.0,
          pos = 8.33,
          neg = 0,
          width = -500
        )

  autosp.param$time.and.scatter <- c(
    "FSC-A", "FSC-H", "FSC-W",
    "SSC (Violet)-A", "SSC (Violet)-H", "SSC (Violet)-W",
    "SSC (Imaging)-A", "SSC (Imaging)-H", "SSC (Imaging)-W",
    "LightLoss (Imaging)-A", "LightLoss (Imaging)-H","LightLoss (Imaging)-W",
    "LightLoss (Violet)-A", "LightLoss (Violet)-H", "LightLoss (Violet)-W"
  )

  autosp.param$non.spectral.channel <- c(
    "Time", "SSC", "FSC", "-H", "-W", "-T",
    "Delta", "Plate",
    "Radial", "Correlation", "Intensity",
    "Eccentricity", "Diffusivity", "Center",
    "Moment", "Size", "LightLoss",
    "Saturated", "Sorted", "Row", "Column",
    "Img", "Protocol", "EventLabel",
    "Region", "Gate", "Index", "Phase",
    "Event", "Drop", "Spectral", "Waveform",
    "Merged", "Flow", "Packet", "Reserved"
  )

  autosp.param$spectral.channel <- "[0-9]\\)-A"

  autosp.param$af.channel <- "V6 (515)-A"

  autosp.param$data.step <- 5e6

  autosp.param$large.gate.scaling.x <- 2.5
  autosp.param$large.gate.scaling.y <- 6

  autosp.param$ribbon.breaks <- c( -1e3, 0, 1e3, 1e4, 1e5, 1e6, 1e7 )

  return( autosp.param )

}

