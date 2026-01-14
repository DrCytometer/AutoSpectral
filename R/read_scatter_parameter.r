# read_scatter_parameter.r

#' @title Read Scatter Parameters
#'
#' @description
#' This function reads scatter parameters from a specified file or uses default
#' scatter parameters if the file is not available.
#'
#' @importFrom utils read.csv
#'
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#'
#' @return A vector of scatter parameters.

read.scatter.parameter <- function( asp )
{
    if ( ! is.null( asp$scatter.parameter.file.name ) &&
         file.exists( asp$scatter.parameter.file.name ) )
        scatter.parameter <- read.csv(
          asp$scatter.parameter.file.name,
          stringsAsFactors = FALSE
        )
    else
        scatter.parameter <- data.frame(
          parameter = asp$default.scatter.parameter,
          stringsAsFactors = FALSE
        )

    return( scatter.parameter$parameter )
}

