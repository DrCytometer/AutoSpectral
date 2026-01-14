# get_autospectral_param.r

#' @title Get AutoSpectral Parameters
#'
#' @description Retrieves autospectral parameters for a specified cytometer.
#'
#' @param cytometer The type of cytometer, default is `aurora`. Supported options
#' include `aurora`, `auroraNL` for Northern Lights, `id7000`, `a8`, `s8`, `a5se`,
#' `opteon`, `mosaic` and `xenith`.
#' @param figures Logical indicating whether to set up directory parameters for
#' figures and tables, default is `TRUE`
#'
#' @return A list of AutoSpectral parameters.
#'
#' @export

get.autospectral.param <- function( cytometer = "aurora", figures = TRUE )
{
  autosp.param <- get.autospectral.param.minimal()

  if ( figures ) {

    autosp.param$figures <- TRUE

    # directory parameters
    autosp.param$figure.scatter.dir.base <- "figure_scatter"

    autosp.param$figure.gate.dir <- "figure_gate"
    autosp.param$figure.af.dir <- "figure_autofluorescence"
    autosp.param$figure.peacoqc.dir <- "figure_peacoQC"
    autosp.param$figure.clean.control.dir <- "figure_clean_controls"
    autosp.param$figure.spectral.ribbon.dir <- "figure_spectral_ribbon"
    autosp.param$figure.spectra.dir <- "figure_spectra"
    autosp.param$figure.similarity.heatmap.dir <- "figure_similarity_heatmap"

    autosp.param$table.spectra.dir <- "table_spectra"

  }

  # cytometer-specific parameters
  get.param.function <- get0( sprintf( "get.autospectral.param.%s", cytometer ) )

  if ( is.null( get.param.function ) )
    stop( "unsupported cytometer", call. = FALSE )

  autosp.param <- get.param.function( autosp.param )

  return( autosp.param )

}

