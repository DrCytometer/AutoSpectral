# get_autospectral_param.r

#' @title Get AutoSpectral Parameters
#'
#' @description Retrieves autospectral parameters for a specified cytometer.
#'
#' @param cytometer The type of cytometer, default is `aurora`. Supported options
#' include `aurora`, `auroraNL` for Northern Lights, `id7000`, `discover` (BD
#' FACSDiscover family), `a8` and `s8` (specific FACSDiscover models), `a5se`,
#' `opteon`, `mosaic`, `xenith` and `cytostellar`. Matching is case-insensitive
#' and supports unambiguous partial matches (e.g. `"aur"` matches `"aurora"`).
#' @param figures Logical indicating whether to set up directory parameters for
#' figures and tables, default is `TRUE`
#'
#' @return A list of AutoSpectral parameters.
#'
#' @export

get.autospectral.param <- function(
    cytometer = "aurora",
    figures = TRUE
) {

  # start with base parameter set
  autosp.param <- get.autospectral.param.minimal()

  if ( figures ) {

    autosp.param$figures <- TRUE

    # directory parameters
    autosp.param$figure.scatter.dir.base <- "figure_scatter"

    autosp.param$figure.gate.dir <- "figure_gate"
    autosp.param$figure.af.dir <- "figure_autofluorescence"
    autosp.param$figure.clean.control.dir <- "figure_clean_controls"
    autosp.param$figure.spectral.ribbon.dir <- "figure_spectral_ribbon"
    autosp.param$figure.spectra.dir <- "figure_spectra"
    autosp.param$figure.similarity.heatmap.dir <- "figure_similarity_heatmap"

    autosp.param$table.spectra.dir <- "table_spectra"

  }

  # cytometer aliases: alias name -> underlying get.autospectral.param.<fun>
  # plus an optional cytometer label override applied after that function runs.
  # A8 and S8 share identical acquisition parameters (get.autospectral.param.discover),
  # but downstream code distinguishes the two models by label, so the specific
  # label is restored here after dispatch.
  cytometer.synonym <- list(
    a8       = list( fun = "discover", label = "FACSDiscover A8" ),
    s8       = list( fun = "discover", label = "FACSDiscover S8" ),
    discover = list( fun = "discover", label = NULL )
  )

  # canonical cytometer names dispatched directly (function name == suffix)
  cytometer.direct <- c(
    "aurora", "auroraNL", "id7000", "a5se", "opteon", "mosaic", "xenith", "cytostellar"
  )

  cytometer.valid <- c( cytometer.direct, names( cytometer.synonym ) )

  # case-insensitive, unambiguous partial matching (match.arg-style)
  match.idx <- pmatch( tolower( cytometer ), tolower( cytometer.valid ) )

  if ( is.na( match.idx ) )
    stop(
      sprintf(
        "unsupported cytometer '%s'; supported values are: %s",
        cytometer, paste( cytometer.valid, collapse = ", " )
      ),
      call. = FALSE
    )

  cytometer.resolved <- cytometer.valid[ match.idx ]

  if ( cytometer.resolved %in% names( cytometer.synonym ) ) {
    cytometer.fun.suffix   <- cytometer.synonym[[ cytometer.resolved ]]$fun
    cytometer.label.override <- cytometer.synonym[[ cytometer.resolved ]]$label
  } else {
    cytometer.fun.suffix   <- cytometer.resolved
    cytometer.label.override <- NULL
  }

  # cytometer-specific parameters
  get.param.function <- get0( sprintf( "get.autospectral.param.%s", cytometer.fun.suffix ) )

  if ( is.null( get.param.function ) )
    stop( "unsupported cytometer", call. = FALSE )

  autosp.param <- get.param.function( autosp.param )

  if ( !is.null( cytometer.label.override ) )
    autosp.param$cytometer <- cytometer.label.override

  message( paste( "Cytometer set to", autosp.param$cytometer ) )

  return( autosp.param )

}
