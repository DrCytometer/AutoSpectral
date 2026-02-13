# unmix_folder.r

#' @title Unmix All FCS Files in a Directory
#'
#' @description
#' This function unmixes all FCS files in a specified directory using the
#' provided spectra and method, and saves the unmixed FCS files to an output
#' directory of the user's choice.
#'
#' @importFrom lifecycle deprecate_warn
#'
#' @param fcs.dir Directory (file path) containing FCS files to be unmixed.
#' @param spectra A matrix containing the spectral data. Fluorophores in rows,
#' detectors in columns.
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#' @param flow.control A list containing flow cytometry control parameters.
#' @param method A character string specifying the unmixing method to use. The
#' default as of version 1.0.0 is now `AutoSpectral` to avoid confusion. To use
#' AutoSpectral unmixing, you must provide at least `af.spectra` to perform
#' autofluorescence extraction (on a per-cell basis). To also optimize
#' fluorophore spectra, provide `spectra.variants`. To perform other types of
#' unmixing, select from the options: `OLS`, `WLS`, `Poisson` or `FastPoisson`.
#' `FastPoisson` requires installation of `AutoSpectralRcpp`.There is also
#' `Automatic`, which switches depending on the inputs provided: it uses
#' `AutoSpectral` for AF extraction if `af.spectra` are provided, and
#' automatically selects `OLS` or `WLS` depending on which is normal for the
#' given cytometer in `asp$cytometer`. This means that files from the ID7000,
#' A8 and S8 will be unmixed using `WLS` while others will be unmixed with `OLS`,
#' if AutoSpectral unmixing is not activated.
#' @param weighted Logical, whether to use ordinary or weighted least squares
#' unmixing as the base algorithm in AutoSpectral unmixing.
#' Default is `FALSE` and will use OLS.
#' @param weights Optional numeric vector of weights: one per fluorescent
#' detector. Default is `NULL`, in which case weighting will be done by
#' channel means. Only used for `WLS`
#' @param af.spectra Spectral signatures of autofluorescences, normalized
#' between 0 and 1, with fluorophores in rows and detectors in columns. Prepare
#' using `get.af.spectra`. Required for `AutoSpectral` unmixing. Default is
#' `NULL` and will thus provoke failure if no spectra are provided and
#' `AutoSpectral` is selected.
#' @param spectra.variants Named list (names are fluorophores) carrying matrices
#' of spectral signature variations for each fluorophore. Prepare using
#' `get.spectral.variants`. Default is `NULL`. Used for
#' AutoSpectral unmixing. Required for per-cell fluorophore optimization.
#' @param output.dir Directory to save the unmixed FCS files
#' (default is `asp$unmixed.fcs.dir`, which is `./AutoSpectral_unmixed`).
#' @param file.suffix A character string to append to the output file name.
#' Default is `NULL`
#' @param include.raw A logical value indicating whether to include raw
#' expression data in the written FCS file. Default is `FALSE` to provide smaller
#' output files.
#' @param include.imaging A logical value indicating whether to include imaging
#' parameters in the written FCS file. Default is `FALSE` to provide smaller
#' output files.
#' @param use.dist0 Logical, controls whether the selection of the optimal AF
#' signature for each cell is determined by which unmixing brings the fluorophore
#' signals closest to 0 (`use.dist0` = `TRUE`) or by which unmixing minimizes the
#' per-cell residual (`use.dist0` = `FALSE`). Default is `TRUE`. Used for
#' AutoSpectral unmixing. The minimization of fluorophore signals can be thought
#' of as a "worst-case" scenario, but it provides more accurate assignments,
#' particularly with large panels.
#' @param divergence.threshold Numeric. Used for `FastPoisson` only.
#' Threshold to trigger reversion towards WLS unmixing when Poisson result
#' diverges for a given point. To be deprecated.
#' @param divergence.handling String. How to handle divergent cells from Poisson
#' IRLS. Options are `NonNeg` (non-negativity will be enforced), `WLS` (revert
#' to WLS initial unmix) or `Balance` (WLS and NonNeg will be averaged).
#' Default is `Balance`. To be deprecated.
#' @param balance.weight Numeric. Weighting to average non-convergent cells.
#' Used for `Balance` option under `divergence.handling`. Default is `0.5`.
#' To be deprecated.
#' @param speed Selector for the precision-speed trade-off in AutoSpectral per-cell
#' fluorophore optimization. Options are `fast`, `medium` and `slow`, with the
#' default being `slow`. As of version 1.0.0, the backend for how this works
#' has changed. Spectral variants and AF signatures are now pre-screened per cell
#' to identify likely candidates, so brute force testing of all variants is no
#' longer required. So, `speed` controls the number of variants to be tested per
#' cell, with `fast` testing a single variant, `medium` testing 3 variants, and
#' `slow` testing 10 variants. While this is now implemented in pure R in
#' `AutoSpectral`, installation of `AutoSpectralRcpp` is strongly encouraged for
#' faster processing.
#' @param parallel Logical, default is `FALSE`. Set to `TRUE` to activate parallel
#' processing for multiple FCS files.
#' @param threads Numeric, default is `NULL`, in which case `asp$worker.process.n`
#' will be used. `asp$worker.process.n` is set by default to be one less than the
#' available cores on the machine. Multi-threading is only used if `parallel` is
#' `TRUE`. If working on a computing cluster, try `parallelly::availableCores()`.
#' @param verbose Logical, controls messaging. Default is `TRUE`. Set to `FALSE`
#' to have it shut up.
#' @param n.variants Number of variants to test per cell. Allows explicit control
#' over the number used, as opposed to `speed`, which selects from pre-defined
#' choices. Providing a numeric value to `n.variants` will override `speed`,
#' allowing up to `n.variants` (or the max available) variants to be tested. The
#' default is `NULL`, in which case `n.variants` will be ignored.
#' @param ... Ignored. Previously used for deprecated arguments such as
#' `calculate.error`.
#'
#' @return None. Saves the unmixed FCS files to the specified output directory.
#'
#' @export

unmix.folder <- function(
    fcs.dir,
    spectra,
    asp,
    flow.control,
    method = "AutoSpectral",
    weighted = FALSE,
    weights = NULL,
    af.spectra = NULL,
    spectra.variants = NULL,
    output.dir = NULL,
    file.suffix = NULL,
    include.raw = FALSE,
    include.imaging = FALSE,
    use.dist0 = TRUE,
    divergence.threshold = 1e4,
    divergence.handling = "Balance",
    balance.weight = 0.5,
    speed = c("slow", "medium", "fast"),
    parallel = FALSE,
    threads = NULL,
    verbose = TRUE,
    n.variants = NULL,
    ...
) {

  # warn regarding deprecated arguments
  dots <- list( ... )

  if ( !is.null( dots$calculate.error ) ) {
    lifecycle::deprecate_warn(
      "0.9.1",
      "unmix.fcs(calculate.error)",
      details = "The 'calculate.error' argument is no longer used."
    )
    dots$calculate.error <- NULL
  }

  # check for other odd stuff being passed
  if ( length( dots ) > 0 ) {
    stop(
      paste(
        "Unknown arguments detected:",
        paste( names( dots ), collapse = ", " )
      )
    )
  }

  # logic for default unmixing with cytometer-based selection
  if ( method == "Automatic" ) {
    if ( !is.null( af.spectra ) ) {
      method <- "AutoSpectral"
      if ( asp$cytometer %in% c( "FACSDiscover S8", "FACSDiscover A8", "ID7000" ) )
        weighted <- TRUE
    } else if ( asp$cytometer %in% c( "FACSDiscover S8", "FACSDiscover A8", "ID7000" ) ) {
      method <- "WLS"
    } else {
      method <- "OLS"
    }
  }

  # include checks on inputs if AutoSpectral unmixing has been selected
  if ( method == "AutoSpectral" ) {
    # check for af.spectra, stop if not
    if ( is.null( af.spectra ) ) {
      stop(
        "For AutoSpectral unmixing, `af.spectra` must be provided.
        See `?get.af.spectra()`.",
        call. = FALSE
      )
    }
    # check that af.spectra is a matrix and has rows
    if ( !is.matrix( af.spectra ) || nrow( af.spectra ) < 2 ) {
      stop(
        "For AutoSpectral unmixing, multiple AF in `af.spectra` must be provided
        in a matrix. See `?get.af.spectra`.",
        call. = FALSE
      )
    }
    # check for variants, warn if not provided
    if ( is.null( spectra.variants ) ) {
      warning(
        "For AutoSpectral unmixing, providing fluorophore variation with `spectra.variants`
        will give better results. See `?get.spectral.variants`.",
        call. = FALSE
      )
    }
    # check for AutoSpectralRcpp
    if ( requireNamespace( "AutoSpectralRcpp", quietly = TRUE ) ) {
      # require 1.0.0 or higher if installed for compatibility
      if ( utils::packageVersion( "AutoSpectralRcpp" ) < package_version( "1.0.0" ) ) {
        stop(
          "Package `AutoSpectralRcpp` >= 1.0.0 is required for this method.
          Please update the package.",
          call. = FALSE
        )
      }
    } else {
      # warn if not available (default to R processing)
      warning(
        "Package `AutoSpectralRcpp` not found. Please install it for faster processing.",
        call. = FALSE
      )
    }

    # set number of variants to test (by `speed` if `n.variants` is not provided)
    if ( length( speed ) > 1 )
      speed <- speed[ 1 ]

    if ( is.null( n.variants ) || !is.numeric( n.variants ) || length( n.variants ) != 1 ) {
      n.variants <- switch(
        speed,
        "slow"   = 10L,
        "medium" = 3L,
        "fast"   = 1L,
        {
          warning(
            paste0(
              "Unrecognized input '",
              speed,
              "' to `speed`. Defaulting to `slow` (n.variants=10)."
            ),
            call. = FALSE
          )
          10L
        }
      )
    }
  }

  # set up, create output folders where FCS files will go
  if ( is.null( output.dir ) )
    output.dir <- asp$unmixed.fcs.dir
  if ( !dir.exists( output.dir ) )
    dir.create( output.dir )

  if ( parallel & is.null( threads ) )
    threads <- asp$worker.process.n

  files.to.unmix <- list.files( fcs.dir, pattern = ".fcs", full.names = TRUE )

  # construct list of arguments
  args.list <- list(
    spectra = spectra,
    asp = asp,
    flow.control = flow.control,
    method = method,
    weighted = weighted,
    weights = weights,
    af.spectra = af.spectra,
    spectra.variants = spectra.variants,
    output.dir = output.dir,
    file.suffix = file.suffix,
    include.raw = include.raw,
    include.imaging = include.imaging,
    use.dist0 = use.dist0,
    divergence.threshold = divergence.threshold,
    divergence.handling = divergence.handling,
    balance.weight = balance.weight,
    speed = speed,
    parallel = parallel,
    threads = threads,
    verbose = verbose,
    n.variants = n.variants
  )

  # Set up parallel processing
  if ( parallel && ( method == "OLS" || method == "WLS" ) ) {
    internal.functions <- c( "unmix.fcs", "unmix.ols", "unmix.wls" )
    exports <- c( "args.list", "files.to.unmix", internal.functions )

    result <- create.parallel.lapply(
      asp,
      exports,
      parallel = parallel,
      threads = threads,
      export.env = environment()
    )
    lapply.function <- result$lapply
  } else {
    lapply.function <- lapply
    result <- list( cleanup = NULL )
  }

  # unmix all files in list
  tryCatch( {
    lapply.function( files.to.unmix, function( f ) {
      do.call( unmix.fcs, c( list( f ), args.list ) )
    } )
  }, finally = {
    # clean up cluster when done if needed
    if ( !is.null( result$cleanup ) ) result$cleanup()
  } )
}
