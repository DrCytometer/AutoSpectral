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
#'   detectors in columns.
#' @param asp The AutoSpectral parameter list. Prepare using
#'   `get.autospectral.param`.
#' @param flow.control A list containing flow cytometry control parameters.
#' @param method A character string specifying the unmixing method to use. The
#'   default as of version 1.0.0 is now `AutoSpectral` to avoid confusion. To use
#'   AutoSpectral unmixing, you must provide at least `af.spectra` to perform
#'   autofluorescence extraction (on a per-cell basis). To also optimize
#'   fluorophore spectra, provide `spectra.variants`. To perform other types of
#'   unmixing, select from the options: `OLS`, `WLS`, `Poisson` or `FastPoisson`.
#'   `FastPoisson` requires installation of `AutoSpectralRcpp`.
#' @param weighted Logical, whether to use ordinary or weighted least squares
#'   unmixing as the base algorithm in AutoSpectral legacy pipeline unmixing.
#'   Default is `FALSE` and will use OLS.
#' @param weights Optional numeric vector of weights (one per fluorescent
#'   detector). Default is `NULL`, in which case weighting will be done by
#'   channel means (Poisson variance). Only used for `WLS`.
#' @param af.spectra Spectral signatures of autofluorescences, normalized
#'   between 0 and 1, with fluorophores in rows and detectors in columns. Prepare
#'   using `get.af.spectra`. Required for `AutoSpectral` unmixing. Default is
#'   `NULL` and will thus provoke failure if no spectra are provided and
#'   `AutoSpectral` is selected.
#' @param spectra.variants Named list (names are fluorophores) carrying matrices
#'   of spectral signature variations for each fluorophore. Prepare using
#'   `get.spectral.variants`. Default is `NULL`. Used for
#'   AutoSpectral unmixing. Required for per-cell fluorophore optimization.
#' @param output.dir A character string specifying the directory to save the
#'   unmixed FCS file. Default is `NULL`, which will use `./AutoSpectral_unmixed`.
#' @param file.suffix A character string to append to the output file name.
#'   Default is `NULL`.
#' @param include.raw A logical value indicating whether to include raw
#'   expression data in the written FCS file. Default is `FALSE` to provide smaller
#'   output files.
#' @param include.imaging A logical value indicating whether to include imaging
#'   parameters in the written FCS file. Default is `TRUE`.
#' @param use.dist0 Legacy pipeline argument. Logical, controls whether the
#'   selection of the optimal AF signature for each cell is determined by the
#'   minimization of potential AF spillover into the fluorophore channels
#'   (`use.dist0` = `TRUE`) or by which unmixing minimizes the per-cell residual
#'   (`use.dist0` = `FALSE`). Default is `TRUE`. Used for legacy AutoSpectral
#'   unmixing.
#' @param divergence.threshold Numeric. Used for `FastPoisson` only.
#'   Threshold to trigger reversion towards WLS unmixing when Poisson result
#'   diverges for a given point. To be deprecated.
#' @param divergence.handling String. How to handle divergent cells from Poisson
#'   IRLS. Options are `NonNeg` (non-negativity will be enforced), `WLS` (revert
#'   to WLS initial unmix) or `Balance` (WLS and NonNeg will be averaged).
#'   Default is `Balance`. To be deprecated.
#' @param balance.weight Numeric. Weighting to average non-convergent cells.
#'   Used for `Balance` option under `divergence.handling`. Default is `0.5`.
#'   To be deprecated.
#' @param speed Selector for the precision-speed trade-off in AutoSpectral per-cell
#'   fluorophore optimization. Options are `fast`, `medium` and `slow`, with the
#'   default being `fast`. As of version 1.0.0, the backend for how this works
#'   has changed. Spectral variants and AF signatures are now pre-screened per cell
#'   to identify likely candidates, so brute force testing of all variants is no
#'   longer required. So, `speed` controls the number of variants to be tested per
#'   cell, with `fast` testing a single variant, `medium` testing 3 variants, and
#'   `slow` testing 10 variants. While this is now implemented in pure R in
#'   `AutoSpectral`, installation of `AutoSpectralRcpp` is strongly encouraged for
#'   faster processing.
#' @param parallel Logical, default is `TRUE`, which enables parallel processing
#'   for per-cell unmixing methods.
#' @param threads Numeric, defaults to a single thread for sequential processing
#'   (`parallel=FALSE`) or all available cores if `parallel=TRUE`.
#' @param verbose Logical, controls messaging. Default is `TRUE`. Set to `FALSE`
#'   to have it shut up.
#' @param n.variants Numeric, used for legacy AutoSpectral pipeline unmixing.
#'   Number of variants to test per cell. Allows explicit control over the
#'   number used, as opposed to `speed`, which selects from pre-defined choices.
#'   Providing a numeric value to `n.variants` will override `speed`, allowing
#'   up to `n.variants` (or the max available) variants to be tested. The default
#'   is `NULL`, in which case `n.variants` will be ignored.
#' @param chunk.size Numeric, number of events to use per chunk of unmixing. Used
#'   to manage memory when processing large FCS files. As a rough guide, you will
#'   need approximately 10x the size of the raw FCS file on disk as available
#'   memory. Default is set at `2e6` events, assuming ~20GB memory available.
#' @param pipeline Character, one of `"joint"` (default) or `"legacy"`. Passed
#'   to `unmix.autospectral.rcpp()`. `"joint"` uses the new covariance-weighted
#'   joint per-cell pipeline; `"legacy"` reproduces the behaviour of
#'   AutoSpectral prior to version 1.6.0.
#' @param n.passes Integer, default `2L`. Number of joint optimisation passes
#'   per cell. Only used when `pipeline = "joint"`.
#' @param n.af.passes Integer, default `1L`. Number of autofluorescence
#'   extraction passes per cell. Only used when `pipeline = "joint"`. Passed
#'   to `unmix.autospectral.rcpp()`.
#' @param cell.weight Logical, default `FALSE`. Applies per-cell detector
#'   weighting to the joint unmixing solve. Only used when
#'   `pipeline = "joint"`. Passed to `unmix.autospectral.rcpp()`. Useful for
#'   ID7000 files.
#' @param noise.floor Numeric, default `125`. Lower clamp on the denominator
#'   of the per-cell detector weights when `cell.weight = TRUE`. Only used
#'   when `pipeline = "joint"`. Passed to `unmix.autospectral.rcpp()`.
#' @param alpha Numeric, default `0.5`. Weighting for balancing residual and
#'   covariance spillover minimization. Only used when `pipeline = "joint"`.
#'   Passed to `unmix.autospectral.rcpp()`.
#' @param collinear.threshold Numeric, default `0.5`. Cosine similarity value
#'   to trigger conflict assessment for collinear fluorophore variants. Only
#'   used when `pipeline = "joint"`. Passed to `unmix.autospectral.rcpp()`.
#' @param joint.pair.resolution Logical, default `TRUE`. Whether to perform
#'   conflict-resolution for collinear fluorophore pairs. Only used when
#'   `pipeline = "joint"`. Passed to `unmix.autospectral.rcpp()`.
#' @param refine.af.quantile Numeric, default `0.5`. Fraction of cells taken
#'   forward for additional AF passes (see `n.af.passes`). Only used when
#'   `pipeline = "joint"`. Passed to `unmix.autospectral.rcpp()`.
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
    method = c("AutoSpectral", "OLS", "WLS", "Poisson", "FastPoisson"),
    weighted = FALSE,
    weights = NULL,
    af.spectra = NULL,
    spectra.variants = NULL,
    output.dir = NULL,
    file.suffix = NULL,
    include.raw = FALSE,
    include.imaging = TRUE,
    use.dist0 = TRUE,
    divergence.threshold = 1e4,
    divergence.handling = "Balance",
    balance.weight = 0.5,
    speed = c("fast", "medium", "slow"),
    parallel = FALSE,
    threads = NULL,
    verbose = TRUE,
    n.variants = NULL,
    chunk.size = 2e6,
    pipeline  = c( "joint", "legacy" ),
    n.passes  = 2L,
    n.af.passes            = 1L,
    cell.weight            = if (asp$cytometer == "ID7000") TRUE else FALSE,
    noise.floor            = 125,
    alpha                  = 0.5,
    collinear.threshold    = 0.5,
    joint.pair.resolution  = TRUE,
    refine.af.quantile     = 0.5,
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

  # handle deprecated "Automatic" method
  if ( identical( method, "Automatic" ) ) {
    lifecycle::deprecate_warn(
      "1.5.5",
      "unmix.fcs(method = 'Automatic')",
      details = paste(
        "Please specify the method explicitly.",
        "Use `method = 'AutoSpectral'` (with `af.spectra`),",
        "`method = 'WLS'` for ID7000/A8/S8, or `method = 'OLS'` otherwise."
      )
    )
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

  method <- match.arg( method )

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
    speed <- match.arg( speed )
    if ( is.null( n.variants ) || !is.numeric( n.variants ) || length( n.variants ) != 1 ) {
      n.variants <- switch(
        speed,
        "slow"   = 10L,
        "medium" = 3L,
        "fast"   = 1L
      )
    }
  }

  # set up, create output folders where FCS files will go
  if ( is.null( output.dir ) ) output.dir <- asp$unmixed.fcs.dir
  if ( !dir.exists( output.dir ) ) dir.create( output.dir )
  if ( is.null( threads ) ) threads <- asp$worker.process.n

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
    n.variants = n.variants,
    chunk.size = chunk.size,
    pipeline  = match.arg( pipeline ),
    n.passes  = n.passes,
    n.af.passes           = n.af.passes,
    cell.weight           = cell.weight,
    noise.floor           = noise.floor,
    alpha                 = alpha,
    collinear.threshold   = collinear.threshold,
    joint.pair.resolution = joint.pair.resolution,
    refine.af.quantile    = refine.af.quantile
  )

  # set up parallel processing
  if ( parallel && ( method == "OLS" || method == "WLS" ) ) {

    # force sequential processing within a cluster
    args.list$parallel <- FALSE
    args.list$threads <- 1
    args.list$verbose <- FALSE

    internal.functions <- c(
      "unmix.fcs", "unmix.ols", "unmix.wls", "unmix.autospectral",
      "unmix.poisson", "readFCS", "writeFCS", "define.keywords"
    )
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

  execute.unmix <- function( lapply.func ) {
    invisible(
      lapply.func( files.to.unmix, function( f ) {
        do.call( unmix.fcs, c( list( f ), args.list ) )
      } )
    )
  }

  # unmix all files in list with fallback for errors in parallel processing
  tryCatch( {
    if ( parallel && !is.null( result$cleanup ) ) {
      # attempt parallel processing
      tryCatch( {
        execute.unmix( lapply.function )
      }, error = function(e) {
        warning( "Parallel setup failed. Falling back to sequential. Error: ", e$message )
        execute.unmix( lapply )
      } )
    } else {
      # sequential processing
      execute.unmix( lapply )
    }
  }, finally = {
    # Clean up cluster when done if needed
    if ( !is.null( result$cleanup ) ) result$cleanup()
  } )

  return( invisible( TRUE ) )
}
