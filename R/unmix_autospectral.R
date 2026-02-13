# unmix_autospectral.r

#' @title Unmix AutoSpectral
#'
#' @description
#' Unmix using the AutoSpectral method to extract autofluorescence and optimize
#' fluorophore signatures at the single cell level.
#'
#' @importFrom lifecycle deprecate_warn
#' @importFrom parallelly availableCores
#'
#' @param raw.data Expression data from raw fcs files. Cells in rows and
#' detectors in columns. Columns should be fluorescent data only and must
#' match the columns in spectra.
#' @param spectra Spectral signatures of fluorophores, normalized between 0
#' and 1, with fluorophores in rows and detectors in columns.
#' @param af.spectra Spectral signatures of autofluorescences, normalized
#' between 0 and 1, with fluorophores in rows and detectors in columns. Prepare
#' using `get.af.spectra`.
#' @param asp The AutoSpectral parameter list.
#' @param spectra.variants Named list (names are fluorophores) carrying matrices
#' of spectral signature variations for each fluorophore. Prepare using
#' `get.spectral.variants`. Default is `NULL`.
#' @param use.dist0 Logical, controls whether the selection of the optimal AF
#' signature for each cell is determined by which unmixing brings the fluorophore
#' signals closest to 0 (`use.dist0` = `TRUE`) or by which unmixing minimizes the
#' per-cell residual (`use.dist0` = `FALSE`). Default is `TRUE`.
#' @param verbose Logical, default `TRUE`. Whether to send messages to the console.
#' @param speed Default is `fast`. Selector for the precision-speed trade-off in
#' AutoSpectral per-cell fluorophore optimization. Options are `slow`, `medium`
#' and `fast`. From v1.0.0, this controls the number of variants tested per cell
#' (and per fluorophore). More variants takes longer, but gives better resolution
#' in some unmixed data. When `speed = fast`, as single variant will be tested;
#' for `medium`, three will be tested and for `slow`, 10 variants will be tested.
#' From AutoSpectral v1.0.0, all options are available in the pure R version.
#' Installation of `AutoSpectralRcpp` is strongly encouraged for speed, though.
#' @param parallel Logical, default is `TRUE`. The new parallel processing
#' should always be faster.
#' @param threads Numeric, default is `NULL`, in which case `asp$worker.process.n`
#' will be used. `asp$worker.process.n` is set by default to be one less than the
#' available cores on the machine. Multi-threading is only used if `parallel` is
#' `TRUE`.
#' @param n.variants Number of variants to test per cell. Allows explicit control
#' over the number used, as opposed to `speed`, which selects from pre-defined
#' choices. Providing a numeric value to `n.variants` will override `speed`,
#' allowing up to `n.variants` (or the max available) variants to be tested. The
#' default is `NULL`, in which case `n.variants` will be ignored.
#' @param ... Ignored. Previously used for deprecated arguments such as
#' `calculate.error`.
#'
#' @return Unmixed data with cells in rows and fluorophores in columns.
#'
#' @export

unmix.autospectral <- function(
    raw.data,
    spectra,
    af.spectra,
    asp,
    spectra.variants = NULL,
    use.dist0 = TRUE,
    verbose = TRUE,
    speed = c("fast", "medium", "slow"),
    parallel = TRUE,
    threads = NULL,
    n.variants = NULL,
    ...
) {

  # warn regarding deprecated arguments
  dots <- list( ... )

  if ( !is.null( dots$calculate.error ) ) {
    lifecycle::deprecate_warn(
      "0.9.1",
      "unmix.autospectral(calculate.error)",
      details = "no longer used"
    )
  }


  #############################################
  ### Autofluorescence Optimization Section ###
  #############################################

  # check for AF in spectra, remove if present
  if ( "AF" %in% rownames( spectra ) )
    spectra <- spectra[ rownames( spectra ) != "AF", , drop = FALSE ]

  if ( is.null( af.spectra ) )
    stop( "Multiple AF spectra must be provided.",
          call. = FALSE )
  if ( nrow( af.spectra ) < 2 )
    stop( "Multiple AF spectra must be provided.",
          call. = FALSE )

  # check for data/spectra column matching
  raw.data.cols <- colnames( raw.data )
  spectra.cols <- colnames( spectra )

  if ( !identical( raw.data.cols, spectra.cols ) ) {

    # ensure both actually have the same columns before reordering
    if ( all( spectra.cols %in% raw.data.cols ) &&
         length( spectra.cols ) == length( raw.data.cols ) ) {
      # reorder raw.data to match the order of spectra
      raw.data <- raw.data[, spectra.cols]
      message( "Columns reordered to match spectra." )
    } else {
      stop( "Column names in spectra and raw.data do not match perfectly;
           cannot reorder by name alone." )
    }
  }

  if ( verbose ) message( "Assigning optimal AF spectrum to each cell" )

  # set up for per-cell AF extraction
  fluorophores <- rownames( spectra )
  fluorophore.n <- nrow( spectra )
  af.idx.in.spectra <- fluorophore.n + 1
  detector.n <- ncol( spectra )
  combined.spectra <- matrix(
    NA_real_,
    nrow = af.idx.in.spectra,
    ncol = detector.n
  )
  colnames( combined.spectra ) <- colnames( spectra )
  fluors.af <- c( fluorophores, "AF" )
  rownames( combined.spectra ) <- fluors.af
  combined.spectra[ 1:fluorophore.n, ] <- spectra
  af.only <- is.null( spectra.variants )
  af.n <- nrow( af.spectra )

  # score the AF spectra per cell to determine the initial best "AF Index"
  if ( use.dist0 ) {
    # use minimization of total fluorophore signal (worst case scenario) as the metric
    af.assignments <- assign.af.fluorophores(
      raw.data,
      spectra,
      af.spectra
    )
  } else {
    # use minimization of/alignment to residuals as the metric
    af.assignments <- assign.af.residuals(
      raw.data,
      spectra,
      af.spectra
    )
  }

  if ( verbose ) message( "Unmixing autofluorescence" )

  # create empty matrix for collection of unmixed data
  unmixed <- matrix(
    0,
    nrow = nrow( raw.data ),
    ncol = fluorophore.n + 2
  )
  colnames( unmixed ) <- c( fluors.af, "AF Index" )

  fitted.af <- matrix(
    0,
    nrow = nrow( raw.data ),
    ncol = ncol( spectra )
  )

  # add AF assignments
  unmixed[ , "AF Index" ] <- af.assignments

  # perform initial unmixing using these AF assignments
  for ( af in seq_len( af.n ) ) {
    # set this AF as the spectrum to use
    combined.spectra[ af.idx.in.spectra, ] <- af.spectra[ af, ]

    # get the cells using this AF
    cell.idx <- which( af.assignments == af )

    if ( length( cell.idx ) > 0 ) {
      # unmix with this AF
      unmixed.af <- unmix.ols.fast(
        raw.data[ cell.idx, , drop = FALSE ],
        combined.spectra
      )
      # update the raw autofluorescence fit if performing fluorophore optimization
      if ( !af.only ) {
        fitted.af[ cell.idx, ] <- unmixed.af[ , "AF", drop = FALSE ] %*%
          af.spectra[ af, , drop = FALSE ]
      }

      # store unmixed data
      unmixed[ cell.idx, fluors.af ] <- unmixed.af
    }
  }

  # if we don't have spectral variants, stop here
  if ( is.null( spectra.variants ) ) return( unmixed )


  #############################################
  ##### Per-Cell Fluorophore Optimization #####
  #############################################

  # set positivity thresholds vector
  pos.thresholds <- rep( Inf, fluorophore.n )
  names( pos.thresholds ) <- fluorophores
  # fill with data
  pos.thresholds[ names( spectra.variants$thresholds ) ] <- spectra.variants$thresholds

  # unpack spectral variants
  variants <- spectra.variants$variants
  delta.list <- spectra.variants$delta.list
  delta.norms <- spectra.variants$delta.norms

  if ( is.null( pos.thresholds ) )
    stop( "Check that spectral variants have been calculated using
          `get.spectra.variants()`",
          call. = FALSE )
  if ( is.null( variants ) )
    stop( "Multiple fluorophore spectral variants must be provided.,
        call. = FALSE" )
  if ( !( length( variants ) > 1 ) )
    stop( "Multiple fluorophore spectral variants must be provided.",
          call. = FALSE )

  # use AF-subtracted raw data as input for fluorophore optimization
  remaining.raw <- raw.data - fitted.af

  # restrict optimization to fluors present in names( variants )
  optimize.fluors <- fluorophores[ fluorophores %in% names( variants ) ]

  # if delta.list and delta.norms are not provided by AutoSpectral (<v1.0.0), calculate
  if ( is.null( delta.list ) ) {
    message(
      paste(
        "For best results, re-calculate `spectra.variants` using AutoSpectral",
        "version 1.0.0 or greater. See `?get.spectra.variants`."
      )
    )
    # calculate deltas for each fluorophore's variants
    delta.list <- lapply( optimize.fluors, function( fl ) {
      variants[[ fl ]] - matrix(
        spectra[ fl, ],
        nrow = nrow( variants[[ fl ]] ),
        ncol = detector.n,
        byrow = TRUE
      )
    } )
    names( delta.list ) <- optimize.fluors

    # precompute delta norms
    delta.norms <- lapply( delta.list, function( d ) {
      sqrt( rowSums( d^2 ) )
    } )
  }

  # check that we have valid spectral variants
  optimize.fluors <- sanitize.optimization.inputs(
    spectra,
    optimize.fluors,
    variants,
    delta.norms
  )

  # if no match, return unmixed with warning
  if ( !( length( optimize.fluors ) > 0 ) ) {
    warning(
      "No matching fluorophores between supplied spectra and spectral variants.
      No spectral optimization performed.",
      call. = FALSE
    )
    return( unmixed )
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

  # set number of threads to use
  if ( parallel ) {
    if ( is.null( threads ) ) threads <- asp$worker.process.n
    if ( threads == 0 ) threads <- parallelly::availableCores()
  } else {
    threads <- 1
  }

  # pre-calculate indices rather than using names
  fluorophores <- which( rownames( combined.spectra ) %in% fluorophores )

  if ( verbose ) message( "Optimizing fluorophore unmixing cell-by-cell" )

  # call per-cell unmixing function
  unmixed[ , fluorophores ] <- optimize.unmix(
    raw.data = remaining.raw,
    unmixed = unmixed[ , fluorophores, drop = FALSE ],
    spectra = spectra,
    pos.thresholds = pos.thresholds,
    optimize.fluors = optimize.fluors,
    variants = variants,
    delta.list = delta.list,
    delta.norms = delta.norms,
    fluorophores = fluorophores,
    asp = asp,
    k = n.variants,
    nthreads = threads,
    parallel = parallel
  )

  return( unmixed )

}

