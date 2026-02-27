# unmix_fcs.r

#' @title Unmix FCS Data
#'
#' @description
#' This function performs spectral unmixing on FCS data using various methods.
#'
#' @importFrom lifecycle deprecate_warn
#' @importFrom parallelly availableCores
#'
#' @param fcs.file A character string specifying the path to the FCS file.
#' @param spectra A matrix containing the spectral data. Fluorophores in rows,
#' detectors in columns.
#' @param asp The AutoSpectral parameter list. Prepare using
#' `get.autospectral.param`.
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
#' @param weights Optional numeric vector of weights (one per fluorescent
#' detector). Default is `NULL`, in which case weighting will be done by
#' channel means (Poisson variance). Only used for `WLS`.
#' @param af.spectra Spectral signatures of autofluorescences, normalized
#' between 0 and 1, with fluorophores in rows and detectors in columns. Prepare
#' using `get.af.spectra`. Required for `AutoSpectral` unmixing. Default is
#' `NULL` and will thus provoke failure if no spectra are provided and
#' `AutoSpectral` is selected.
#' @param spectra.variants Named list (names are fluorophores) carrying matrices
#' of spectral signature variations for each fluorophore. Prepare using
#' `get.spectral.variants`. Default is `NULL`. Used for
#' AutoSpectral unmixing. Required for per-cell fluorophore optimization.
#' @param output.dir A character string specifying the directory to save the
#' unmixed FCS file. Default is `NULL`, which will use `./AutoSpectral_unmixed`.
#' @param file.suffix A character string to append to the output file name.
#' Default is `NULL`.
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
#' default being `fast`. As of version 1.0.0, the backend for how this works
#' has changed. Spectral variants and AF signatures are now pre-screened per cell
#' to identify likely candidates, so brute force testing of all variants is no
#' longer required. So, `speed` controls the number of variants to be tested per
#' cell, with `fast` testing a single variant, `medium` testing 3 variants, and
#' `slow` testing 10 variants. While this is now implemented in pure R in
#' `AutoSpectral`, installation of `AutoSpectralRcpp` is strongly encouraged for
#' faster processing.
#' @param parallel Logical, default is `TRUE`, which enables parallel processing
#' for per-cell unmixing methods.
#' @param threads Numeric, defaults to a single thread for sequential processing
#' (`parallel=FALSE`) or all available cores if `parallel=TRUE`.
#' @param verbose Logical, controls messaging. Default is `TRUE`. Set to `FALSE`
#' to have it shut up.
#' @param n.variants Number of variants to test per cell. Allows explicit control
#' over the number used, as opposed to `speed`, which selects from pre-defined
#' choices. Providing a numeric value to `n.variants` will override `speed`,
#' allowing up to `n.variants` (or the max available) variants to be tested. The
#' default is `NULL`, in which case `n.variants` will be ignored.
#' @param chunk.size Numeric, number of events to use per chunk of unmixing. Used
#' to manage memory when processing large FCS files. As a rough guide, you will
#' need approximately 10x the size of the raw FCS file on disk as available
#' memory. Default is set at `2e6` events, assuming ~20GB memory available.
#' @param ... Ignored. Used to catch deprecated arguments.
#'
#' @return None. The function writes the unmixed FCS data to a file.
#'
#' @export

unmix.fcs <- function(
    fcs.file,
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
    speed = c("fast", "medium", "slow"),
    parallel = TRUE,
    threads = if ( parallel ) 0 else 1,
    verbose = TRUE,
    n.variants = NULL,
    chunk.size = 2e6,
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
      message(
        paste(
          "For AutoSpectral unmixing, providing fluorophore variation with",
          "`spectra.variants` may give better results. See `?get.spectral.variants`."
        )
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
    if ( is.null( n.variants ) ) {
      n.variants <- switch( match.arg( speed ), "fast" = 1L, "medium" = 3L, "slow" = 10L )
    }

    # remove any existing "AF" parameter in spectra
    if ( "AF" %in% rownames( spectra ) )
      spectra <- spectra[ rownames( spectra ) != "AF", , drop = FALSE ]
  }

  # create output folder if it doesn't exist
  if ( is.null( output.dir ) ) output.dir <- asp$unmixed.fcs.dir
  if ( !dir.exists( output.dir ) ) dir.create( output.dir )

  # import FCS, without warnings for fcs 3.2
  if ( verbose ) message( "Reading FCS metadata: ", fcs.file )
  import.meta <- readFCS( fcs.file, return.keywords = TRUE, start.row = 1, end.row = 1 )
  fcs.keywords <- import.meta$keywords
  total.events <- as.numeric( fcs.keywords[[ "$TOT" ]] )
  original.param <- colnames( import.meta$data )

  # determine original file name
  file.name <- if ( !is.null( fcs.keywords$`$FIL` ) ) fcs.keywords$`$FIL` else basename( fcs.file )
  if ( !grepl( "\\.fcs$", file.name, ignore.case = TRUE ) )  file.name <- paste0( file.name, ".fcs" )

  # deal with manufacturer peculiarities in writing fcs files
  if ( asp$cytometer %in% c( "ID7000", "Mosaic" ) ) {
    file.name <- sub(
      "([ _])Raw(\\.fcs$|\\s|$)",
      paste0("\\1", method, "\\2"),
      file.name,
      ignore.case = TRUE
    )
  } else if ( grepl( "Discover", asp$cytometer ) ) {
    #file.name <- fcs.keywords$FILENAME
    #file.name <- sub( ".*\\/", "", file.name )
    file.name <- sub( ".fcs", paste0( " ", method, ".fcs" ), file.name )
  } else {
    file.name <- sub( ".fcs", paste0( " ", method, ".fcs" ), file.name )
  }

  if ( !is.null( file.suffix ) )
    file.name <- sub( ".fcs", paste0( " ", file.suffix, ".fcs" ), file.name )

  # extract spectral data
  spectral.channel <- colnames( spectra )
  other.channels <- setdiff( original.param, spectral.channel )

  # remove height and width if present
  for ( ch in spectral.channel[ grepl( "-A$", spectral.channel ) ] ) {
    base <- sub( "-A$", "", ch )
    other.channels <- setdiff( other.channels, paste0( base, c( "-H", "-W" ) ) )
  }
  if ( grepl("Discover", asp$cytometer ) && !include.imaging ) {
    other.channels <- intersect( other.channels, asp$time.and.scatter )
  }

  # set multithreading
  if ( parallel & is.null( threads ) )
    threads <- asp$worker.process.n
  else if ( parallel & threads == 0 )
    threads <- parallelly::availableCores()

  # apply unmixing using selected method ---------------
  # pre-allocate the full results matrix
  extra.cols <- if ( method == "AutoSpectral" ) 2 else 0
  cols.n <- length( other.channels ) + nrow( spectra ) + extra.cols
  final.matrix <- matrix( 0, nrow = total.events, ncol = cols.n )

  # how many chunks do we need?
  chunk.n <- ceiling( total.events / chunk.size )

  # track if we've set the column names
  colnames.set <- FALSE

  # define weights if needed
  if ( ( weighted || method %in% c( "WLS", "Poisson", "FastPoisson" ) ) && is.null( weights ) ) {
    # Poisson-like weighting based on mean expression (variance = mean in a Poisson distribution)
    weight.sample <- readFCS( fcs.file, start.row = 1, end.row = min( 1e5, total.events ) )
    # weights are inverse of signal (more signal, more noise, less reliable)
    weights <- 1 / pmax( abs( colMeans( weight.sample[ , spectral.channel ] ) ), 1e-6 )
    rm( weight.sample )
  }

  # unmix in chunks for big files
  for ( i in 1:chunk.n ) {
    s.row <- ( (i - 1) * chunk.size ) + 1
    e.row <- min( i * chunk.size, total.events )

    if ( verbose )
      message( sprintf( "Processing chunk %d/%d (Events %d to %d)", i, chunk.n, s.row, e.row ) )

    # read in only events from this chunk
    chunk.data <- readFCS(
      fcs.file,
      return.keywords = FALSE,
      start.row = s.row,
      end.row = e.row
    )
    chunk.spectral <- chunk.data[ , spectral.channel, drop = FALSE ]
    chunk.other <- chunk.data[ , other.channels, drop = FALSE ]

    # unmix this chunk of data with the selected unmixing method
    unmixed.chunk <- switch(
      method,
      "OLS" = unmix.ols( chunk.spectral, spectra ),
      "WLS" = unmix.wls( chunk.spectral, spectra, weights ),
      "AutoSpectral" = {
        if ( requireNamespace("AutoSpectralRcpp", quietly = TRUE ) &&
             "unmix.autospectral.rcpp" %in% ls( getNamespace( "AutoSpectralRcpp" ) ) ) {
          tryCatch(
            AutoSpectralRcpp::unmix.autospectral.rcpp(
              raw.data = chunk.spectral,
              spectra = spectra,
              af.spectra = af.spectra,
              spectra.variants = spectra.variants,
              use.dist0 = use.dist0,
              verbose = verbose,
              speed = speed,
              parallel = parallel,
              threads = threads,
              n.variants = n.variants
            ),
            error = function( e ) {
              warning(
                "AutoSpectralRcpp unmixing failed, falling back to standard AutoSpectral: ",
                e$message,
                call. = FALSE
              )
              unmix.autospectral(
                raw.data = chunk.spectral,
                spectra = spectra,
                af.spectra = af.spectra,
                asp = asp,
                spectra.variants = spectra.variants,
                use.dist0 = use.dist0,
                verbose = verbose,
                speed = speed,
                parallel = parallel,
                threads = threads,
                n.variants = n.variants
              )
            }
          )
        } else {
          warning( "AutoSpectralRcpp not available, falling back to standard AutoSpectral" )
          unmix.autospectral(
            raw.data = chunk.spectral,
            spectra = spectra,
            af.spectra = af.spectra,
            asp = asp,
            spectra.variants = spectra.variants,
            use.dist0 = use.dist0,
            verbose = verbose,
            speed = speed,
            parallel = parallel,
            threads = threads,
            n.variants = n.variants
          )
        }
      },
      "Poisson" = unmix.poisson( chunk.spectral, spectra, asp, weights ),
      "FastPoisson" = {
        if ( requireNamespace("AutoSpectralRcpp", quietly = TRUE ) &&
             "unmix.poisson.fast" %in% ls( getNamespace( "AutoSpectralRcpp" ) ) ) {
          tryCatch(
            AutoSpectralRcpp::unmix.poisson.fast(
              raw.data = chunk.spectral,
              spectra = spectra,
              weights = weights,
              maxit = asp$rlm.iter.max,
              tol = 1e-6,
              n_threads = threads,
              divergence.threshold = divergence.threshold,
              divergence.handling = divergence.handling,
              balance.weight = balance.weight
            ),
            error = function( e ) {
              warning(
                "FastPoisson failed, falling back to standard Poisson: ",
                e$message,
                call. = FALSE
              )
              unmix.poisson(
                raw.data = chunk.spectral,
                spectra = spectra,
                asp = asp,
                initial.weights = weights,
                parallel = parallel,
                threads = threads
              )
            }
          )
        } else {
          warning( "AutoSpectralRcpp not available, falling back to standard Poisson.",
                   call. = FALSE )
          unmix.poisson(
            raw.data = chunk.spectral,
            spectra = spectra,
            asp = asp,
            initial.weights = weights,
            parallel = parallel,
            threads = threads
          )
        }
      },
      stop( "Unknown method" )
    )

    if ( !colnames.set ) {
      colnames( final.matrix ) <- c( other.channels, colnames( unmixed.chunk ) )
      colnames.set <- TRUE
    }

    # store in Pre-allocated Matrix
    final.matrix[ s.row:e.row, 1:ncol( chunk.other ) ] <- chunk.other
    final.matrix[ s.row:e.row, ( ncol( chunk.other ) + 1 ):cols.n ] <- unmixed.chunk

    # cleanup memory immediately
    rm( chunk.data, chunk.spectral, chunk.other, unmixed.chunk )
    if (i %% 5 == 0) gc()
  }

  # fix any NA values (e.g., plate location with S8)
  if ( anyNA( final.matrix ) ) final.matrix[ is.na( final.matrix ) ] <- 0

  # append "-A" to fluorophore and AF channel names
  colnames( final.matrix ) <- ifelse(
    colnames( final.matrix ) %in% c( rownames( spectra ), "AF" ),
    paste0( colnames( final.matrix ), "-A" ),
    colnames( final.matrix )
  )

  # update keywords
  new.keywords <- define.keywords(
    fcs.keywords,
    final.matrix,
    original.param,
    spectra,
    af.spectra,
    flow.control,
    asp,
    method,
    file.name,
    weights,
    spectral.channel
  )

  # save file
  if ( verbose ) message( paste( "Writing:", file.name ) )
  writeFCS( final.matrix, new.keywords, file.name, output.dir )
}
