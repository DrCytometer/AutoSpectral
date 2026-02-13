# unmix_fcs.r

#' @title Unmix FCS Data
#'
#' @description
#' This function performs spectral unmixing on FCS data using various methods.
#'
#' @importFrom flowCore read.FCS keyword exprs flowFrame parameters
#' @importFrom flowCore write.FCS parameters<- keyword<-
#' @importFrom lifecycle deprecate_warn
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
#' default being `slow`. As of version 1.0.0, the backend for how this works
#' has changed. Spectral variants and AF signatures are now pre-screened per cell
#' to identify likely candidates, so brute force testing of all variants is no
#' longer required. So, `speed` controls the number of variants to be tested per
#' cell, with `fast` testing a single variant, `medium` testing 3 variants, and
#' `slow` testing 10 variants. While this is now implemented in pure R in
#' `AutoSpectral`, installation of `AutoSpectralRcpp` is strongly encouraged for
#' faster processing.
#' @param parallel Logical, default is `TRUE`, which enables parallel processing
#' for per-cell unmixing methods.
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
    speed = c("slow", "medium", "fast"),
    parallel = TRUE,
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

  # create output folder if it doesn't exist
  if ( is.null( output.dir ) )
    output.dir <- asp$unmixed.fcs.dir
  if ( !dir.exists( output.dir ) )
    dir.create( output.dir )

  if ( is.null( threads ) )
    threads <- asp$worker.process.n

  # import FCS, without warnings for fcs 3.2
  fcs.data <- suppressWarnings(
    flowCore::read.FCS(
      fcs.file,
      transformation = FALSE,
      truncate_max_range = FALSE,
      emptyValue = FALSE
    )
  )

  # get keyword information
  fcs.keywords <- flowCore::keyword( fcs.data )
  file.name <- flowCore::keyword( fcs.data, "$FIL" )

  # deal with manufacturer peculiarities in writing fcs files
  if ( asp$cytometer %in% c( "ID7000", "Mosaic" ) ) {
    file.name <- sub(
      "([ _])Raw(\\.fcs$|\\s|$)",
      paste0("\\1", method, "\\2"),
      file.name,
      ignore.case = TRUE
    )
  } else if ( grepl( "Discover", asp$cytometer ) ) {
    file.name <- fcs.keywords$FILENAME
    file.name <- sub( ".*\\/", "", file.name )
    file.name <- sub( ".fcs", paste0( " ", method, ".fcs" ), file.name )
  } else {
    file.name <- sub( ".fcs", paste0( " ", method, ".fcs" ), file.name )
  }

  if ( !is.null( file.suffix ) )
    file.name <- sub( ".fcs", paste0( " ", file.suffix, ".fcs" ), file.name )

  # extract exprs
  fcs.exprs <- flowCore::exprs( fcs.data )
  rm( fcs.data ) # free up memory
  original.param <- colnames( fcs.exprs )

  # extract spectral data
  spectral.channel <- colnames( spectra )
  spectral.exprs <- fcs.exprs[ , spectral.channel, drop = FALSE ]

  other.channels <- setdiff( colnames( fcs.exprs ), spectral.channel )

  # remove height and width if present
  suffixes <- c( "-H", "-W" )
  for ( ch in spectral.channel[ grepl( "-A$", spectral.channel ) ] ) {
    base <- sub( "-A$", "", ch )
    other.channels <- setdiff( other.channels, paste0( base, suffixes ) )
  }
  other.exprs <- fcs.exprs[ , other.channels, drop = FALSE ]

  if ( !include.raw )
    rm( fcs.exprs ) # free up memory

  # remove imaging parameters if desired
  if ( grepl( "Discover", asp$cytometer ) & !include.imaging )
    other.exprs <- other.exprs[ , asp$time.and.scatter ]

  # define weights if needed
  if ( weighted | method == "WLS"| method == "Poisson"| method == "FastPoisson" ) {
    if ( is.null( weights ) ) {
      # Poisson-like weighting based on mean expression
      # (variance = mean in a Poisson distribution)
      # enforce non-zero
      weights <- pmax( abs( colMeans( spectral.exprs ) ), 1e-6 )
      # weights are inverse of signal (more signal, more noise, less reliable)
      weights <- 1 / weights
    }
  }

  # apply unmixing using selected method ---------------
  unmixed.data <- switch(
    method,
    "OLS" = unmix.ols( spectral.exprs, spectra ),
    "WLS" = unmix.wls( spectral.exprs, spectra, weights ),
    "AutoSpectral" = {
      if ( requireNamespace("AutoSpectralRcpp", quietly = TRUE ) &&
           "unmix.autospectral.rcpp" %in% ls( getNamespace( "AutoSpectralRcpp" ) ) ) {
        tryCatch(
          AutoSpectralRcpp::unmix.autospectral.rcpp(
            raw.data = spectral.exprs,
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
              raw.data = spectral.exprs,
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
          raw.data = spectral.exprs,
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
    "Poisson" = unmix.poisson( spectral.exprs, spectra, asp, weights ),
    "FastPoisson" = {
      if ( requireNamespace("AutoSpectralRcpp", quietly = TRUE ) &&
           "unmix.poisson.fast" %in% ls( getNamespace( "AutoSpectralRcpp" ) ) ) {
        tryCatch(
          AutoSpectralRcpp::unmix.poisson.fast(
            raw.data = spectral.exprs,
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
              raw.data = spectral.exprs,
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
          raw.data = spectral.exprs,
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

  # add back raw exprs and others columns as desired
  if ( include.raw ) {
    unmixed.data <- cbind( fcs.exprs, unmixed.data )
  } else {
    unmixed.data <- cbind( other.exprs, unmixed.data )
  }

  rm( spectral.exprs, other.exprs ) # free up memory

  # fix any NA values (e.g., plate location with S8)
  if ( anyNA( unmixed.data ) )
    unmixed.data[ is.na( unmixed.data ) ] <- 0

  # update keywords----------
  # identify non-parameter keywords
  non.param.keys <- fcs.keywords[ !grepl( "^\\$?P\\d+", names( fcs.keywords ) ) ]
  if ( asp$cytometer == "Mosaic" )
    non.param.keys <- non.param.keys[ !grepl( "^\\$?CH\\d+", names( non.param.keys ) ) ]

  # build lookup
  pN.keys <- grep( "^\\$?P\\d+N$", names( fcs.keywords ), value = TRUE )
  param.lookup <- lapply( pN.keys, function( k ) {
    p.idx <- sub( "\\$?P(\\d+)N", "\\1", k )
    matches <- grep( paste0( "^\\$?P", p.idx, "(?:[A-Z]+)$" ), names( fcs.keywords ),
                     value = TRUE )
    stats::setNames( fcs.keywords[ matches], matches )
  } )
  # name the list by parameter name
  names( param.lookup ) <- sapply( pN.keys, function( k ) fcs.keywords[[ k ]])

  # keywords for new parameters
  param.keywords <- list()
  n.param <- ncol( unmixed.data )

  # check all parameters and update as needed
  for ( i in seq_len( n.param ) ) {
    p.name <- colnames( unmixed.data )[ i ]

    if ( p.name %in% original.param ) {
      # retain keywords from original file if present
      old.entry <- param.lookup[[ p.name ]]
      if ( !is.null( old.entry ) ) {
        # update index to current parameter number
        names( old.entry ) <- sub( "^\\$P\\d+", paste0( "$P", i ), names( old.entry ) )
        param.keywords <- c( param.keywords, old.entry )

      } else {
        # fallback if missing
        param.keywords[[ paste0( "$P", i, "N" )]] <- p.name
        param.keywords[[ paste0( "$P", i, "S" )]] <- p.name
      }

    } else {
      # keywords for new unmixed parameters
      bit.depth <- if ( !is.null( asp$bit.depth ) ) asp$bit.depth else "32"

      param.keywords[[ paste0( "$P", i, "N" ) ]] <- p.name
      param.keywords[[ paste0( "$P", i, "B" ) ]] <- as.character( bit.depth )
      param.keywords[[ paste0( "$P", i, "E" ) ]] <- "0,0"
      param.keywords[[ paste0( "$P", i, "R" ) ]] <- as.character( asp$expr.data.max )
      param.keywords[[ paste0( "$P", i, "DISPLAY" ) ]] <- "LOG"
      param.keywords[[ paste0( "$P", i, "TYPE" ) ]] <- "Fluorescence"

      # exception for AF.Index
      if ( p.name == "AF.Index" ) {
        param.keywords[[ paste0( "$P", i, "DISPLAY" ) ]] <- "LIN"
        param.keywords[[ paste0( "$P", i, "TYPE" ) ]] <- "AF_Index"
      }

      # assign $PnS (stain) based on flow.control
      f.idx <- match( p.name, flow.control$fluorophore )
      marker <- if ( !is.na( f.idx ) ) flow.control$antigen[ f.idx ] else ""
      param.keywords[[ paste0( "$P", i, "S" ) ]] <- as.character( marker )
    }
  }

  # combine new keywords with original keywords
  new.keywords <- utils::modifyList(
    utils::modifyList( non.param.keys, param.keywords ),
    list(
      "$FIL" = file.name,
      "$PAR" = as.character( n.param ),
      "$UNMIXINGMETHOD" = method,
      "$AUTOSPECTRAL" = as.character( utils::packageVersion( "AutoSpectral" ) ),
      # add AutoSpectralRcpp's version if available
      "$AUTOSPECTRALRCPP" = if ( requireNamespace( "AutoSpectralRcpp", quietly = TRUE ) ) {
        as.character( utils::packageVersion( "AutoSpectralRcpp" ) )
      } else {
        "0"
      }
    )
  )

  # weighting
  if ( !is.null( weights ) ) {
    # add weights to a new keyword in correct format
    weights.str <- paste(
      c(
        length( spectral.channel ),
        spectral.channel,
        formatC(
          weights,
          digits = 8,
          format = "fg"
        )
      ),
      collapse = ","
    )
    new.keywords[[ "$WEIGHTS" ]] <- weights.str
  }

  # spectra
  # TBD switch to using SPILL slot
  fluor.n <- nrow( spectra )
  detector.n <- ncol( spectra )
  fluorophores <- paste0( rownames( spectra ), "-A" )
  vals <- as.vector( t( spectra ) )

  # add spillover/spectra to a new keyword in correct format
  formatted.vals <- formatC( vals, digits = 8, format = "fg", flag = "#" )
  spill.string <- paste(
    c( fluor.n, detector.n, fluorophores, colnames( spectra ), formatted.vals ),
    collapse = ","
  )
  new.keywords[[ "$SPECTRA" ]] <- spill.string
  new.keywords[[ "$FLUOROCHROMES" ]] <- paste( fluorophores, collapse = "," )

  # add AF spectra if used
  if ( !is.null( af.spectra ) ) {
    af.n <- nrow( af.spectra )
    vals <- as.vector( t( af.n ) )
    formatted.vals <- formatC( vals, digits = 8, format = "fg", flag = "#" )
    af.string <- paste(
      c( af.n, detector.n, rownames( af.n ), colnames( af.spectra ), formatted.vals ),
      collapse = ","
    )
    new.keywords[[ "$AUTOFLUORESCENCE" ]] <- af.string
  }

  ### define new FCS file
  # append "-A" to fluorophore and AF channel names
  fluor.orig <- colnames( unmixed.data )
  colnames( unmixed.data ) <-
    ifelse( fluor.orig %in% c( rownames( spectra ), "AF" ),
            paste0( fluor.orig, "-A" ),
            fluor.orig )

  # create the flowFrame for writing the FCS file
  flow.frame <- suppressWarnings( flowCore::flowFrame( unmixed.data ) )
  param.desc <- flowCore::parameters( flow.frame )@data$desc

  # add marker names to description
  for ( i in seq_len( n.param ) ) {
    orig.name <- fluor.orig[ i ]
    # get the marker from flow.control
    f.idx <- match( orig.name, flow.control$fluorophore )
    if ( !is.na( f.idx ) )
      param.desc[ i ] <- as.character( flow.control$antigen[ f.idx ] )
  }

  # write the parameter to the flowFrame
  flowCore::parameters( flow.frame )@data$desc <- param.desc
  keyword( flow.frame ) <- new.keywords

  # save file ---------
  write.FCS(
    flow.frame,
    filename = file.path( output.dir, file.name )
  )

}
