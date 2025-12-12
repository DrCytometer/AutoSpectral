# unmix_fcs.r

#' @title Unmix FCS Data
#'
#' @description
#' This function performs spectral unmixing on FCS data using various methods.
#'
#' @importFrom flowCore read.FCS keyword exprs flowFrame parameters
#' @importFrom flowCore write.FCS parameters<-
#' @importFrom utils packageVersion modifyList
#'
#' @param fcs.file A character string specifying the path to the FCS file.
#' @param spectra A matrix containing the spectral data.
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#' @param flow.control A list containing flow cytometry control parameters.
#' @param method A character string specifying the unmixing method to use.
#' The default is `Automatic`, which uses `AutoSpectral` for AF extraction if
#' af.spectra are provided and automatically selects `OLS` or `WLS` depending
#' on which is normal for the given cytometer in `asp$cytometer`. This means
#' that files from the ID7000, A8 and S8 will be unmixed using `WLS` while
#' others will be unmixed with `OLS`. Any option can be set manually.
#' Manual options are `OLS`, `WLS`, `AutoSpectral`, `Poisson` and `FastPoisson`.
#' Default is `OLS`. `FastPoisson` requires installation of `AutoSpectralRcpp`.
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
#' unmixed FCS file. Default is `NULL`.
#' @param file.suffix A character string to append to the output file name.
#' Default is `NULL`.
#' @param include.raw A logical value indicating whether to include raw
#' expression data in the written FCS file. Default is `FALSE`.
#' @param include.imaging A logical value indicating whether to include imaging
#' parameters in the written FCS file. Default is `FALSE`.
#' @param calculate.error Logical, whether to calculate the RMSE unmixing model
#' accuracy and include it as a keyword in the FCS file.
#' @param use.dist0 Logical, controls whether the selection of the optimal AF
#' signature for each cell is determined by which unmixing brings the fluorophore
#' signals closest to 0 (`use.dist0` = `TRUE`) or by which unmixing minimizes the
#' per-cell residual (`use.dist0` = `FALSE`). Default is `TRUE`. Used for
#' AutoSpectral unmixing.
#' @param divergence.threshold Numeric. Used for `FastPoisson` only.
#' Threshold to trigger reversion towards WLS unmixing when Poisson result
#' diverges for a given point.
#' @param divergence.handling String. How to handle divergent cells from Poisson
#' IRLS. Options are `NonNeg` (non-negativity will be enforced), `WLS` (revert
#' to WLS intial unmix) or `Balance` (WLS and NonNeg will be averaged).
#' Default is `Balance`
#' @param balance.weight Numeric. Weighting to average non-convergent cells.
#' Used for `Balance` option under `divergence.handling`. Default is `0.5`.
#' @param speed Selector for the precision-speed trade-off in AutoSpectral per-cell
#' fluorophore optimization. Options are the default `fast`, which selects the
#' best spectral fit per cell by updating the predicted values for each
#' fluorophore independently without repeating the unnmixing, `medium` which uses
#' a Woodbury-Sherman-Morrison rank-one updating of the unnmixing matrix for
#' better results and a moderate slow-down, or `slow`, which explicitly
#' recomputes the unmixing matrix for each variant for maximum precision. The
#' `fast` method is only one available in the `AutoSpectral` package and will be
#' slow in the pure R implementation. Installation of `AutoSpectralRcpp` is
#' strongly encouraged.
#' @param parallel Logical, default is `TRUE`, which enables parallel processing
#' for per-cell unmixing methods.
#' @param threads Numeric, default is `NULL`, in which case `asp$worker.process.n`
#' will be used. `asp$worker.process.n` is set by default to be one less than the
#' available cores on the machine. Multi-threading is only used if `parallel` is
#' `TRUE`.
#'
#' @return None. The function writes the unmixed FCS data to a file.
#'
#' @export

unmix.fcs <- function( fcs.file, spectra, asp, flow.control,
                       method = "Automatic",
                       weighted = FALSE,
                       weights = NULL,
                       af.spectra = NULL,
                       spectra.variants = NULL,
                       output.dir = NULL,
                       file.suffix = NULL,
                       include.raw = FALSE,
                       include.imaging = FALSE,
                       calculate.error = FALSE,
                       use.dist0 = TRUE,
                       divergence.threshold = 1e4,
                       divergence.handling = "Balance",
                       balance.weight = 0.5,
                       speed = "fast",
                       parallel = TRUE,
                       threads = NULL ) {

  if ( is.null( output.dir ) )
    output.dir <- asp$unmixed.fcs.dir

  if ( is.null( threads ) )
    threads <- asp$worker.process.n

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

  # import fcs, without warnings for fcs 3.2
  fcs.data <- suppressWarnings(
    flowCore::read.FCS( fcs.file, transformation = FALSE,
                        truncate_max_range = FALSE, emptyValue = FALSE )
  )

  fcs.keywords <- flowCore::keyword( fcs.data )
  file.name <- flowCore::keyword( fcs.data, "$FIL" )
  RMSE <- NULL

  # deal with manufacturer peculiarities in writing fcs files
  if ( asp$cytometer %in% c( "ID7000", "Mosaic" ) ) {
    file.name <- sub( "([ _])Raw(\\.fcs$|\\s|$)", paste0("\\1", method, "\\2"),
                      file.name, ignore.case = TRUE )
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
  rm( fcs.data )
  original.param <- colnames( fcs.exprs )

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
    rm( fcs.exprs )

  # remove imaging parameters if desired
  if ( grepl( "Discover", asp$cytometer ) & !include.imaging )
    other.exprs <- other.exprs[ , asp$time.and.scatter ]

  # define weights if needed
  if ( weighted | method == "WLS"| method == "Poisson"| method == "FastPoisson" ) {
    if ( is.null( weights ) ) {
      weights <- pmax( abs( colMeans( spectral.exprs ) ), 1e-6 )
      weights <- 1 / weights
    }
  }

  # apply unmixing using selected method ---------------
  unmixed.data <- switch( method,
                         "OLS" = unmix.ols( spectral.exprs, spectra ),
                         "WLS" = unmix.wls( spectral.exprs, spectra, weights ),
                         "NNLS" = {
                           if (!requireNamespace("RcppML", quietly = TRUE)) {
                             stop(
                               "The package 'RcppML' is required for Non Negative Least Squares.",
                               "Install it with : install.packages('RcppML')",
                               call. = FALSE)
                           }
                           unmix.nnls( spectral.exprs, spectra )
                         },
                         "OLS_pos" = unmix.ols_pos( spectral.exprs, spectra ),
                         "WLS_qtopbot" = unmix.wls_qtopbot( spectral.exprs, spectra, weights ),
                         "AutoSpectral" = {
                           if ( requireNamespace("AutoSpectralRcpp", quietly = TRUE ) &&
                                "unmix.autospectral.rcpp" %in% ls( getNamespace( "AutoSpectralRcpp" ) ) ) {
                             tryCatch(
                               AutoSpectralRcpp::unmix.autospectral.rcpp( raw.data = spectral.exprs,
                                                                          spectra = spectra,
                                                                          af.spectra = af.spectra,
                                                                          spectra.variants = spectra.variants,
                                                                          weighted = weighted,
                                                                          weights = weights,
                                                                          calculate.error = calculate.error,
                                                                          use.dist0 = use.dist0,
                                                                          verbose = asp$verbose,
                                                                          parallel = parallel,
                                                                          threads = threads,
                                                                          speed = speed ),
                               error = function( e ) {
                                 warning( "AutoSpectralRcpp unmixing failed, falling back to standard AutoSpectral: ", e$message )
                                 unmix.autospectral( raw.data = spectral.exprs,
                                                     spectra = spectra,
                                                     af.spectra = af.spectra,
                                                     spectra.variants = spectra.variants,
                                                     weighted = weighted,
                                                     weights = weights,
                                                     calculate.error = calculate.error,
                                                     use.dist0 = use.dist0,
                                                     verbose = asp$verbose )
                               }
                             )
                           } else {
                             warning( "AutoSpectralRcpp not available, falling back to standard AutoSpectral" )
                             unmix.autospectral( raw.data = spectral.exprs,
                                                 spectra = spectra,
                                                 af.spectra = af.spectra,
                                                 spectra.variants = spectra.variants,
                                                 weighted = weighted,
                                                 weights = weights,
                                                 calculate.error = calculate.error,
                                                 use.dist0 = use.dist0,
                                                 verbose = asp$verbose )
                           }
                         },
                         "Poisson" = unmix.poisson( spectral.exprs, spectra, asp, weights ),
                         "FastPoisson" = {
                           if ( requireNamespace("AutoSpectralRcpp", quietly = TRUE ) &&
                               "unmix.poisson.fast" %in% ls( getNamespace( "AutoSpectralRcpp" ) ) ) {
                             tryCatch(
                               AutoSpectralRcpp::unmix.poisson.fast( spectral.exprs,
                                                                     spectra,
                                                                     weights = weights,
                                                                     maxit = asp$rlm.iter.max,
                                                                     tol = 1e-6,
                                                                     n_threads = threads,
                                                                     divergence.threshold = divergence.threshold,
                                                                     divergence.handling = divergence.handling,
                                                                     balance.weight = balance.weight ),
                               error = function( e ) {
                                 warning( "FastPoisson failed, falling back to standard Poisson: ", e$message )
                                 unmix.poisson( spectral.exprs, spectra, asp, weights,
                                                parallel, threads )
                               }
                             )
                           } else {
                             warning( "AutoSpectralRcpp not available, falling back to standard Poisson." )
                             unmix.poisson( spectral.exprs, spectra, asp, weights,
                                            parallel, threads )
                           }
                         },
                         stop( "Unknown method" )
  )

  # calculate model accuracy if desired
  if ( calculate.error & method != "AutoSpectral" ) {
    # for AutoSpectral unmixing, get error directly from the function call
    residual <- rowSums( ( spectral.exprs - ( unmixed.data %*% spectra ) )^2 )
    RMSE <- sqrt( mean( residual ) )
  }

  if ( method == "AutoSpectral" & calculate.error ) {
    RMSE <- unmixed.data$RMSE
    unmixed.data <- unmixed.data$unmixed.data
  }

  # add back raw exprs and others columns as desired
  if ( include.raw ) {
    unmixed.data <- cbind( fcs.exprs, unmixed.data )
  } else {
    unmixed.data <- cbind( other.exprs, unmixed.data )
  }

  rm( spectral.exprs, other.exprs )

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
    setNames( fcs.keywords[ matches], matches )
  } )
  # name the list by parameter name
  names( param.lookup ) <- sapply( pN.keys, function( k ) fcs.keywords[[ k ]])

  # keywords for new parameters
  param.keywords <- list()
  n.param <- ncol( unmixed.data )

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

  # combine
  new.keywords <- modifyList(
    modifyList( non.param.keys, param.keywords ),
    list(
      "$FIL" = file.name,
      "$PAR" = as.character( n.param ),
      "$UNMIXINGMETHOD" = method,
      "$AUTOSPECTRAL" = as.character( packageVersion( "AutoSpectral" ) )
    )
  )

  # RMSE
  if ( calculate.error & !is.null( RMSE ) )
    new.keywords[[ "$RMSE" ]] <- RMSE

  # weighting
  if ( !is.null( weights ) ) {
    weights.str <- paste( c( length( spectral.channel ),
                           spectral.channel,
                           formatC( weights, digits = 8, format = "fg" ) ), collapse = "," )
    new.keywords[[ "$WEIGHTS" ]] <- weights.str
  }

  # spectra
  # TBD switch to using SPILL slot
  fluor.n <- nrow( spectra )
  detector.n <- ncol( spectra )
  vals <- as.vector( t( spectra ) )
  formatted.vals <- formatC( vals, digits = 8, format = "fg", flag = "#" )
  spill.string <- paste(
    c( fluor.n, detector.n, rownames( spectra ), colnames( spectra ), formatted.vals ),
    collapse = ","
  )
  new.keywords[[ "$SPECTRA" ]] <- spill.string
  new.keywords[[ "$FLUOROCHROMES" ]] <- paste( rownames( spectra ), collapse = "," )

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

  # define new FCS file
  flow.frame <- suppressWarnings( flowCore::flowFrame( unmixed.data ) )
  param.desc <- flowCore::parameters( flow.frame )@data$desc
  for ( i in seq_len( n.param ) ) {
    p.name <- colnames( unmixed.data )[ i ]
    # Get the marker from flow.control
    f.idx <- match( p.name, flow.control$fluorophore )
    if ( !is.na( f.idx ) )
      param.desc[ i ] <- as.character( flow.control$antigen[ f.idx ] )
  }
  flowCore::parameters( flow.frame )@data$desc <- param.desc
  keyword( flow.frame ) <- new.keywords

  # save file ---------
  write.FCS( flow.frame, filename = file.path( output.dir, file.name ) )

}
