# unmix_fcs.r

#' @title Unmix FCS Data
#'
#' @description
#' This function performs spectral unmixing on FCS data using various methods.
#'
#' @importFrom flowCore read.FCS keyword exprs flowFrame parameters pData
#' @importFrom flowCore write.FCS parameters<-
#' @importFrom Biobase AnnotatedDataFrame
#' @importFrom utils packageVersion
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
                       speed = "fast" ){

  if ( is.null( output.dir ) )
    output.dir <- asp$unmixed.fcs.dir

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
  if ( asp$cytometer == "ID7000" ) {
    file.name <- sub( "Raw", method, file.name )

  } else if ( asp$cytometer == "FACSDiscover S8" | asp$cytometer == "FACSDiscover A8" ) {
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

  spectral.channel <- colnames( spectra )
  spectral.exprs <- fcs.exprs[ , spectral.channel, drop = FALSE ]

  other.channels <- setdiff( colnames( fcs.exprs ), spectral.channel )
  other.exprs <- fcs.exprs[ , other.channels, drop = FALSE ]

  if ( !include.raw )
    rm( fcs.exprs )

  # remove imaging parameters if desired
  if ( asp$cytometer == "FACSDiscover S8" | asp$cytometer == "FACSDiscover A8" &
       !include.imaging ) {
    other.exprs <- other.exprs[ , asp$time.and.scatter ]
  }

  # define weights if needed
  # if A8 or S8, pull detector reliability info from FCS file
  # else, use empirical Poisson variance
  if ( weighted | method == "WLS"| method == "Poisson"| method == "FastPoisson" ) {
    if ( is.null( weights ) ) {

      if ( asp$cytometer == "FACSDiscover S8" | asp$cytometer == "FACSDiscover A8" ) {
        qspe <- fcs.keywords[[ "BDSPECTRAL QSPE" ]]
        qspe.values <- strsplit( qspe, ",")[[ 1 ]]
        n.channels <- as.numeric( qspe.values[ 1 ] )
        channel.names <- qspe.values[ 2:( n.channels + 1 ) ]
        weights <- as.numeric( qspe.values[ ( n.channels + 2 ):( n.channels * 2 + 1 ) ] )
        names( weights ) <- channel.names

        if ( all( spectral.channel %in% names( weights ) ) ) {
          weights <- 1 / weights[ spectral.channel ]
        } else {
          # fallback to empirical weighting
          channel.var <- colMeans( spectral.exprs )
          weights <- 1 / ( channel.var + 1e-6 )
        }

      } else {
        # weights are inverse of channel variances (mean if Poisson)
        channel.var <- colMeans( spectral.exprs )
        weights <- 1 / ( channel.var + 1e-6 )
      }
    }
  }

  # apply unmixing using selected method
  unmixed.data <- switch( method,
                         "OLS" = unmix.ols( spectral.exprs, spectra ),
                         "WLS" = unmix.wls( spectral.exprs, spectra, weights ),
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
                                                                          parallel = TRUE,
                                                                          threads = asp$worker.process.n,
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
                                                                     n_threads = asp$worker.process.n,
                                                                     divergence.threshold = divergence.threshold,
                                                                     divergence.handling = divergence.handling,
                                                                     balance.weight = balance.weight ),
                               error = function( e ) {
                                 warning( "FastPoisson failed, falling back to standard Poisson: ", e$message )
                                 unmix.poisson( spectral.exprs, spectra, asp )
                               }
                             )
                           } else {
                             warning( "AutoSpectralRcpp not available, falling back to standard Poisson." )
                             unmix.poisson( spectral.exprs, spectra, asp, weights )
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

  if ( include.raw ) {
    # add back raw exprs and others
    unmixed.data <- cbind( fcs.exprs, unmixed.data )
  } else {
    # add back other columns
    unmixed.data <- cbind( other.exprs, unmixed.data )
  }

  rm( spectral.exprs, other.exprs )

  # fix any NA values (e.g., plate location with S8)
  if ( anyNA( unmixed.data ) )
    unmixed.data[ is.na( unmixed.data ) ] <- 0

  # define new FCS file
  flow.frame <- suppressWarnings( flowCore::flowFrame( unmixed.data ) )

  # set max range values
  params <- pData( parameters( flow.frame ) )
  params$maxRange[ !( params$name %in% c( "Time", "TIME" ) ) ] <- asp$expr.data.max

  # add antigen labels based on match between param name and fluorophore
  fluor.match.idx <- match( params$name, flow.control$fluorophore )

  params$desc[ !is.na( fluor.match.idx ) ] <-
    flow.control$antigen[ fluor.match.idx[ !is.na( fluor.match.idx ) ] ]

  other.match.idx <- match( params$name, other.channels )

  params$desc[ !is.na( other.match.idx ) ] <- NA

  parameters( flow.frame ) <- Biobase::AnnotatedDataFrame( params )

  # update keywords
  keyword( flow.frame ) <- fcs.keywords
  keyword( flow.frame )[[ "$FIL" ]] <- file.name
  keyword( flow.frame )[[ "$UNMIXINGMETHOD" ]] <- method
  keyword( flow.frame )[[ "$AUTOSPECTRAL" ]] <- as.character( packageVersion( "AutoSpectral" ) )

  # add RMSE as a keyword
  if ( calculate.error & !is.null( RMSE ) )
    keyword( flow.frame )[[ "$RMSE" ]] <- RMSE

  # add weighting values as a keyword
  if ( !is.null( weights ) ) {
        weights <- paste( c( length( spectral.channel ),
                             spectral.channel,
            formatC( weights, digits = 8, format = "fg" ) ), collapse = "," )

        keyword( flow.frame )[[ "$WEIGHTS" ]] <- weights
  }

  # write spectral matrix to FCS file
  fluor.n <- nrow( spectra )
  detector.n <- ncol( spectra )
  vals <- as.vector( t( spectra ) )
  formatted.vals <- formatC( vals, digits = 8, format = "fg", flag = "#" )
  spill.string <- paste(
    c( fluor.n, detector.n, rownames( spectra ), colnames( spectra ), formatted.vals ),
    collapse = ","
  )

  keyword( flow.frame )[[ "$SPECTRA" ]] <- spill.string
  keyword( flow.frame )[[ "$FLUOROCHROMES" ]] <- paste( rownames( spectra ), collapse = "," )

  # write parameter names
  for ( i in seq_along( colnames( unmixed.data ) ) ) {
    keyword( flow.frame )[[ paste0( "$P", i, "N" ) ]] <- colnames( unmixed.data )[ i ]
    keyword( flow.frame )[[ paste0( "$P", i, "S" ) ]] <- colnames( unmixed.data )[ i ]
  }

  write.FCS( flow.frame, filename = file.path( output.dir, file.name ) )

}
