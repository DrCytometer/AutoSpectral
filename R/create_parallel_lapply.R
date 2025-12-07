# create_parallel_lapply.r

#' @title Create Parallel Lapply
#'
#' @description
#' Sets up parallel processing lapply function for Windows, Mac OS or Linux.
#'
#' @importFrom parallelly makeClusterPSOCK availableCores
#' @importFrom parallel clusterSetRNGStream clusterEvalQ clusterCall stopCluster
#' @importFrom parallel clusterExport parLapply mclapply
#'
#' @param asp The AutoSpectral parameter list.
#' @param exports The vector of variables and functions to pass to the clusters.
#' @param parallel Logical, controls whether parallel processing is used. Default
#' is `TRUE`.
#' @param threads Numeric, number of threads to use for parallel processing.
#' Default is `NULL` which will revert to `asp$worker.process.n` if
#' `parallel=TRUE`.
#' @param export.env The environment containing other functions and global
#' variables potentially needed by the clusters. Default is `parent.frame()`.
#' @param dev.mode Logical, allows testing of function while in development.
#' Default is `FALSE`.
#' @param package.path File.path to the R package files for AutoSpectral to
#' permit loading of the functions via `devtools::load_all()` while in `dev.mode`.
#' Default is `NULL`.
#'
#' @return An lapply function, either based on `parLapply` for Windows or `mcLapply`
#' on Mac OS and Linux. If parallel backend initialization fails, sequential
#' `lapply` is returned.



create.parallel.lapply <- function( asp,
                                    exports,
                                    parallel = TRUE,
                                    threads = NULL,
                                    export.env = parent.frame(),
                                    dev.mode = FALSE,
                                    package.path = NULL ) {
  lapply.function <- NULL
  cleanup <- NULL
  if ( is.null( threads ) ) threads <- asp$worker.process.n
  if ( threads == 0 ) threads <- parallelly::availableCores()
  os <- Sys.info()[[ "sysname" ]]

  if ( dev.mode && is.null( package.path ) ) {
    package.path <- getwd()
    if ( !file.exists( file.path( package.path, "DESCRIPTION" ) ) ) {
      message( "dev.mode=TRUE but no DESCRIPTION file found in directory.")
      message( "Falling back to sequential processing." )
      parallel <- FALSE
    }
  }

  if ( parallel ) {
    if ( os == "Windows" ) {
      for ( var in exports ) {
        if ( !exists( var, envir = export.env ))
          message( "WARNING: Variable '", var, "' not found in export environment!" )
      }

      backend <- tryCatch(
        parallel.backend( asp, exports, threads, export.env = export.env,
                          dev.mode = dev.mode, package.path = package.path ),
        error = function( e ) {
          message( "Windows parallel backend failed. Falling back to sequential processing." )
          message( "Cause: ", e$message )
          NULL
        }
      )

      if ( is.null( backend ) || is.null( backend$cl ) ) {
        message( "Using sequential processing." )
        lapply.function <- lapply
      } else {
        lapply.function <- backend$lapply
        cleanup <- function() parallel::stopCluster( backend$cl )
      }

    } else {
      message( "Using mcLapply with ", threads, " cores..." )
      lapply.function <- tryCatch( {
        function( x, FUN, ... ) {
          RNGkind( "L'Ecuyer-CMRG" )
          set.seed( asp$gate.downsample.seed )
          parallel::mclapply( x, FUN, ...,
                              mc.cores = threads,
                              mc.preschedule = FALSE )
        }
      }, error = function( e ) {
        message( "mclapply failed. Falling back to sequential processing." )
        message( "Cause: ", e$message )
        lapply
      } )
    }
  } else {
    message( "Parallel processing disabled. Using sequential processing." )
    set.seed( asp$gate.downsample.seed )
    lapply.function <- lapply
  }

  return( list( lapply = lapply.function, cleanup = cleanup ) )
}
