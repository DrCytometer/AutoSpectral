# create_parallel_lapply.r

#' @title Create Parallel Lapply
#'
#' @description
#' Sets up parallel processing lapply function for Windows, Mac OS or Linux.
#'
#' @importFrom parallelly makeClusterPSOCK availableCores
#' @importFrom parallel clusterSetRNGStream clusterEvalQ clusterCall stopCluster
#' @importFrom parallel clusterExport parLapply mclapply
#' @importFrom RhpcBLASctl blas_set_num_threads omp_set_num_threads
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
#' @param allow.mclapply.mac Logical, if `TRUE` permits `mclapply()` forking on
#' Mac OS. Default `FALSE` forces PSOCK cluster use on Mac to prevent
#' multithreaded Accelerate BLAS from crashing parallels when matrix ops are used.
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
                                    package.path = NULL,
                                    allow.mclapply.mac = FALSE ) {

  if ( is.null( threads ) ) threads <- asp$worker.process.n
  if ( is.null( threads ) || !is.numeric( threads ) ) threads <- 1
  threads <- as.integer( threads )
  if ( threads == 0 ) {
    threads <- parallelly::availableCores()
  }

  os <- Sys.info()[[ "sysname" ]]
  lapply.function <- NULL
  cleanup <- NULL

  if ( !parallel ) {
    message( "Parallel processing disabled. Using sequential processing." )
    set.seed( asp$gate.downsample.seed )
    lapply.function <- lapply

    return(
      list(
        lapply = lapply.function,
        cleanup = cleanup
      )
    )
  }

  # Helper: wrap user FUN so that inside mclapply each child enforces single-thread BLAS
  make.mclapply.wrapper <- function( mc.cores ) {
    function( x, FUN, ... ) {
      wrapper.FUN <- function( ... ) {
        # try to restrict BLAS threads inside forked child
        try( {
          if ( requireNamespace( "RhpcBLASctl", quietly = TRUE ) ) {
            RhpcBLASctl::blas_set_num_threads( 1 )
            RhpcBLASctl::omp_set_num_threads( 1 )
          } else {
            Sys.setenv( OMP_NUM_THREADS = "1",
                        OPENBLAS_NUM_THREADS = "1",
                        VECLIB_MAXIMUM_THREADS = "1" )
          }
        }, silent = TRUE )

        FUN( ... )
      }
      RNGkind( "L'Ecuyer-CMRG" )
      set.seed( asp$gate.downsample.seed )
      parallel::mclapply( x, wrapper.FUN, mc.cores = mc.cores, mc.preschedule = FALSE, ... )
    }
  }

  if ( dev.mode && is.null( package.path ) ) {
    package.path <- getwd()
    if ( !file.exists( file.path( package.path, "DESCRIPTION" ) ) ) {
      message( "dev.mode=TRUE but no DESCRIPTION file found in directory.")
      message( "Falling back to sequential processing." )
      parallel <- FALSE
    }
  }

  # Platform-specific behavior
  if ( os == "Windows" || ( os == "Darwin" && !allow.mclapply.mac ) ) {
    # Windows and macOS default: use PSOCK cluster backend for safety with BLAS
    for ( var in exports ) {
      if ( !exists( var, envir = export.env ) ) {
        message( "WARNING: Variable '", var, "' not found in export environment!" )
      }
    }

    backend <- tryCatch( {
      parallel.backend( asp, exports = exports, threads = threads,
                        export.env = export.env, dev.mode = dev.mode,
                        package.path = package.path )
    }, error = function( e ) {
      message( "PSOCK parallel backend failed. Falling back to sequential processing." )
      message( "Cause: ", e$message )
      NULL
    } )

    if ( is.null( backend ) || is.null( backend$lapply ) ) {
      message( "Using sequential processing." )
      lapply.function <- lapply
    } else {
      lapply.function <- backend$lapply
      cleanup <- backend$cleanup
    }
  } else {
    # Assume Linux or other UNIX-like: allow mclapply but enforce single-thread BLAS per child
    message( "Using mcLapply with ", threads, " cores (fork-safe)." )
    lapply.function <- make.mclapply.wrapper( mc.cores = threads )
    # No cluster to cleanup for mclapply
    cleanup <- NULL
  }

  return( list( lapply = lapply.function, cleanup = cleanup ) )
}
