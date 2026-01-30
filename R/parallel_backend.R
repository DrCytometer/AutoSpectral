# parallel_backend.r

#' @title Parallel Backend
#'
#' @description
#' Sets up parallel processing backend using parallel and parallelly to establish
#' PSOCK clusters on Windows. Sets up clusters, defines parLapply-based function.
#'
#' @importFrom parallelly makeClusterPSOCK
#' @importFrom parallel clusterSetRNGStream clusterEvalQ clusterCall stopCluster
#' @importFrom parallel clusterExport parLapply
#' @importFrom RhpcBLASctl blas_set_num_threads omp_set_num_threads
#'
#' @param asp The AutoSpectral parameter list.
#' @param exports The vector of variables and functions to pass to the clusters.
#' @param threads Numeric, number of threads to use for parallel processing.
#' @param export.env The environment containing other functions and global
#' variables potentially needed by the clusters.
#' @param dev.mode Logical, allows testing of function while in development.
#' Default is `FALSE`.
#' @param package.path File.path to the R package files for AutoSpectral to
#' permit loading of the functions via `devtools::load_all()` while in `dev.mode`.
#'
#' @return An lapply function, either based on parLapply if parallel backend
#' initialization was successful or using sequential lapply if not.
#'
#' @export

parallel.backend <- function(
    asp,
    exports,
    threads,
    export.env = parent.frame(),
    dev.mode = FALSE,
    package.path = NULL
  ) {

  # Check if parallelly is available
  if ( !requireNamespace( "parallelly", quietly = TRUE ) ) {
    message( "Package 'parallelly' not found. Install with: install.packages('parallelly')" )
    message( "Falling back to sequential processing." )
    return( list( cl = NULL, lapply = lapply ) )
  }

  # check threads
  if ( is.null( threads ) || !is.numeric( threads ) || threads < 1 ) {
    threads <- 1
  } else {
    threads <- as.integer( threads )
  }

  cl <- tryCatch( {
    parallelly::makeClusterPSOCK(
      threads,
      autoStop = FALSE,
      verbose = FALSE
    )
  }, error = function( e ) {
    message( "Failed to create PSOCK cluster. Falling back to sequential processing." )
    message( "Cause: ", e$message )
    return( NULL )
  } )

  if ( is.null( cl ) ) {
    return( list( cl = NULL, lapply = lapply ) )
  }

  # Set reproducible RNG
  RNGkind( "L'Ecuyer-CMRG" )
  parallel::clusterSetRNGStream( cl, asp$gate.downsample.seed )

  # Load packages and code on all workers
  load.result <- tryCatch( {
    if ( dev.mode && !is.null( package.path )) {
      # Load dependencies
      parallel::clusterEvalQ( cl, {
        library( flowCore )
      } )

      # Get all R files
      r.files <- list.files(
        file.path( package.path, "R" ),
        pattern = "\\.R$",
        full.names = TRUE
        )

      # Source all files into global environment
      parallel::clusterCall( cl, function( files ) {
        for ( f in files ) {
          sys.source( f, envir = .GlobalEnv )
        }
      }, files = r.files )

      # run thread setup
      parallel::clusterEvalQ(cl, {
        try( {
          if ( requireNamespace( "RhpcBLASctl", quietly = TRUE ) ) {
            RhpcBLASctl::blas_set_num_threads( 1 )
            RhpcBLASctl::omp_set_num_threads( 1 )
          } else {
            Sys.setenv(
              OMP_NUM_THREADS = "1",
              OPENBLAS_NUM_THREADS = "1",
              VECLIB_MAXIMUM_THREADS = "1"
            )
          }
        }, silent = TRUE)
      } )

    } else {
      # standard, non-dev mode
      # load AutoSpectral onto clusters
      parallel::clusterEvalQ( cl, {
        library( AutoSpectral )
      } )
      # run thread setup
      parallel::clusterEvalQ(cl, {
        try( {
          if ( requireNamespace( "RhpcBLASctl", quietly = TRUE ) ) {
            RhpcBLASctl::blas_set_num_threads( 1 )
            RhpcBLASctl::omp_set_num_threads( 1 )
          } else {
            Sys.setenv(
              OMP_NUM_THREADS = "1",
              OPENBLAS_NUM_THREADS = "1",
              VECLIB_MAXIMUM_THREADS = "1"
            )
          }
        }, silent = TRUE)
      } )
    }
    TRUE
  }, error = function( e ) {
    try( parallel::stopCluster( cl ), silent = TRUE )
    message( "Failed to load packages on workers. Falling back to sequential processing." )
    message( "Error: ", e$message )
    return( NULL )
  } )

  if ( is.null( load.result ) ) {
    return( list( cl = NULL, lapply = lapply ) )
  }

  # Export variables to workers
  export.result <- tryCatch( {
    if ( length( exports ) > 0 ) {
      parallel::clusterExport( cl, varlist = exports, envir = export.env )
    }
    TRUE
  }, error = function( e ) {
    try( parallel::stopCluster( cl ), silent = TRUE )
    message( "Failed to export variables to workers. Falling back to sequential processing." )
    message( "Error: ", e$message )
    return( NULL )
  } )

  if ( is.null( export.result ) ) {
    return( list( cl = NULL, lapply = lapply,cleanup = NULL ) )
  }

  # Wrapper lapply function
  lapply.fun <- function( x, FUN, ... ) {
    parallel::parLapply( cl, x, FUN, ... )
  }

  cleanup <- function() {
    # Stop cluster gracefully; ignore errors
    try( parallel::stopCluster( cl ), silent = TRUE )
  }

  list( cl = cl, lapply = lapply.fun, cleanup = cleanup )
}
