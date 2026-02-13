# optimize_unmix.r

#' @title Optimize Spectral Unmixing
#'
#' @description
#' Parallel backend for per-cell spectral optimization in R.
#'
#' @param raw.data Numeric matrix (cells x detectors)
#' @param unmixed Numeric matrix (cells x fluors)
#' @param spectra Numeric matrix (fluors x detectors)
#' @param pos.thresholds Numeric vector (n fluors)
#' @param optimize.fluors Character vector of fluorophores present in variants
#' @param variants List of variant matrices per fluorophore
#' @param delta.list List of delta matrices per fluorophore
#' @param delta.norms List of delta norms per fluorophore
#' @param fluorophores Character vector of fluorophore names
#' @param asp The AutoSpectral parameter list.
#' @param k Integer, number of variants to test
#' @param nthreads Integer, number of threads
#' @param parallel Logical, whether to use parallel processing
#'
#' @return Unmixed data with cells in rows and fluorophores in columns.
#'
#' @export

optimize.unmix <- function(
    raw.data,
    unmixed,
    spectra,
    pos.thresholds,
    optimize.fluors,
    variants,
    delta.list,
    delta.norms,
    fluorophores,
    asp,
    k = 10L,
    nthreads = 1L,
    parallel = TRUE
) {

  # how many cells are there?
  cell.n <- nrow( raw.data )

  # set up parallel backend
  if ( parallel ) {

    result <- create.parallel.lapply( # call from AutoSpectral
      asp,
      # modify exports as needed
      exports = c(
        "raw.data", "unmixed", "unmix.ols.fast",
        "spectra", "pos.thresholds", "optimize.fluors",
        "variants", "delta.list", "delta.norms", "fluorophores",
        "k"
      ),
      parallel = TRUE,
      threads = nthreads,
      export.env = environment()
    )
    lapply.function <- result$lapply
  } else {
    lapply.function <- lapply
    result <- list( cleanup = NULL )
  }

  # should add error catch, return original unmixing for that cell
  # loop over each cell, optimizing fluorophore spectra
  unmixed.opt <- tryCatch(
    expr = {
      lapply.function( seq_len( cell.n ), function( cell ) {

        # get cell's data
        cell.raw <- raw.data[ cell, , drop = FALSE ]
        cell.unmixed <- unmixed[ cell, , drop = FALSE ]

        # set baseline spectra
        cell.spectra.final <- spectra

        # check whether this cell has any fluorophores present
        pos.fluors <- as.numeric( cell.unmixed ) >= pos.thresholds
        # early exit if no fluorophores present
        if ( !any( pos.fluors[ fluorophores ] ) ) return( cell.unmixed )

        # remove absent fluorophores for optimization, unmix
        cell.spectra.curr <- cell.spectra.final[ which( pos.fluors ), , drop = FALSE ]
        cell.unmixed <- unmix.ols.fast( cell.raw, cell.spectra.curr )

        # set baseline unmixed and residuals
        resid <- cell.raw - ( cell.unmixed %*% cell.spectra.curr )
        error.final <- sum( abs( resid ) )

        ########################################
        ### Fluorophore Optimization Section ###
        ########################################

        # restrict optimization to fluorophores we have variants for
        ### TBD: use indexing rather than %in%
        fluors.to.sort <- optimize.fluors[
          optimize.fluors %in% names( pos.fluors )[ pos.fluors ] ]

        if ( length( fluors.to.sort ) > 0 ) {
          # sort by abundance to optimize brightest fluors first (error is proportional to signal)
          fluor.order <- sort( cell.unmixed[ , fluors.to.sort ], decreasing = TRUE )

          for ( fl in names( fluor.order ) ) {
            fl.variants <- variants[[ fl ]]
            delta.fl <- delta.list[[ fl ]]
            delta.norm  <- delta.norms[[ fl ]]

            # score variants
            # Note: here we score based on the L2 error (sum of the squared residual)
            # Later, we evaluate changes based on the L1 error (sum of the absolute residual)
            # This mismatch may be helpful in avoiding overfitting, but should be investigated further.
            # Using the L2 error is closer to the optimization performed by least squares,
            # but in all my testing, it performs a bit worse. This is not something I understand.
            joint.score <- as.numeric( delta.fl %*% t( resid ) ) * cell.unmixed[ , fl ]
            joint.score <- joint.score / delta.norm / sqrt( sum( resid^2 ) )

            # select k variants up to the max we have available
            k.eff <- min( k, length( joint.score ) )
            topK <- order( joint.score, decreasing = TRUE )[ seq_len( k.eff ) ]

            # test the top k scoring variants
            for ( var in topK ) {
              # supplant the base spectrum with this variant
              cell.spectra.curr[ fl, ] <- fl.variants[ var, ]

              # reunmix with this variant
              trial.unmix <- unmix.ols.fast( cell.raw, cell.spectra.curr )

              # assess the residual error with this variant
              trial.resid <- cell.raw - ( trial.unmix %*% cell.spectra.curr )
              trial.error <- sum( abs( trial.resid ) )

              # accept change if residual is lower
              if ( trial.error < error.final ) {
                error.final <- trial.error
                cell.spectra.final[ fl, ] <- cell.spectra.curr[ fl, ]
                resid <- trial.resid
              } else {
                # reject if not
                cell.spectra.curr[ fl, ] <- cell.spectra.final[ fl, ]
              }
            }
          }
        }


        ##############################################################
        ### Final Unmix Using Optimized Spectra (All Fluorophores) ###
        ##############################################################

        cell.unmixed <- unmix.ols.fast( cell.raw, cell.spectra.final )

        return( cell.unmixed )
      } )
    },
    finally = {
      if ( !is.null( result$cleanup ) ) result$cleanup()
    }
  )

  # combine data into a matrix
  return( do.call( rbind, unmixed.opt ) )
}
