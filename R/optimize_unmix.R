# optimize_unmix.r

#' @title Optimize Spectral Unmixing
#'
#' @description
#' Parallel backend for per-cell spectral optimization in R.
#'
#' @param raw.data Numeric matrix (cells x detectors)
#' @param unmixed Numeric matrix (cells x fluors)
#' @param unmix.fun Unmixing function (fast OLS or WLS)
#' @param combined.spectra Numeric matrix (fluors x detectors)
#' @param weights Numeric vector (n detectors)
#' @param pos.thresholds Numeric vector (n fluors)
#' @param af.idx Integer vector (cells), AF spectrum assignments
#' @param af.spectra Numeric matrix (n af_variants x n detectors)
#' @param optimize.fluors Character vector of fluorophores present in variants
#' @param variants List of variant matrices per fluorophore
#' @param delta.list List of delta matrices per fluorophore
#' @param delta.norms List of delta norms per fluorophore
#' @param fluorophores Character vector of fluorophore names
#' @param af.idx.in.spectra Integer, row index of AF in spectra
#' @param asp The AutoSpectral parameter list.
#' @param k Integer, number of variants to test
#' @param cell.weighting Logical, use cell-specific weights
#' @param cell.weight.regularize Logical, regularize cell weights
#' @param nthreads Integer, number of threads
#' @param parallel Logical, whether to use parallel processing
#'
#' @return Unmixed data with cells in rows and fluorophores in columns.
#'
#' @export

optimize.unmix <- function(
    raw.data,
    unmixed,
    unmix.fun,
    combined.spectra,
    weights,
    pos.thresholds,
    af.idx,
    af.spectra,
    optimize.fluors,
    variants,
    delta.list,
    delta.norms,
    fluorophores,
    af.idx.in.spectra,
    asp,
    k = 10L,
    cell.weighting = FALSE,
    cell.weight.regularize = FALSE,
    nthreads = 1L,
    parallel = TRUE,
    refine.af = FALSE
) {

  # this would probably be faster if we moved AF to position 1

  cell.n <- nrow( raw.data )

  # set up parallel backend
  if ( parallel ) {

    result <- create.parallel.lapply( # call from AutoSpectral
      asp,
      # modify exports as needed
      exports = c(
        "raw.data", "unmixed", "unmix.fun",
        "combined.spectra", "weights", "pos.thresholds",
        "af.idx", "af.spectra", "optimize.fluors",
        "variants", "delta.list", "delta.norms", "fluorophores",
        "af.idx.in.spectra", "k", "cell.weighting", "cell.weight.regularize",
        "refine.af"
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

  # loop over each cell, optimizing fluorophore spectra
  unmixed.opt <- tryCatch(
    expr = {
      lapply.function( seq_len( cell.n ), function( cell ) {

        # get cell's data
        cell.raw <- raw.data[ cell, , drop = FALSE ]
        cell.unmixed <- unmixed[ cell, , drop = FALSE ]

        # set weights in a cell-specific manner
        cell.weights <- weights

        if ( cell.weighting ) {
          # use cell-specific weighting (Poisson-like)
          w <- abs( cell.raw )
          w[ w < 1e-6 ] <- 1e-6
          cell.weights <- 1 / w

          if ( cell.weight.regularize ) {
            # regularize weight towards weights for full data
            cell.weights <- ( cell.weights + weights ) * 0.5
          }
        }

        # determine which AF has been selected
        cell.af.idx <- af.idx[ cell ]
        cell.af <- af.spectra[ cell.af.idx, ]

        # calculate delta matrix and norm for this AF spectrum
        af.delta <- sweep( af.spectra, 2, cell.af, "-" )
        af.delta.norm <- sqrt( rowSums( af.delta^2 ) )

        # set baseline spectra
        cell.spectra.final <- combined.spectra
        cell.spectra.final[ af.idx.in.spectra, ] <- cell.af

        # check whether this cell has any fluorophores present
        pos.fluors <- as.numeric( cell.unmixed ) >= pos.thresholds

        # remove absent fluorophores for optimization, unmix
        cell.spectra.curr <- cell.spectra.final[ which( pos.fluors ), , drop = FALSE ]
        cell.unmixed <- unmix.fun( cell.raw, cell.spectra.curr, cell.weights )

        # set baseline unmixed and residuals
        resid <- cell.raw - ( cell.unmixed %*% cell.spectra.curr )
        error.final <- sum( abs( resid ) )


        #############################################
        ### Autofluorescence Optimization Section ###
        #############################################

        # score AF spectra for alignment with residual
        joint.score <- as.numeric( af.delta %*% t( resid ) ) * cell.unmixed[ , "AF" ]
        joint.score <- joint.score / af.delta.norm / sqrt( sum( resid^2 ) )

        # pick top k candidates to test
        k.eff <- min( k, length( joint.score ) )
        topK <- order( joint.score, decreasing = TRUE )[ seq_len( k.eff ) ]

        # test the top k scoring variants
        for ( var in topK ) {
          # supplant the base spectrum with this variant
          cell.spectra.curr[ "AF", ] <- af.spectra[ var, ]

          # reunmix with this variant
          trial.unmix <- unmix.fun( cell.raw, cell.spectra.curr, cell.weights )

          # assess the residual error with this variant
          trial.resid <- cell.raw - ( trial.unmix %*% cell.spectra.curr )
          trial.error <- sum( abs( trial.resid ) )

          # accept change if residual is lower
          if ( trial.error < error.final ) {
            error.final <- trial.error
            cell.spectra.final[ af.idx.in.spectra, ] <- af.spectra[ var, ]
            resid <- trial.resid
          } else {
            # reject if not
            cell.spectra.curr[ "AF", ] <- cell.spectra.final[ af.idx.in.spectra, ]
          }
        }

        # recompute baseline unmix after AF correction
        cell.unmixed <- unmix.fun( cell.raw, cell.spectra.final, cell.weights )

        # reassess fluorophore positivity--this may not be really necessary
        pos.fluors <- as.numeric( cell.unmixed ) >= pos.thresholds

        # if it is necessary, we can again drop out cells with only AF remaining
        if ( !any( pos.fluors[ fluorophores ] ) ) return( cell.unmixed )

        # reset baseline spectra
        cell.spectra.curr <- cell.spectra.final[ which( pos.fluors ), , drop = FALSE ]

        # reset lower dimensional unmixing and targets
        cell.unmixed <- unmix.fun( cell.raw, cell.spectra.curr, cell.weights )
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
              trial.unmix <- unmix.fun( cell.raw, cell.spectra.curr, cell.weights )

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


        #############################################
        ### Autofluorescence Optimization Round 2 ###
        #############################################

        if ( refine.af ) {
          # score AF spectra for alignment with residual
          joint.score <- as.numeric( af.delta %*% t( resid ) ) * cell.unmixed[ , "AF" ]
          joint.score <- joint.score / af.delta.norm / sqrt( sum( resid^2 ) )

          # pick top k candidates to test
          k.eff <- min( k, length( joint.score ) )
          topK <- order( joint.score, decreasing = TRUE )[ seq_len( k.eff ) ]

          # test the top k scoring variants
          for ( var in topK ) {
            # supplant the base spectrum with this variant
            cell.spectra.curr[ "AF", ] <- af.spectra[ var, ]

            # reunmix with this variant
            trial.unmix <- unmix.fun( cell.raw, cell.spectra.curr, cell.weights )

            # assess the residual error with this variant
            trial.resid <- cell.raw - ( trial.unmix %*% cell.spectra.curr )
            trial.error <- sum( abs( trial.resid ) )

            # accept change if residual is lower
            if ( trial.error < error.final ) {
              error.final <- trial.error
              cell.spectra.final[ af.idx.in.spectra, ] <- af.spectra[ var, ]
              resid <- trial.resid
            } else {
              # reject if not
              cell.spectra.curr[ "AF", ] <- cell.spectra.final[ af.idx.in.spectra, ]
            }
          }
        }


        ##############################################################
        ### Final Unmix Using Optimized Spectra (All Fluorophores) ###
        ##############################################################

        cell.unmixed <- unmix.fun( cell.raw, cell.spectra.final, cell.weights )

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
