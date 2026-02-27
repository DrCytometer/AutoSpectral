# get_af_spectra.r

#' @title Get Autofluorescence Spectra
#'
#' @description
#' Extracts autofluorescence spectra from an unstained samples. Intended for use
#' with `unmix.autospectral`. Uses FlowSOM (EmbedSOM) clustering for rapid
#' identification of cells with similar AF profiles.
#'
#' @importFrom FlowSOM SOM
#'
#' @param unstained.sample Path and file name for a unstained sample FCS file.
#' The sample type and processing (protocol) method should match the fully
#' stained samples to which the AF will be applied, ideally.
#' @param asp The AutoSpectral parameter list.
#' @param spectra Spectral signatures of fluorophores, normalized between 0
#' and 1, with fluorophores in rows and detectors in columns.
#' @param som.dim Number of x and y dimensions for the SOM. Default is `10`.
#' @param figures Logical, whether to plot the spectral traces and heatmap for
#' the AF signatures. Default is `TRUE`.
#' @param plot.dir Directory (folder) where the plots will be saved. Default is
#' `NULL`, which inherits from `asp$figure.af.dir`.
#' @param table.dir Directory (folder) where the spectra csv file will be saved.
#' Default is `NULL`, which inherits from `asp$table.af.dir`.
#' @param title Title for the output spectral plots and csv file. Default is
#' `Autofluorescence spectra`.
#' @param verbose Logical, controls messaging. Default is `TRUE`.
#' @param refine Logical, default is `FALSE`. Controls whether to perform a second
#' round of autofluorescence measurement on "problem cells", which are those
#' with the highest spillover, as defined by `problem.quantile`. When `FALSE`,
#' behavior is identical to versions of AutoSpectral prior to 1.0.0. If you are
#' working with samples containing complex autofluorescence, e.g., tissues or
#' tumors, using `refine=TRUE` will improve autofluorescence extraction in the
#' unmixing at the cost of an increase in unmixing time. The increase in time
#' will depend on the method used to assign autofluorescence spectra per cell
#' (residual based assignment is very fast) and whether you have installed
#' `AutoSpectralRcpp`, which will speed up assignment and unmixing.
#' @param problem.quantile Numeric, default `0.99`. The quantile for determining
#' which cells will be considered "problematic" after unmixing with per-cell AF
#' extraction. Cells in the `problem.quantile` or above with respect to total
#' signal in the fluorophore (non-AF) channels after per-cell AF extraction will
#' be used to determine additional autofluorescence spectra, using a second round
#' of clustering and modulation of the previously selected autofluorescence
#' spectra. A value of `0.99` means the top 1% of cells, those farthest from zero,
#' will be selected for further investigation.
#'
#' @return A matrix of autofluorescence spectra.
#'
#' @export
#'
#' @references
#' Van Gassen S et al. (2015). "FlowSOM: Using self-organizing maps for
#' visualization and interpretation of cytometry data." \emph{Cytometry Part A},
#' 87(7), 636-645. \doi{10.1002/cyto.a.22625}
#' Wehrens R, Kruisselbrink J (2018). “Flexible Self-Organizing Maps in kohonen
#' 3.0.” \emph{Journal of Statistical Software}, \emph{87}(7), 1-18.
#' \doi{10.18637/jss.v087.i07}

get.af.spectra <- function(
    unstained.sample,
    asp,
    spectra,
    som.dim = 10,
    figures = TRUE,
    plot.dir = NULL,
    table.dir = NULL,
    title = "Autofluorescence spectra",
    verbose = TRUE,
    refine = FALSE,
    problem.quantile = 0.99
) {

  # set up output folders
  if ( is.null( plot.dir ) )
    plot.dir <- asp$figure.af.dir
  if ( !dir.exists( plot.dir ) )
    dir.create( plot.dir )
  if ( is.null( table.dir ) )
    table.dir <- asp$table.spectra.dir
  if ( !dir.exists( table.dir ) )
    dir.create( table.dir )
  if ( is.null( title ) )
    title <- asp$af.file.name

  # check for AF in spectra, remove if present
  if ( "AF" %in% rownames( spectra ) )
    spectra <- spectra[ rownames( spectra ) != "AF", , drop = FALSE ]

  spectral.channels <- colnames( spectra )

  # import unstained sample
  unstained.ff <- readFCS( unstained.sample, return.keywords = TRUE )

  # get filename for plotting
  file.name <- unstained.ff$keywords[[ "$FIL" ]]

  # extract raw spectral data
  unstained.exprs <- unstained.ff$data[ , spectral.channels ]

  # unmix to get information about how AF propagates into the fluorophore space
  unmixed.no.af <- unmix.ols.fast( unstained.exprs, spectra )

  # use raw and unmixed data as input for better clustering (more info is better)
  cluster.data <- cbind( unstained.exprs, unmixed.no.af )

  # get clusters of AF from unstained
  if ( verbose ) message( "Creating a self-organizing map of the autofluorescence" )

  cell.n <- nrow( cluster.data )
  # reduce som dimensions if few cells are provided
  if ( cell.n < 100 ) {
    stop( paste( "Inadequate cell numbers provided:", cell.n ) )
  } else if ( cell.n < 500 ) {
    # use a minimum of 2 dimensions, or around 3 cells per node
    som.dim <- max( 2, floor( sqrt( cell.n / 3 ) ) )
  }

  # map to SOM (no metaclustering)
  set.seed( 42 )
  map <- FlowSOM::SOM(
    cluster.data,
    xdim = som.dim,
    ydim = som.dim,
    silent = TRUE
  )

  # flip the sign of any negatively vectored AF spectra (occurs on S8, A8)
  af.spectra <- t(
    apply( map$codes[ , spectral.channels ], 1, function( x ) {
      max.x <- ifelse( max( abs( x ) ) > max( x ), min( x ), max( x ) )
      x / max.x } )
  )

  # unlikely, but remove any NAs
  af.spectra <- as.matrix( stats::na.omit( af.spectra ) )

  # add mean as the first row
  mean.af <- colMeans( af.spectra )
  af.spectra <- rbind( mean.af, af.spectra )

  # set names
  rownames( af.spectra ) <- paste0( "AF", 1:nrow( af.spectra ) )

  if ( figures ) {
    if ( verbose ) message( "Plotting autofluorescence spectra" )

    # error handling so plotting never causes the function to abort
    tryCatch(
      expr = {
        # plot the base AF spectra as traces (too many for full refined set)
        spectral.trace(
          spectral.matrix = af.spectra,
          asp = asp,
          title = paste( title, "Autofluorescence spectra" ),
          plot.dir = plot.dir,
          split.lasers = FALSE
        )
        # as a heatmap
        spectral.heatmap( af.spectra, title, plot.dir )
      },
      error = function( e ) {
        message( "Error in plotting AF spectra: ", e$message )
        return( NULL )
      }
    )
  }


  if ( refine ) {
    ###################################################################
    ### Refine Autofluorescence Extraction: identify problem events ###
    ###################################################################

    ### unmix with per-cell AF extraction based on assignments ###

    if ( verbose ) message( "Identifying best-fitting AF: first pass" )
    # allow Rcpp if available
    # assign AF by fluorophore minimization
    af.assignments <- assign.af.fluorophores(
      raw.data = unstained.exprs,
      spectra = spectra,
      af.spectra = af.spectra
    )

    # add dummy AF column (to be filled in)
    af.abundance <- rep( 0, nrow( unstained.exprs ) )
    unmixed <- cbind( af.abundance, unmixed.no.af )
    fluor.idx <- 2:( nrow( spectra ) + 1 )

    # track residuals per cell
    residuals <- matrix(
      0,
      nrow = nrow( unstained.exprs ),
      ncol = ncol( spectra )
    )

    # set base spectra
    combined.spectra <- matrix(
      NA_real_,
      nrow = nrow( spectra ) + 1,
      ncol = ncol( spectra )
    )
    combined.spectra[ fluor.idx, ] <- spectra

    # track projected fluorophore signal
    proj.fluor <- matrix(
      0,
      nrow = nrow( unstained.exprs ),
      ncol = ncol( spectra )
    )

    ### unmix with per-cell AF extraction based on assignments ###
    for ( af in seq_len( nrow( af.spectra ) ) ) {
      # set this AF as the spectrum to use
      combined.spectra[ 1, ] <- af.spectra[ af, ]

      # get the cells using this AF
      cell.idx <- which( af.assignments == af )

      # unmix all used AF spectra
      if ( length( cell.idx ) > 0 ) {
        # unmix with this AF
        unmixed[ cell.idx, ] <- unmix.ols.fast(
          unstained.exprs[ cell.idx, , drop = FALSE ],
          combined.spectra
        )

        # track residuals
        residuals[ cell.idx, ] <- unstained.exprs[ cell.idx, , drop = FALSE ] -
          ( unmixed[ cell.idx, , drop = FALSE ] %*% combined.spectra )

        # project fluorophore signal
        proj.fluor[ cell.idx, ] <- unmixed[ cell.idx, fluor.idx, drop = FALSE ] %*%
          combined.spectra[ fluor.idx, , drop = FALSE ]
      }
    }

    # error is fluorophore signal + residuals
    error <- residuals + proj.fluor


    ### run a second round of clustering on cells that are far from zero ###
    if ( verbose ) message( "Calculating target spectra for problematic cells" )

    # in an unstained sample, any fluorophore signal is error.
    if ( length( fluor.idx ) > 1) {
      error.magnitude <- sqrt( rowSums( unmixed[, fluor.idx ]^2 ) )
    } else {
      error.magnitude <- abs( unmixed[ , fluor.idx ] )
    }

    # filter for problematic cells
    threshold <- stats::quantile( error.magnitude, problem.quantile )
    problem.idx <- which( error.magnitude > threshold )
    problem.cell.n <- length( problem.idx )

    # if we don't have many cells, step down the threshold until we do
    while ( problem.cell.n < 500 ) {
      problem.quantile <- problem.quantile - 0.05

      # if we hit 50%, adjust SOM dimensions and exit
      if ( problem.quantile < 0.5 ) {
        # use a minimum of 2 dimensions, or around 3 cells per node
        som.dim <- max( 2, floor( sqrt( problem.cell.n / 3 ) ) )
        threshold <- stats::quantile( error.magnitude, problem.quantile )
        problem.idx <- which( error.magnitude > threshold )
        problem.cell.n <- length( problem.idx )
        break
      }

      # update for the next iteration check
      threshold <- stats::quantile( error.magnitude, problem.quantile )
      problem.idx <- which( error.magnitude > threshold )
      problem.cell.n <- length( problem.idx )
    }

    if ( problem.cell.n > 10 ) {
      # current AF estimate
      af.abundance <- unmixed[ problem.idx, 1 ]
      af.abundance[ af.abundance == 0 ] <- 1e-6

      # what are the error signatures?
      spill.ratios <- error[ problem.idx, ] / af.abundance

      if ( verbose ) {
        message(
          paste(
            "Clustering",
            length( problem.idx ),
            "cells with highest spillover error"
          )
        )
      }

      # set safe SOM dimensions
      som.dim <- max( 2, floor( sqrt( problem.cell.n / 3 ) ) )

      # cluster only the problematic data
      set.seed( 42 )
      map.error <- FlowSOM::SOM(
        spill.ratios,
        xdim = som.dim,
        ydim = som.dim,
        silent = TRUE
      )

      # for each cluster, we find which base AFs were present and update them
      cluster.ids <- unique( map.error$mapping[ , 1 ] )

      modulated.list <- lapply( cluster.ids, function( cl) {
        cl.sub.idx <- which( map.error$mapping[ , 1] == cl )
        global.idx <- problem.idx[ cl.sub.idx ]

        # get the median correction pattern for this cluster
        median.ratio <- apply(
          spill.ratios[ cl.sub.idx, , drop = FALSE ],
          2,
          stats::median
        )

        # identify all unique base AFs that fell into this error cluster
        contributing.af.ids <- unique( af.assignments[ global.idx ] )

        # create a new version of each contributing AF using this cluster's ratio
        new.specs <- lapply( contributing.af.ids, function( id ) {
          base.spec <- af.spectra[ id, ]
          updated <- base.spec * ( 1 + median.ratio )
          return( updated / max( updated ) ) # Re-normalize
        } )

        return( do.call( rbind, new.specs ) )
      } )

      modulated.af.spectra <- do.call( rbind, modulated.list )
      af.spectra <- rbind( af.spectra, modulated.af.spectra )

      # unlikely, but check for NAs
      af.spectra <- as.matrix( stats::na.omit( af.spectra ) )

      # add names
      rownames( af.spectra ) <- paste0( "AF", 1:nrow( af.spectra ) )

      if ( figures ) {
        ### before and after AF extraction plotting ###
        if ( verbose ) message( "Identifying best-fitting AF: second pass" )

        # reassign AF using full set of AF spectra
        af.assignments <- assign.af.fluorophores(
          raw.data = unstained.exprs,
          spectra = spectra,
          af.spectra = af.spectra
        )

        # create a copy of the data for AF extraction with the full set of AF spectra
        unmixed.second <- unmixed

        # unmix per assigned AF
        for ( af in seq_len( nrow( af.spectra ) ) ) {
          # set this AF as the spectrum to use
          combined.spectra[ 1, ] <- af.spectra[ af, ]

          # get the cells using this AF
          cell.idx <- which( af.assignments == af )

          if ( length( cell.idx ) > 0 ) {
            # unmix with this AF
            unmixed.second[ cell.idx, ] <- unmix.ols.fast(
              unstained.exprs[ cell.idx, , drop = FALSE ],
              combined.spectra
            )
          }
        }

        # plot the changes for the worst pair of channels
        if ( verbose ) message( "Plotting impact of AF extraction" )

        # which channels are the worst affected by AF when we don't extract it?
        channel.sd <- apply( unmixed.no.af, 2, stats::sd )
        worst.channels <- order( channel.sd, decreasing = TRUE )[ 1:2 ]
        worst.channels <- colnames( unmixed.no.af )[ worst.channels ]


        # error handling so plotting never causes the function to abort
        tryCatch(
          expr = {
            # plot unmixed data before and after AF extraction
            create.biplot(
              unmixed.no.af,
              x.dim = worst.channels[ 1 ],
              y.dim = worst.channels[ 2 ],
              asp = asp,
              title = paste( file.name, "_", title, "_No_AF_Extraction" ),
              output.dir = plot.dir
            )
            create.biplot(
              unmixed,
              x.dim = worst.channels[ 1 ],
              y.dim = worst.channels[ 2 ],
              asp = asp,
              title = paste0( file.name, "_", title, "_PerCell_AF_Extraction_First_Pass" ),
              output.dir = plot.dir
            )
            create.biplot(
              unmixed.second,
              x.dim = worst.channels[ 1 ],
              y.dim = worst.channels[ 2 ],
              asp = asp,
              title = paste0( file.name, "_", title, "_PerCell_AF_Extraction_Second_Pass" ),
              output.dir = plot.dir
            )
          },
          error = function( e ) {
            message( "Error in plotting AF extraction: ", e$message )
            return( NULL )
          }
        )
      }

    } else {
      message( "Insufficient error-prone cells found. Skipping modulation." )
    }
  }

  # save as CSV file for later use
  if ( is.null( title ) )
    af.file.name <- paste0( file.name, "_", asp$af.file.name, ".csv" )
  else
    af.file.name <- paste0( file.name, "_", title, ".csv" )

  utils::write.csv( af.spectra, file = file.path( table.dir, af.file.name ) )

  # plotting as heatmaps and signatures
  if ( figures ) {
    ### plotting of spectra ###
    if ( verbose ) message( "Plotting autofluorescence variation" )

    # error handling so plotting never causes the function to abort
    tryCatch(
      expr = {
        # plot the full set as a spectral variation plot
        spectral.variant.plot(
          af.spectra,
          mean.af,
          title = paste( title, "Autofluorescence variation" ),
          save = TRUE,
          plot.dir = plot.dir
        )
      },
      error = function( e ) e
    )
  }

  return( af.spectra )
}
