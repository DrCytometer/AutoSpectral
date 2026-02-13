# get_fluorophore_spectra.r

#' @title Get Fluorophore Spectra
#'
#' @description
#' This function retrieves the fluorophore spectra for flow cytometry data,
#' optionally using cleaned expression data.
#' It also plots and saves the spectra, and performs cosine similarity QC for
#' controls.
#'
#' @param flow.control A list containing flow cytometry control data.
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#' @param use.clean.expr Logical indicating whether to use cleaned expression
#' data, default is `TRUE`
#' @param af.spectra Optional autofluorescence spectra to include.
#' @param title Optional prefix for plot titles, default is `NULL`, which gives
#' "Initial" when `use.clean.expr` is `FALSE` and "Clean" when `use.clean.expr`
#' is `TRUE`.
#' @param figures Logical, default is `TRUE`. Whether to produce plots of the
#' fluorophore spectra and cosine similarity.
#'
#' @return A matrix with the fluorophore spectra.
#'
#' @export

get.fluorophore.spectra <- function(
    flow.control,
    asp,
    use.clean.expr = TRUE,
    af.spectra = NULL,
    title = NULL,
    figures = TRUE
  ) {

  # empty collection vector
  spectra.zero <- rep( 0, flow.control$spectral.channel.n )
  names( spectra.zero ) <- flow.control$spectral.channel

  if ( is.null( title ) & use.clean.expr )
    title <- "Clean"
  else if ( is.null( title ) )
    title <- "Initial"

  # iterate only over non-negative samples
  fluorophore.samples <- flow.control$fluorophore[
    !grepl( "negative", flow.control$fluorophore, ignore.case = TRUE ) ]

  fluorophore.channels <- flow.control$channel[
    !grepl( "negative", flow.control$fluorophore, ignore.case = TRUE ) ]

  # check for data
  if ( use.clean.expr ) {
    if ( is.null( flow.control$clean.expr ) )
      stop( "Cleaned control data could not be found.
            Run `clean.controls()` first." )
  } else {
    if ( is.null( flow.control$expr.data ) )
      stop( "Control data could not be found in `flow.control`.
            Run `define.flow.control()` first." )
  }

  # extract fluophore spectra using original "dirty" data
  if ( !use.clean.expr ) {
    message(
      paste0( "\033[34m",
              "Using original `dirty` expression data to get spectra",
              "\033[0m"
      )
    )

    # identify data corresponding to fluorophores (not negatives)
    fluorophore.event.samples <- flow.control$event.sample[
      flow.control$event.sample %in% fluorophore.samples ]
    fluorophore.event.samples <- droplevels( fluorophore.event.samples )

    expr.data <- flow.control$expr.data[
      flow.control$event.sample %in% fluorophore.event.samples, ]

    # loop over all of these
    marker.spectra <- lapply( fluorophore.samples, function( samp ) {

      message( paste0( "\033[32m", "Processing ", samp, "\033[0m" ) )

      peak.channel <- fluorophore.channels[ fluorophore.samples == samp ]

      peak.channel.expr <- expr.data[
        which( fluorophore.event.samples == samp ),
        peak.channel ]

      fluor.spectra.coef <- spectra.zero

      for ( channel in flow.control$spectral.channel ) {
        if ( channel == peak.channel ) {
          fluor.spectra.coef[ channel ] <- 1.0
        } else {
          channel.expr <- expr.data[
            which( fluorophore.event.samples == samp ), channel ]

          # fit robust linear model
          spectra.model.result <- fit.robust.linear.model(
            peak.channel.expr, channel.expr,
            peak.channel, channel, asp$rlm.iter.max )

          fluor.spectra.coef[ channel ] <- spectra.model.result[ 2 ]
        }
      }

      # normalize fluor.spectra.coef (L-infinity)
      fluor.spectra.coef <- fluor.spectra.coef / max( fluor.spectra.coef )

      fluor.spectra.coef

    } )

    marker.spectra <- do.call( rbind, marker.spectra )
    rownames( marker.spectra ) <- make.unique( fluorophore.samples )

  } else {
    # extract fluorophore spectra using "cleaned" data
    message(
      paste0( "\033[34m",
              "Using cleaned expression data to get spectra",
              "\033[0m"
      )
    )

    # identify data corresponding to fluorophores (not negatives)
    fluorophore.event.samples <- flow.control$clean.event.sample[
      flow.control$clean.event.sample %in% fluorophore.samples ]
    fluorophore.event.samples <- droplevels( fluorophore.event.samples )

    expr.data <- flow.control$clean.expr[
      flow.control$clean.event.sample %in% fluorophore.event.samples, ]

    # loop over all of these
    marker.spectra <- lapply( fluorophore.samples, function( samp ) {

      message( paste0( "\033[32m", "Processing ", samp, "\033[0m" ) )

      peak.channel <- fluorophore.channels[ fluorophore.samples == samp ]

      peak.channel.expr <- expr.data[
        which( fluorophore.event.samples == samp ),
        peak.channel ]

      fluor.spectra.coef <- spectra.zero

      for ( channel in flow.control$spectral.channel ) {
        if ( channel == peak.channel ) {
          fluor.spectra.coef[ channel ] <- 1.0
        } else {
          channel.expr <- expr.data[
            which( fluorophore.event.samples == samp ),
            channel ]

          # fit robust linear model
          spectra.model.result <- fit.robust.linear.model(
            peak.channel.expr, channel.expr,
            peak.channel, channel, asp$rlm.iter.max )

          fluor.spectra.coef[ channel ] <- spectra.model.result[ 2 ]
        }
      }

      # normalize fluor.spectra.coef (L-infinity)
      fluor.spectra.coef <- fluor.spectra.coef / max( fluor.spectra.coef )

      fluor.spectra.coef

    } )

    marker.spectra <- do.call( rbind, marker.spectra )
    rownames( marker.spectra ) <- make.unique( fluorophore.samples )
  }

  # plot spectra
  if ( figures ) {
    message( paste0( "\033[32m", "Plotting figures", "\033[0m" ) )

    fluorophore.spectra.plot <- marker.spectra

    if ( !is.null( af.spectra ) )
      fluorophore.spectra.plot <- rbind( fluorophore.spectra.plot, af.spectra )

    # error handling so plotting never causes the function to abort
    tryCatch(
      expr = {
        spectral.trace(
          spectral.matrix = fluorophore.spectra.plot,
          asp = asp,
          title = paste( title, asp$spectra.file.name, sep = "_" ),
          plot.dir = asp$figure.spectra.dir,
          split.lasers = TRUE,
          figure.spectra.line.size = asp$figure.spectra.line.size,
          figure.spectra.point.size = asp$figure.spectra.point.size
        )

        spectral.heatmap(
          fluorophore.spectra.plot,
          title = paste( title, asp$spectra.file.name, sep = "_" ),
          plot.dir = asp$figure.spectra.dir
        )

        cosine.similarity.plot(
          fluorophore.spectra.plot,
          filename = asp$similarity.heatmap.file.name,
          title,
          output.dir = asp$figure.similarity.heatmap.dir,
          figure.width = asp$figure.similarity.width,
          figure.height = asp$figure.similarity.height
        )

        hotspot.matrix <- calculate.hotspot.matrix( marker.spectra )

        create.heatmap(
          hotspot.matrix,
          title = paste( title, "Hotspot_Matrix", sep = "_" ),
          legend.label = expression( "Hotspot Matrix"^"TM" ),
          triangular = TRUE,
          plot.dir = asp$figure.similarity.heatmap.dir,
          color.palette = "inferno",
          figure.width = asp$figure.similarity.width,
          figure.height = asp$figure.similarity.height
        )

        # calculate OLS unmixing matrix using singular value decomposition
        sv <- svd( t( marker.spectra ) )
        unmixing.matrix <- sv$v %*% ( t( sv$u ) / sv$d )
        colnames( unmixing.matrix ) <- colnames( marker.spectra )
        rownames( unmixing.matrix ) <- rownames( marker.spectra )

        # plot the unmixing matrix as a heatmap
        spectral.heatmap(
          spectra = unmixing.matrix,
          title = paste( title, "unmixing_matrix", sep = "_" ),
          plot.dir = asp$figure.spectra.dir,
          legend.label = "Coefficients",
          color.palette = "mako"
        )
      },
      error = function( e ) {
        message( "Error in plotting fluorophore spectra: ", e$message )
        return( NULL )
      }
    )

  }

  # save the spectra as a CSV file
  if ( !is.null( asp$table.spectra.dir ) ) {
    utils::write.csv(
      marker.spectra,
      file = file.path(
        asp$table.spectra.dir,
        paste0( title, "_", asp$spectra.file.name, ".csv" )
      )
    )
  }

  # cosine similarity QC for controls
  similarity.matrix <- cosine.similarity( marker.spectra )

  unique.similarity <- similarity.matrix * lower.tri( similarity.matrix )

  similarity.idx <- which(
    unique.similarity > asp$similarity.warning.n,
    arr.ind = TRUE
  )

  similarity.error <- nrow( similarity.idx ) > 0

  if ( similarity.error ) {
    fluor1 <- rownames( similarity.matrix )[ similarity.idx[ , 1 ] ]
    fluor2 <- colnames( similarity.matrix )[ similarity.idx[ , 2 ] ]
    similarity.values <- similarity.matrix[ similarity.idx ]

    similarity.qc <- data.frame(
      Fluor1 = fluor1,
      Fluor2 = fluor2,
      Similarity = similarity.values
    )

    print( similarity.qc )

    message(
    "\033[31mSimilarity over 0.95 detected for one or more pairs of fluorophores.

    Check the printed table in the console for problematic combinations.
    If both Fluor1 and Fluor2 are fluorophores,
    manually inspect the controls to confirm they have been prepared correctly.
    Check the fcs_control_file to be sure you have set it up properly.

    If one of the pair is AF, the other likely has minimal signal.
    In this case, run clean.controls and set use.clean.expr to TRUE.
    If you have already done that, manually inspect the control for real signal.
             \033[0m" )
  }

  return( marker.spectra )
}

