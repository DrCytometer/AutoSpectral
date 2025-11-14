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
#'
#' @return A matrix with the fluorophore spectra.
#'
#' @export

get.fluorophore.spectra <- function( flow.control, asp, use.clean.expr = TRUE,
                                     af.spectra = NULL, title = NULL )
{
  spectra.zero <- rep( 0, flow.control$spectral.channel.n )
  names( spectra.zero ) <- flow.control$spectral.channel

  if ( !is.null( title ) )
    title <- paste( title, asp$spectra.file.name, sep = "_" )
  else if ( use.clean.expr )
    title <- paste( "Clean", asp$spectra.file.name, sep = "_" )
  else
    title <- paste( "Initial", asp$spectra.file.name, sep = "_" )

  # iterate only over non-negative samples
  fluorophore.samples <- flow.control$fluorophore[ ! grepl( "negative",
                                                          flow.control$fluorophore,
                                                          ignore.case = TRUE ) ]

  fluorophore.channels <- flow.control$channel[ ! grepl( "negative",
                                                       flow.control$fluorophore,
                                                       ignore.case = TRUE ) ]
  # check for data
  if ( use.clean.expr ) {
    if ( is.null( flow.control$clean.expr ) )
      stop( "Cleaned control data could not be found.
            Run `clean.controls` first." )
  } else {
    if ( is.null( flow.control$expr.data ) )
      stop( "Control data could not be found in `flow.control`.
            Run `define.flow.control` first." )
  }

  # extract fluophore spectra using original "dirty" data
  if ( !use.clean.expr ) {
    fluorophore.event.samples <- flow.control$event.sample[ flow.control$event.sample %in%
                                                             fluorophore.samples ]
    fluorophore.event.samples <- droplevels( fluorophore.event.samples )

    expr.data <- flow.control$expr.data[ flow.control$event.sample %in%
                                           fluorophore.event.samples, ]

    marker.spectra <- lapply( fluorophore.samples, function( samp ) {

      message( paste( "\033[32m", "Processing", samp, "\033[0m" ) )

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

      # normalize fluor.spectra.coef
      fluor.spectra.coef <- fluor.spectra.coef / max( fluor.spectra.coef )

      fluor.spectra.coef

    } )

    marker.spectra <- do.call( rbind, marker.spectra )
    rownames( marker.spectra ) <- fluorophore.samples

  } else {
    # extract fluophore spectra using "cleaned" data
    fluorophore.event.samples <- flow.control$clean.event.sample[ flow.control$clean.event.sample %in%
                                                                   fluorophore.samples ]
    fluorophore.event.samples <- droplevels( fluorophore.event.samples )

    expr.data <- flow.control$clean.expr[ flow.control$clean.event.sample %in%
                                            fluorophore.event.samples, ]

    marker.spectra <- lapply( fluorophore.samples, function( samp ) {

      message( paste("\033[32m", "Processing", samp, "\033[0m" ) )

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

      # normalize fluor.spectra.coef
      fluor.spectra.coef <- fluor.spectra.coef / max( fluor.spectra.coef )

      fluor.spectra.coef

    } )

    marker.spectra <- do.call( rbind, marker.spectra )
    rownames( marker.spectra ) <- fluorophore.samples
  }

  # plot spectra
  if ( asp$figures ) {
    message( paste("\033[32m", "Plotting figures", "\033[0m" ) )

    fluorophore.spectra.plot <- marker.spectra

    if ( !is.null( af.spectra ) )
      fluorophore.spectra.plot <- rbind( fluorophore.spectra.plot, af.spectra )

    spectral.trace( spectral.matrix = fluorophore.spectra.plot,
                    title = title, plot.dir = asp$figure.spectra.dir,
                    split.lasers = TRUE,
                    asp$figure.spectra.line.size,
                    asp$figure.spectra.point.size )
    spectral.heatmap( fluorophore.spectra.plot, title,
                      plot.dir = asp$figure.spectra.dir )
  }

  if ( !is.null( asp$table.spectra.dir ) )
    write.csv( fluorophore.spectra.plot,
              file = file.path( asp$table.spectra.dir,
                               sprintf( "%s.csv", title ) ) )

  # cosine similarity QC for controls
  if ( asp$figures )
    cosine.similarity.plot( fluorophore.spectra.plot,
                            filename = asp$similarity.heatmap.file.name,
                            title,
                            output.dir = asp$figure.similarity.heatmap.dir,
                            figure.width = asp$figure.similarity.width,
                            figure.height = asp$figure.similarity.height )

  similarity.matrix <- cosine.similarity( fluorophore.spectra.plot )

  unique.similarity <- similarity.matrix * lower.tri( similarity.matrix )

  similarity.idx <- which( unique.similarity > asp$similarity.warning.n,
                           arr.ind = TRUE )

  similarity.error <- nrow( similarity.idx ) > 0

  if ( similarity.error ) {
    fluor1 <- rownames( similarity.matrix )[ similarity.idx[ , 1 ] ]
    fluor2 <- colnames( similarity.matrix )[ similarity.idx[ , 2 ] ]
    similarity.values <- similarity.matrix[ similarity.idx ]

    similarity.qc <- data.frame( Fluor1 = fluor1, Fluor2 = fluor2,
                                Similarity = similarity.values )

    print( similarity.qc )

    message( "\033[31m Similarity over 0.95 detected for one or more pairs of fluorophores.

    Check the table below for problematic combinations.
    If both Fluor1 and Fluor2 are fluorophores,
    manually inspect the controls to confirm they have been prepared correctly.
    Check the fcs_control_table to be sure you have set it up properly.

    If one of the pair is AF, the other likely has minimal signal.
    In this case, run clean.controls and set use.clean.expr to TRUE.
    If you have already done that, manually inspect the control for real signal.
         \033[0m" )
  }

  return( marker.spectra )
}

