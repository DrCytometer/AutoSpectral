# spectral_trace.r

#' @title Plot Fluorophore Spectra Traces
#'
#' @description
#' This function plots the fluorophore spectra as normalized traces, optionally
#' splitting by excitation lasers, and saves the plots as JPEG files.
#'
#' @importFrom ggplot2 ggplot aes geom_path geom_point labs theme_minimal
#' @importFrom ggplot2 element_text facet_wrap ggsave theme scale_color_viridis_d
#' @importFrom tidyr pivot_longer
#' @importFrom utils read.csv
#' @importFrom stats setNames
#'
#' @param spectral.matrix Matrix or dataframe containing spectral data. This
#' should be in format fluorophores x detectors. Row names will be used as the
#' fluorophore names. Column names will be used as the detectors (channels).
#' @param asp The AutoSpectral parameter list defined using
#' `get.autospectral.param`.
#' @param title Title for the plot. Default is `Fluorophore_Spectra`
#' @param plot.dir Directory to save the plot files. Default is `NULL`, in
#' which case the current working directory will be used.
#' @param split.lasers Logical indicating whether to create a second plot split
#' by excitation lasers. Default is `TRUE`.
#' @param figure.spectra.line.size Numeric. Defines the width of the trace lines
#' on the plot. Default is `1`.
#' @param figure.spectra.point.size Numeric. Defines the size of the points
#' on the plot (one per detector). Default is `1`.
#' @param color.palette Optional character string defining the viridis color
#' palette to be used for the fluorophore traces. Default is `NULL`, in which
#' case default R Brewer colors will be assigned automatically. Options are the
#' viridis color options: `magma`, `inferno`, `plasma`, `viridis`, `cividis`,
#' `rocket`, `mako` and `turbo`.
#' @param show.legend Logical. If `TRUE`, figure legend will be included.
#' @param plot.width Optional numeric to manually set the plot width. Default
#' is `NULL`.
#' @param plot.height Optional numeric to manually set the plot width. Default
#' is `NULL`.
#' @param save Logical, if `TRUE`, saves a JPEG file to the `plot.dir`.
#' Otherwise, the plot will simply be created in the Viewer.
#'
#' @return Saves the plot(s) as JPEG files in the specified directory.
#' @export

spectral.trace <- function( spectral.matrix,
                            asp,
                            title = "Fluorophore_Spectra",
                            plot.dir = NULL,
                            split.lasers = TRUE,
                            figure.spectra.line.size = 1,
                            figure.spectra.point.size = 1,
                            color.palette = NULL,
                            show.legend = TRUE,
                            plot.width = NULL,
                            plot.height = NULL,
                            save = TRUE ) {


  # get excitation laser based on peak emission detector
  peak.detectors <- colnames( spectral.matrix )[ max.col( spectral.matrix ) ]

  # convert to data.frame for plotting
  fluor.spectra.plotting <- data.frame( spectral.matrix, check.names = FALSE )
  fluor.spectra.plotting$Fluorophore <- rownames( fluor.spectra.plotting )

  if ( is.null( plot.dir ) )
    plot.dir <- getwd()

  data.path <- system.file( "extdata", "cytometer_database.csv",
                            package = "AutoSpectral" )

  # can only do splitting if we know the laser
  if ( data.path == "" ) {
    warning( "cytometer_database.csv.csv not found in extdata." )
    split.lasers <- FALSE
  }

  cytometer.database <- read.csv( data.path )
  cytometer.database[ cytometer.database == "" ] <- NA

  if ( asp$cytometer == "Aurora" ) {
    if ( asp$cytometer.version == "NL" ) {
      detectors <- setNames(
        cytometer.database$NorthernLights, cytometer.database$NorthernLights_laser )
    } else {
      detectors <- setNames(
        cytometer.database$Aurora, cytometer.database$Aurora_laser )
    }
  } else if ( asp$cytometer == "ID7000" ) {
    detectors <- setNames(
      cytometer.database$ID7000, cytometer.database$ID7000_laser )
  } else if ( asp$cytometer == "FACSDiscover A8" ) {
    detectors <- setNames(
      cytometer.database$Discover, cytometer.database$Discover_laser )
  } else if ( asp$cytometer == "FACSDiscover S8" ) {
    detectors <- setNames(
      cytometer.database$Discover, cytometer.database$Discover_laser )
  } else if ( asp$cytometer == "Opteon" ) {
    detectors <- setNames(
      cytometer.database$Opteon, cytometer.database$Opteon_laser )
  } else if ( asp$cytometer == "Mosaic" ) {
    detectors <- setNames(
      cytometer.database$Mosaic, cytometer.database$Mosaic_laser )
  } else if ( asp$cytometer == "Xenith" ) {
    detectors <- setNames(
      cytometer.database$Xenith, cytometer.database$Xenith_laser )
  } else if ( asp$cytometer == "Symphony" ) {
    detectors <- setNames(
      cytometer.database$A5SE, cytometer.database$A5SE_laser )
  } else {
    warning( "Unsupported cytometer" )
    split.lasers <- FALSE
  }

  detectors <- detectors[ !is.na( detectors ) ]

  laser.order <- unique( names( detectors ) )

  laser.idx <- match( peak.detectors, detectors )

  fluor.spectra.plotting$Laser <- names( detectors )[ laser.idx ]
  fluor.spectra.plotting$Laser[ is.na( fluor.spectra.plotting$Laser ) ] <- "Violet"
  fluor.spectra.plotting$Laser <- factor(
    fluor.spectra.plotting$Laser,
    levels = laser.order
  )

  if ( is.null( plot.width ) )
    plot.width <- max( ( ( ncol( fluor.spectra.plotting ) - 1 ) / 64 * 12 ), 3 )

  if ( is.null( plot.height ) )
    plot.height <- 5 + round( nrow( fluor.spectra.plotting ) / 8, 0 )

  fluor.spectra.long <- tidyr::pivot_longer(
    fluor.spectra.plotting,
    -c( Fluorophore, Laser ),
    names_to = "Detector",
    values_to = "Intensity"
  )

  fluor.spectra.long$Detector <-  factor(
    fluor.spectra.long$Detector,
    levels = unique( fluor.spectra.long$Detector ),
    ordered = TRUE
  )

  spectra.plot <- ggplot(
    fluor.spectra.long,
    aes(
      x = Detector,
      y = Intensity,
      group = Fluorophore,
      color = Fluorophore
      )
    ) +
    geom_path( linewidth = figure.spectra.line.size ) +
    geom_point( size = figure.spectra.point.size ) +
    labs( title = title,
         x = "Detector",
         y = "Normalized Intensity" ) +
    theme_minimal() +
    theme( axis.text.x = element_text( angle = 45, hjust = 1 )  ) +
    theme( legend.position = "bottom" )

  if ( !is.null( color.palette ) )
    spectra.plot <- spectra.plot +
      scale_color_viridis_d( option = color.palette )

  if ( !show.legend )
    spectra.plot <- spectra.plot + theme( legend.position = "none" )

  if ( save ) {
    ggsave( file.path( plot.dir, sprintf( "%s.jpg", title )),
            spectra.plot,
            width = plot.width, height = plot.height,
            limitsize = FALSE )
  } else {
    return( spectra.plot )
  }

  if ( split.lasers ) {
    # get number of lasers used
    laser.n <- length( unique( fluor.spectra.plotting$Laser ) )

    plot.height <- ( plot.height - 1 ) * laser.n

    spectra.plot.split <- ggplot(
      fluor.spectra.long,
      aes(
        x = Detector,
        y = Intensity,
        group = Fluorophore,
        color = Fluorophore
        )
      ) +
      geom_path( linewidth = figure.spectra.line.size ) +
      geom_point( size = figure.spectra.point.size ) +
      facet_wrap( ~ Laser, nrow = laser.n ) +
      labs(
        title = title,
        x = "Detector",
        y = "Normalized Intensity"
      ) +
      theme_minimal() +
      theme( axis.text.x = element_text( angle = 45, hjust = 1 )  ) +
      theme( legend.position = "bottom" )

    if ( !is.null( color.palette ) )
      spectra.plot.split <- spectra.plot.split +
      scale_color_viridis_d( option = color.palette )

    if ( !show.legend )
      spectra.plot.split <- spectra.plot.split + theme( legend.position = "none" )

    if ( save ) {
      ggsave(
        file.path( plot.dir, sprintf( "%s by laser.jpg", title ) ),
        spectra.plot.split,
        width = plot.width,
        height = plot.height,
        limitsize = FALSE
      )
    } else {
      return( spectra.plot.split )
    }
  }
}
