# spectral_reference_plot.r

#' @title Spectral Reference Plot
#'
#' @description
#' This function plots a comparison between the user's fluorophore spectrum,
#' from `spectra`, and a known reference spectrum (if available) for the same
#' fluorophore. Quality control is performed using cosine similarity between the
#' two spectral profiles.
#'
#' @importFrom ggplot2 ggplot aes geom_line scale_x_continuous labs scale_linetype_manual
#' @importFrom ggplot2 theme_minimal theme element_text scale_color_manual
#'
#' @param spectra Matrix or dataframe containing spectral data. This
#' should be in format fluorophores x detectors. Row names will be used as the
#' fluorophore names. Column names will be used as the detectors (channels).
#' @param asp The AutoSpectral parameter list. Used to determine which cytometer
#' produced the data.
#' @param qc.threshold Numeric, default `0.98`. The similarity value to trigger
#' a QC failure warning.
#' @param experiment.control.color Color for the line representing the user's
#' fluorophore spectrum from the single-stained reference control. Default is
#' `black`.
#' @param library.reference.color Color for the line representing the library
#' reference standard fluorophore spectrum. Default is `blue`.
#' @param experiment.line.type Line style for the line representing the user's
#' fluorophore spectrum from the single-stained reference control. Default is
#' `solid`.
#' @param library.line.type Line style for the line representing the library
#' reference standard fluorophore spectrum. Default is `dotted`.
#' @param pass.color Color to label similarity values above the `qc.threshold`,
#' i.e., fluorophores passing QC. Default is `darkgreen`.
#' @param fail.color Color to label similarity values below the `qc.threshold`,
#' i.e., fluorophores failing QC. Default is `red`.
#' @param linewidth Width of the line for the spectral traces. Default is `1`.
#' @param plot.dir Directory where the files will be saved.
#' Default is `./figure_spectra`.
#' @param filename Name for the output PDF file. Default is `spectral_qc_report.pdf`.
#'
#' @return None. Plots are saved to a PDF in `plot.dir`.
#'
#' @export

spectral.reference.plot <- function(
  spectra,
  asp,
  qc.threshold = 0.98,
  experiment.control.color = "black",
  library.reference.color = "blue",
  experiment.line.type = "solid",
  library.line.type = "dotted",
  pass.color = "darkgreen",
  fail.color = "red",
  linewidth = 1,
  plot.dir = "./figure_spectra",
  filename = "spectral_qc_report.pdf"
) {

  # which cytometer is being used?
  cyt <- asp$cytometer
  # small switch to consolidate A8/S8 to same reference library
  if ( grepl( "Discover", cyt ) ) cyt <- "Discover"

  # read in the reference spectra for this cytometer
  #database.path <- system.file( "extdata", package = "AutoSpectral" )
  database.path <- getwd()
  ref.file.name <- paste0( cyt, "_spectral_reference_library.csv" )

  # set fluors as row names when reading in the reference spectra
  reference.fluors <- utils::read.csv(
    file.path( database.path, ref.file.name ),
    row.names = 1,
    check.names = FALSE
  )

  # filter to the detectors present (reduce from max)
  detector.names <- colnames( spectra )
  matched.cols <- match( detector.names, colnames( reference.fluors ) )
  reference.fluors <- reference.fluors[ , matched.cols ]

  # re-normalize reference fluor spectra
  reference.fluors <- t( apply(
    reference.fluors,
    1,
    function( x ) x / max( x, na.rm = TRUE )
  ) )

  # match to fluorophores in spectra, excluding AF
  fluorophores <- rownames( spectra )[ rownames( spectra ) != "AF" ]
  ref.plots <- list()

  # summary data frame of cosine similarity QC values
  summary.df <- data.frame(
    Fluorophore = fluorophores,
    Similarity = NA,
    Status = "No Reference Available",
    stringsAsFactors = FALSE
  )

  # for loop, generate all plots
  for ( i in seq_along( fluorophores ) ) {

    f <- fluorophores[ i ]
    simil.value <- NA
    simil.color <- "black"
    lib.spec <- rep( NA, ncol( spectra ) )

    # check if we have a reference for this fluor on this cytometer
    if ( f %in% rownames( reference.fluors ) ) {
      lib.spec <- reference.fluors[f, ]

      simil.value <- cosine.similarity( rbind( spectra[ f, ], lib.spec ) )[ 1, 2 ]

      summary.df$Similarity[ i ] <- round( simil.value, 4 )

      if ( !is.na( simil.value ) && simil.value < qc.threshold ) {
        simil.color <- fail.color
        summary.df$Status[ i ] <- "FAIL"
      } else {
        simil.color <- pass.color
        summary.df$Status[ i ] <- "PASS"
      }
    }

    plot.df <- data.frame(
      detector.n = seq_len( ncol( spectra ) ),
      detector = factor( detector.names, levels = detector.names ),
      experiment = spectra[ f, ],
      library = lib.spec
    )

    qc.text <- if ( is.na( simil.value ) ) {
      "No Reference Available"
    } else {
      paste( "Cosine Similarity:", round( simil.value, 4 ) )
    }

    # create and store plot
    ref.plots[[ f ]] <- ggplot( plot.df, aes( x = detector.n ) ) +
      geom_line(
        aes( y = experiment, color = "Experiment", linetype = "Experiment" ),
        linewidth = linewidth,
        na.rm = TRUE
      ) +
      geom_line(
        aes( y = library, color = "Library Reference", linetype = "Library Reference" ),
        linewidth = linewidth,
        na.rm = TRUE
      ) +
      scale_x_continuous(
        breaks = plot.df$detector.n,
        labels = plot.df$detector
      ) +
      scale_color_manual(
        name = f,
        values = c(
          "Experiment" = experiment.control.color,
          "Library Reference" = library.reference.color
        )
      ) +
      scale_linetype_manual(
        name = f,
        values = c(
          "Experiment" = experiment.line.type,
          "Library Reference" = library.line.type
        )
      ) +
      labs(
        title = f,
        subtitle = qc.text,
        x = "Detector",
        y = "Intensity"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text( angle = 45, hjust = 1 ),
        plot.subtitle = element_text( color = simil.color, face = "bold" ),
        legend.position = "top"
      )
  }

  # console output
  cat( "\n--- Spectral QC Compared to Library References ---\n" )
  print( summary.df )
  cat( "---------------------------\n" )


  # save to pdf
  if ( !dir.exists( plot.dir ) ) dir.create( plot.dir, recursive = TRUE )
  pdf.path <- file.path( plot.dir, filename )

  # A4 landscape PDF for QC summary
  grDevices::pdf( pdf.path, width = 11.7, height = 8.3 )

  # split into multiple pages if there are many fluorophores
  n.fluors <- nrow( summary.df )
  pages <- split( 1:n.fluors, ceiling( seq_along( 1:n.fluors ) / 20 ) )

  for( idx in seq_along( pages ) ) {
    grid::grid.newpage()
    grid::grid.text(
      "Spectral QC Summary Report",
      y = 0.95,
      gp = grid::gpar( fontsize = 16, fontface = "bold" )
    )

    # subheader for multi-page reports
    if( length( pages ) > 1 ) {
      grid::grid.text(
        paste( "Page", idx, "of", length( pages ) ),
        y = 0.92,
        gp = grid::gpar( fontsize = 10 )
      )
    }

    y.pos <- seq( 0.85, by = -0.035, length.out = length( pages[[ idx ]] ) + 1 )
    cols <- c( 0.2, 0.5, 0.8 )
    headers <- colnames( summary.df )

    # draw headers
    for( j in 1:3 ) {
      grid::grid.text(
        headers[ j ],
        x = cols[ j ],
        y = y.pos[ 1 ],
        gp = grid::gpar( fontface = "bold" )
      )
    }

    # draw rows for this chunk
    for( i in seq_along( pages[[ idx ]] ) ) {
      row.idx <- pages[[ idx ]][ i ]
      row.y <- y.pos[ i+1 ]
      row.color <- switch(
        summary.df$Status[ row.idx ],
        "PASS" = pass.color,
        "FAIL" = fail.color,
        "black"
      )

      grid::grid.text(
        summary.df$Fluorophore[ row.idx ],
        x = cols[ 1 ],
        y = row.y,
        gp = grid::gpar( fontsize = 10 )
      )
      grid::grid.text(
        ifelse( is.na( summary.df$Similarity[ row.idx ] ), "N/A", summary.df$Similarity[ row.idx ] ),
        x = cols[ 2 ],
        y = row.y,
        gp = grid::gpar( fontsize = 10 )
      )
      grid::grid.text(
        summary.df$Status[ row.idx ],
        x = cols[ 3 ],
        y = row.y,
        gp = grid::gpar( col = row.color, fontsize = 10, fontface = "bold" )
      )
    }
  }

  # add the plots
  for ( p in ref.plots ) print( p )
  grDevices::dev.off()

  message( paste( "QC Report saved to:", pdf.path ) )

}
