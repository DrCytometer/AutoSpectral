# af_qc_plot.r

#' @title Autofluorescence QC Plot
#'
#' @description
#' This function plots a comparison between the user's fluorophore spectrum,
#' from `spectra`, and an identified autofluorescence spectrum that may represent
#' contamination. Quality control is performed using cosine similarity between the
#' two spectral profiles.
#'
#' @importFrom ggplot2 ggplot aes geom_line scale_x_continuous labs scale_linetype_manual
#' @importFrom ggplot2 theme_minimal theme element_text scale_color_manual
#'
#' @param af.spectra Matrix or dataframe containing autofluorescence spectral 
#' signatures. This should be in format AFn x detectors. Row names will be used 
#' as the names.
#' @param spectra Matrix or dataframe containing spectral data. This
#' should be in format fluorophores x detectors. Row names will be used as the
#' fluorophore names. Column names will be used as the detectors (channels).
#' @param qc.table Dataframe containing the pre-screened problematic AF spectra
#' and the fluorophores to which they are similar. Format: columes named `AF`,
#' `Fluorophore` and `Similarity` containing, respectively, the names of the AF
#' spectra, the names of the similar fluorophores and the cosine similarity values.
#' @param af.color Color for the line representing the user's autofluorescence 
#' spectrum from the unstained reference control. Default is `black`.
#' @param fluor.color Color for the line representing the similar fluorophore 
#' spectrum. Default is `blue`.
#' @param af.line.type Line style for the line representing the user's
#' AF spectrum. Default is `solid`.
#' @param fluor.line.type Line style for the line representing thefluorophore 
#' spectrum. Default is `dotted`.
#' @param linewidth Width of the line for the spectral traces. Default is `1`.
#' @param plot.dir Directory where the files will be saved.
#' Default is `./figure_autofluorescence`.
#' @param filename Name for the output PDF file. Default is 
#' `autofluorescence_qc_report.pdf`.
#' 
#' @seealso
#' * [qc.af.spectra()]
#' * [get.af.spectra()]
#'
#' @return None. Plots are saved to a PDF in `plot.dir`.
#'
#' @export

af.qc.plot <- function(
    af.spectra,
    spectra,
    qc.table,
    af.color = "black",
    fluor.color = "blue",
    af.line.type = "solid",
    fluor.line.type = "dotted",
    linewidth = 1,
    plot.dir = "./figure_autofluorescence",
    filename = "autofluorescence_qc_report.pdf"
  ) {
  
  detector.names <- colnames( spectra )
  
  af.plots <- list()
  
  # for loop, generate all plots
  for ( i in seq_len( nrow( qc.table ) ) ) {
    
    af <- qc.table$AF[ i ]
    af.spectrum <- af.spectra[ af, ]
    simil.value <- qc.table$Similarity[ i ]
    fluor <- qc.table$Fluorophore[ i ]
    fluor.spectrum <- spectra[ fluor, ]
    
    plot.df <- data.frame(
      detector.n = seq_len( ncol( spectra ) ),
      detector = factor( detector.names, levels = detector.names ),
      AF_Intensity = as.numeric( af.spectrum ),
      Fluorophore_Intensity = as.numeric( fluor.spectrum )
    )
    
    qc.text <- paste( "Cosine Similarity:", round( simil.value, 4 ) )
    title.text <- paste( af, "vs.", fluor )
    
    # create and store plot
    af.plots[[ i ]] <- ggplot( plot.df, aes( x = detector.n ) ) +
      geom_line(
        aes( y = AF_Intensity, color = "AF", linetype = "AF" ),
        linewidth = linewidth,
        na.rm = TRUE
      ) +
      geom_line(
        aes( y = Fluorophore_Intensity, color = "Fluorophore", linetype = "Fluorophore" ),
        linewidth = linewidth,
        na.rm = TRUE
      ) +
      scale_x_continuous(
        breaks = plot.df$detector.n,
        labels = plot.df$detector
      ) +
      scale_color_manual(
        name = title.text,
        values = c(
          "AF" = af.color,
          "Fluorophore" = fluor.color
        )
      ) +
      scale_linetype_manual(
        name = title.text,
        values = c(
          "AF" = af.line.type,
          "Fluorophore" = fluor.line.type
        )
      ) +
      labs(
        title = title.text,
        subtitle = qc.text,
        x = "Detector",
        y = "Intensity"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text( angle = 45, hjust = 1 ),
        plot.subtitle = element_text( color = "black", face = "bold" ),
        legend.position = "top"
      )
  }
  
  # save to pdf
  if ( !dir.exists( plot.dir ) ) dir.create( plot.dir, recursive = TRUE )
  pdf.path <- file.path( plot.dir, filename )
  
  # A4 landscape PDF for QC summary
  grDevices::pdf( pdf.path, width = 11.7, height = 8.3 )
  
  # split into multiple pages if there are many plots
  n.contaminants <- nrow( qc.table )
  pages <- split( 1:n.contaminants, ceiling( seq_along( 1:n.contaminants ) / 20 ) )
  
  for( idx in seq_along( pages ) ) {
    grid::grid.newpage()
    grid::grid.text(
      "AF QC Summary Report",
      y = 0.95,
      gp = grid::gpar( fontsize = 16, fontface = "bold" )
    )
    
    # subheader for multi-page reports
    if ( length( pages ) > 1 ) {
      grid::grid.text(
        paste( "Page", idx, "of", length( pages ) ),
        y = 0.92,
        gp = grid::gpar( fontsize = 10 )
      )
    }
    
    y.pos <- seq( 0.85, by = -0.035, length.out = length( pages[[ idx ]] ) + 1 )
    cols <- c( 0.2, 0.5, 0.8 )
    headers <- colnames( qc.table )
    
    # draw headers
    for ( j in 1:3 ) {
      grid::grid.text(
        headers[ j ],
        x = cols[ j ],
        y = y.pos[ 1 ],
        gp = grid::gpar( fontface = "bold" )
      )
    }
    
    # draw rows for this chunk
    for ( i in seq_along( pages[[ idx ]] ) ) {
      row.idx <- pages[[ idx ]][ i ]
      row.y <- y.pos[ i+1 ]
      
      grid::grid.text(
        qc.table$AF[ row.idx ],
        x = cols[ 1 ],
        y = row.y,
        gp = grid::gpar( fontsize = 10 )
      )
      grid::grid.text(
        qc.table$Fluorophore[ row.idx ],
        x = cols[ 2 ],
        y = row.y,
        gp = grid::gpar( col = "black", fontsize = 10, fontface = "bold" )
      )
      grid::grid.text(
        ifelse( is.na( qc.table$Similarity[ row.idx ] ), "N/A", qc.table$Similarity[ row.idx ] ),
        x = cols[ 3 ],
        y = row.y,
        gp = grid::gpar( fontsize = 10 )
      )
    }
  }
  
  # add the plots
  for ( p in af.plots ) print( p )
  grDevices::dev.off()
  
  message( paste( "AF QC Report saved to:", pdf.path ) )
}