# figure_palettes.R
#
# Colour mappings shared by every figure of the bead-vs-cell manuscript, so a
# given particle type, dye class or cytometer has the same colour everywhere.
# One RColorBrewer palette per categorical variable:
#   particle type  cells = grey; bead types = "Set1"
#   dye class      "Dark2"
#   cytometer      "Set2"
# Colours are assigned by name in the fixed orders below, so they do not change
# when a figure shows only some of the levels.
#   unstained / negative events  black (manuscript.colors()$negative)
# Not covered here (deliberately): laser colours (default.laser.colors in
# figure_helpers.R; they follow wavelength) and sequential/diverging scales
# (viridis for antibody amount, Brewer RdBu / PuOr in the Figure 3 heatmaps).

manuscript.palettes <- list(
  bead      = "Set1",
  class     = "Dark2",
  cytometer = "Set2"
)
manuscript.cell.color <- "grey60"
manuscript.negative.color <- "black"
manuscript.bead.order <- c( "BD", "UltraComp", "UltraCompSpectral", "BioLegend", "SpectraComp" )
manuscript.class.order <- c( "Small", "Brilliant", "Brilliant_tandem", "Real", "Protein", "Protein_tandem" )
manuscript.cytometer.order <- c( "Aurora", "A8", "ID7000", "Xenith" )

#' Named colour vectors for the manuscript
#' @return list(bead, particle (= Cells + bead), class, cytometer, negative)
manuscript.colors <- function() {
  if ( !requireNamespace( "RColorBrewer", quietly = TRUE ) )
    stop( "Install RColorBrewer: install.packages(\"RColorBrewer\")", call. = FALSE )
  pick <- function( levels, palette )
    stats::setNames( RColorBrewer::brewer.pal( max( 3, length( levels ) ), palette )[ seq_along( levels ) ],
                     levels )
  bead <- pick( manuscript.bead.order, manuscript.palettes$bead )
  list(
    bead      = bead,
    particle  = c( Cells = manuscript.cell.color, bead ),
    class     = pick( manuscript.class.order, manuscript.palettes$class ),
    cytometer = pick( manuscript.cytometer.order, manuscript.palettes$cytometer ),
    negative  = manuscript.negative.color
  )
}
