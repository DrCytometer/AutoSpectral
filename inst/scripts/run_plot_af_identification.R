# run_plot_af_identification.R
#
# Example script for generating the "Supplementary Figure 1" autofluorescence
# identification/extraction workflow figure.

library(AutoSpectral)
source( "plot_af_identification.R" )

asp <- get.autospectral.param()

spectra <- read.spectra( "Cells AF Removed autospectral spectra.csv", remove.af = TRUE )

results <- af.identification.plot(
  unstained.sample = file.path( "./SSC", "G2 WT Lung_Samples.fcs" ),
  spectra          = spectra,
  asp              = asp,
  
  n.cells       = 10000L,
  sample.method = "random",
  
  # supply a precomputed af.spectra here to skip get.af.spectra() on repeat
  # runs, e.g. af.spectra = results$af.spectra
  af.spectra           = NULL,
  get.af.spectra.args  = list( refine = FALSE, deduplicate = FALSE ),
  
  panel.b.x.dim = "UV9-A",   # panel.b.y.dim defaults to asp$af.channel ("V7-A" on Aurora)
  panel.c.y.dim = "BUV615",
  
  ssc.channel = "SSC-B-A",
  
  synthetic.fluorophores = c( "BUV395", "BV750", "APC" ),
  
  max.points = 5e4,
  parallel   = TRUE,
  
  panel.width  = 4,
  panel.height = 4,
  
  output.dir = "figure_af_identification",
  save       = TRUE,
  file.type  = "jpg",   # "jpg" | "tiff" | "png" | "pdf" -- use "pdf" for the manuscript
  verbose    = TRUE
)

# results$panels$a .. $h hold the individual panel ggplot objects;
# results$composite is the assembled figure.