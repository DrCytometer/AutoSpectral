# run_plot_spectra_standard_workflow.R
#
# Example script for generating the "standard" (manual-gating, e.g.
# SpectroFlo-style) workflow comparison figure, for the same fluorophore(s)
# as run_plot_spectra_automated_steps.R.

library(AutoSpectral)
source( "plot_spectra_standard_workflow.R" )

asp <- get.autospectral.param()

control.dir      <- "./Controls"
control.def.file <- "fcs_control_file.csv"

fluors.to.plot <- c( "BUV805", "PE", "PE-Cy7" )

results <- plot.spectra.standard.workflow(
  control.dir       = control.dir,
  control.def.file  = control.def.file,
  asp               = asp,
  fluorophores      = fluors.to.plot,
  
  # panel A: octagon gate
  octagon.width.factor = 3,           # bigger = larger octagon around the
  # density peak (lymphocyte region)
  gate.color            = "darkgoldenrod1",
  density.palette       = "rainbow",  # viridis option, or any other value to
  # fall back to asp$density.palette.base.color
  
  # panel B: brightest / negative event selection
  n.bright.events        = 2000L,     # "positive" population size
  negative.quantile      = 0.25,      # upper quantile bound of the negative gate
  negative.quantile.min  = 0.01,      # lower quantile bound (excludes extreme dim tail)
  x.min.quantile         = 0.005,     # x-axis lower limit for the panel B histogram
  selection.fill.color   = "steelblue",
  selection.line.color   = "black",
  negative.bracket.color = "#377EB8",
  positive.bracket.color = "#E41A1C",
  
  # panel C: negative (black) vs positive (cosine-coloured) biplots
  negative.point.color = "black",
  event.point.size     = 3,
  
  # panel D: background subtraction
  subtraction.trace.colors = c(
    "Negative"                      = "grey40",
    "Positive (pre-subtraction)"    = "#1B9E77",
    "Final (background-subtracted)" = "#D95F02"
  ),
  
  max.points   = 5e4,
  panel.width  = 4,
  panel.height = 4,
  
  output.dir = "figure_standard_workflow",
  save       = TRUE,
  file.type  = "jpg",   # "jpg" | "tiff" | "png" | "pdf" -- use "pdf" for the manuscript
  verbose    = TRUE
)

# Each element of `results` (one per fluorophore) holds the individual panel
# ggplot objects plus the assembled `composite`, e.g.:
#   results[["BUV395"]]$composite
#   results[["BUV395"]]$cosine.panel