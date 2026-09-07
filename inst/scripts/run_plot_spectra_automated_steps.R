# run_plot_spectra_automated_steps.R
#
# Example script for generating the manuscript figure illustrating the
# get.spectra.automated() pipeline steps for one or more single-stained
# controls.

asp.dir <- "/Users/oliverburton/Bioinformatics/AutoSpectral"
devtools::load_all(asp.dir)
asp.rcpp.dir <- "/Users/oliverburton/Bioinformatics/AutoSpectralRcpp"
devtools::load_all(asp.rcpp.dir)

asp <- get.autospectral.param()

control.dir      <- "./raw_data/SSC_cells_beads"
create.control.file(control.dir, asp)
control.def.file <- "workflow_control_file.csv"
control.table <- read.csv(control.def.file)
control.table$fluorophore
# need to add in bead controls--create a folder just for this
# need new control file too

# One or more fluorophore names, exactly as they appear in the `fluorophore`
# column of the control definition file.
fluors.to.plot <- unique(control.table$fluorophore)[1:5]

results <- plot.spectra.automated.steps(
  control.dir       = control.dir,
  control.def.file  = control.def.file,
  asp               = asp,
  fluorophores      = fluors.to.plot,
  clean.positive.point.size = 5,
  output.dir = "figure_automated_pipeline",
  n.spectral = 50
)

std.results <- plot.spectra.standard.workflow(
  control.dir       = control.dir,
  control.def.file  = control.def.file,
  asp               = asp,
  fluorophores      = fluors.to.plot,
  clean.positive.point.size = 5,
  output.dir = "figure_standard_workflow"
)

legacy.results <- plot.spectra.legacy.steps(
  control.dir       = control.dir,
  control.def.file  = control.def.file,
  asp               = asp,
  fluorophores      = fluors.to.plot,
  clean.positive.point.size = 5,
  output.dir = "figure_legacy_workflow"
)


