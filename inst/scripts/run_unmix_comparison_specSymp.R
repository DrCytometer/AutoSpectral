# run_unmix_comparison_example.R
#
# Example driver for the two-stage unmixing comparison workflow.
# Assumes the AutoSpectral package (with these two new files added) is loaded.

asp.dir <- "/Users/oliverburton/Bioinformatics/AutoSpectral"
asp.rcpp.dir <- "/Users/oliverburton/Bioinformatics/AutoSpectralRcpp"
devtools::load_all(asp.rcpp.dir)
devtools::load_all(asp.dir)
# or
library(AutoSpectral)
source("setup_unmix_comparison.R")
source("plot_unmix_comparison.R")
source("test_unmix_comparison.R")
source("compare_unmix_folders.R")

# --- one folder per unmixed result set to compare ---------------------------
# names become the x-axis labels on every plot
folders <- c(
  Low_Effort = "./FlowJo Low Effort Unmix/Reference Group_renamed",
  Operator1     = "./SpectroFlo_cells_poor/Unmixed/Reference Group",
  Operator2  = "./SpectroFlo_cells_better/Unmixed/Reference Group",
  Operator2_beads     = "./SpectroFlo_beads/Unmixed/SSC_cells",
  AutoSpectral     = "./AutoSpectral/AutoSpectral_unmixed",
  AutoSpectral_beads  = "./AutoSpectral_beads/unmixed_cell_controls"
)

setup.unmix.comparison(
  folders          = folders,
  fluorophore.list = "./AutoSpectral/cell_control_file.csv",
  output.dir       = "./unmix_comparison_setup"
)

# --- stage 2: compute all six metrics, write plots + results CSV ------------
setup.files <- c(
  Low_Effort = "./unmix_comparison_setup/Low_Effort_unmix_comparison_setup.csv",
  Operator1    = "./unmix_comparison_setup/Operator1_unmix_comparison_setup.csv",
  Operator2    = "./unmix_comparison_setup/Operator2_unmix_comparison_setup.csv",
  Operator2_beads = "./unmix_comparison_setup/Operator2_beads_unmix_comparison_setup.csv",
  AutoSpectral = "./unmix_comparison_setup/AutoSpectral_unmix_comparison_setup.csv",
  AutoSpectral_beads    = "./unmix_comparison_setup/AutoSpectral_beads_unmix_comparison_setup.csv"
)

output <- compare.unmix.folders(
  folders     = folders,
  setup.files = setup.files,
  plot.dir    = "./figure_unmix_comparison",
  output.csv  = "unmix_comparison_results.csv",
  summary.csv = "unmix_comparison_summary.csv"
)

stats.fluor    <- test.unmix.comparison( output$results, output$summary )
stats.channels <- test.unmix.comparison.channels( output$results )

# re-plot the same output with a log axis, larger points, no code re-run
plot.unmix.comparison(
  results    = output$results,
  summary    = output$summary,
  setup.files = setup.files,
  log.scale  = FALSE,
  point.size = 2,
  plot.dir = "comparison_linear",
  plot.per.fluorophore = TRUE
)
plot.unmix.comparison(
  results    = output$results,
  summary    = output$summary,
  setup.files = setup.files,
  log.scale  = TRUE,
  point.size = 2,
  plot.dir = "comparison_log",
  plot.per.fluorophore = TRUE
)

# fluorophore-level stats on the log (ratio) scale
stats.fluor <- test.unmix.comparison(
  output$results, output$summary,
  log.transform = TRUE,
  output.csv        = "unmix_comparison_stats_channels_log.csv",
  detail.csv        = "unmix_comparison_stats_channels_detail_log.csv",
)

stats.fluor.lin <- test.unmix.comparison(
  output$results, output$summary,
  log.transform = FALSE,
  output.csv        = "unmix_comparison_stats_channels_linear.csv",
  detail.csv        = "unmix_comparison_stats_channels_detail_linear.csv",
)

  

















