# run_compensaid_comparison_example.R
#
# Example driver: runs CompensAID on the fully-stained sample(s) produced by
# each unmixing method/operator already compared in
# run_unmix_comparison_example.R, and summarises flagged marker
# combinations per operator as a barchart. Assumes the AutoSpectral package
# (with compare_compensaid_folders.R and plot_compensaid_comparison.R
# added) is loaded, and that setup.unmix.comparison() has already been run
# so the setup CSVs referenced below exist.

asp.dir <- "/Users/oliverburton/Bioinformatics/AutoSpectral"
asp.rcpp.dir <- "/Users/oliverburton/Bioinformatics/AutoSpectralRcpp"
devtools::load_all(asp.rcpp.dir)
devtools::load_all(asp.dir)
source("~/Bioinformatics/AutoSpectral/inst/scripts/plot_compensaid_comparison.R")
source("~/Bioinformatics/AutoSpectral/inst/scripts/compare_compensaid_folders.R")
source("/Volumes/T7 Shield/AutoSpectral figure data/Revision/Unmixing comparison/SpecSymp analysis/plot_compensaid_comparison.R")
source("/Volumes/T7 Shield/AutoSpectral figure data/Revision/Unmixing comparison/SpecSymp analysis/plot_compensaid_dotplot_row.R")
source("/Volumes/T7 Shield/AutoSpectral figure data/Revision/Unmixing comparison/SpecSymp analysis/compare_compensaid_folders.R")
source("/Volumes/T7 Shield/AutoSpectral figure data/Revision/Unmixing comparison/SpecSymp analysis/test_compensaid_comparison.R")
source("~/Bioinformatics/AutoSpectral/inst/scripts/plot_compensaid_dotplot_row.R")
source("~/Bioinformatics/AutoSpectral/inst/scripts/compare_compensaid_folders.R")
source("~/Bioinformatics/AutoSpectral/inst/scripts/compensaid_patch.R")
source("~/Bioinformatics/AutoSpectral/inst/scripts/compensaid_patch_2.R")
source("~/Bioinformatics/AutoSpectral/inst/scripts/compensaid_patch_3.R")
# these patches reduce unnecessary calls, cutting overhead while preserving hash-
# identical results in testing
# computational time is cut by ~2.5x

# --- reuse the same folders/setup.files as the single-stain comparison -----
folders <- c(
  `FlowJo Unmixing Wizard` = "./FlowJo Low Effort Unmix/Samples_renamed/",
  Operator1     = "./SpectroFlo_cells_poor/Unmixed/Samples",
  Operator2  = "./SpectroFlo_cells_better/Unmixed/Samples",
  `Operator2 beads`     = "./SpectroFlo_beads/Unmixed/Samples",
  AutoSpectral       = "./AutoSpectral/Unmixed_samples",
  `AutoSpectral beads`  = "./AutoSpectral_beads/Unmixed_samples"
)

setup.files <- c(
  `FlowJo Unmixing Wizard` = "./unmix_comparison_setup/FlowJo Unmixing Wizard_unmix_comparison_setup.csv",
  Operator1    = "./unmix_comparison_setup/Operator1_unmix_comparison_setup.csv",
  Operator2    = "./unmix_comparison_setup/Operator2_unmix_comparison_setup.csv",
  `Operator2 beads` = "./unmix_comparison_setup/Operator2_beads_unmix_comparison_setup.csv",
  AutoSpectral = "./unmix_comparison_setup/AutoSpectral_unmix_comparison_setup.csv",
  `AutoSpectral beads`    = "./unmix_comparison_setup/AutoSpectral_beads_unmix_comparison_setup.csv"
)

# --- fully-stained sample(s) per folder, paths relative to each folder -----
# one or more real multi-color specimens per operator/method -- NOT the
# single-stain controls used above. Edit these filenames to match what is
# actually in each folder; a folder can list more than one file if you want
# to score multiple replicate specimens per operator.
sample.files <- list(
  `FlowJo Unmixing Wizard` = c(
    Rep1 = "A1 Spleen_WT_001_Samples.fcs",
    Rep2 = "A2 Spleen_WT_002_Samples.fcs",
    Rep3 = "A3 Spleen_GFP_003_Samples.fcs",
    Rep4 = "A4 Spleen_GFP_004_Samples.fcs"
  ),
  Operator1 = c(
    Rep1 = "A1 Spleen_WT_001.fcs",
    Rep2 = "A2 Spleen_WT_002.fcs",
    Rep3 = "A3 Spleen_GFP_003.fcs",
    Rep4 = "A4 Spleen_GFP_004.fcs"
  ),
  Operator2 = c(
    Rep1 = "A1 Spleen_WT_001.fcs",
    Rep2 = "A2 Spleen_WT_002.fcs",
    Rep3 = "A3 Spleen_GFP_003.fcs",
    Rep4 = "A4 Spleen_GFP_004.fcs"
  ),
  `Operator2 beads`= c(
    Rep1 = "/Spleen/A1 Spleen_WT_001.fcs",
    Rep2 = "/Spleen/A2 Spleen_WT_002.fcs",
    Rep3 = "/Spleen/A3 Spleen_GFP_003.fcs",
    Rep4 = "/Spleen/A4 Spleen_GFP_004.fcs"
  ),
  AutoSpectral = c(
    Rep1 = "A1 Spleen_WT_001_Samples AutoSpectral.fcs",
    Rep2 = "A2 Spleen_WT_002_Samples AutoSpectral.fcs",
    Rep3 = "A3 Spleen_GFP_003_Samples AutoSpectral.fcs",
    Rep4 = "A4 Spleen_GFP_004_Samples AutoSpectral.fcs"
  ),
  `AutoSpectral beads` = c(
    Rep1 = "A1 Spleen_WT_001_Samples AutoSpectral.fcs",
    Rep2 = "A2 Spleen_WT_002_Samples AutoSpectral.fcs",
    Rep3 = "A3 Spleen_GFP_003_Samples AutoSpectral.fcs",
    Rep4 = "A4 Spleen_GFP_004_Samples AutoSpectral.fcs"
  )
)

# --- run CompensAID per folder, write results/summary CSVs + barchart ------
output <- compare.compensaid.folders(
  folders      = folders,
  setup.files  = setup.files,
  sample.files = sample.files,
  plot.dir               = "./figure_compensaid_comparison",
  output.csv             = "compensaid_comparison_results.csv",
  replicate.summary.csv  = "compensaid_comparison_replicate_summary.csv",
  pair.summary.csv       = "compensaid_comparison_pair_summary.csv",
  summary.csv            = "compensaid_comparison_summary.csv",
  normalize.plot         = TRUE
)

saveRDS(output, file = "CompensAID_spleen_results.rds")

# output$results: one row per (folder, sample.file, primary.fluorophore,
#                  secondary.fluorophore) with ssi, flagged, severity
# output$summary: one row per folder with n.tested, n.flagged, frac.flagged,
#                  and one count column per severity band (Severe/Moderate/Mild)

# re-plot with the folder order fixed to match the single-stain comparison
# above, plus the normalized (percent-flagged) panel, without re-running
# CompensAID
plot.compensaid.comparison(
  #results       = output$results,
  summary       = output$summary,
  folder.levels = names( folders ),
  normalize     = TRUE,
  plot.dir      = "./figure_compensaid_comparison"
)

# marker combinations CompensAID actually flagged, worst case across
# samples where a folder has more than one, for manual follow-up
flagged.pairs <- output$results[ output$results$flagged, ]
flagged.pairs <- flagged.pairs[
  order( flagged.pairs$folder, flagged.pairs$ssi ),
]
utils::write.csv(
  flagged.pairs, "compensaid_flagged_pairs_spleen.csv", row.names = FALSE
)
