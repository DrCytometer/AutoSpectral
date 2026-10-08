# Figure 8 additional unmixing
# take each person's SpectroFlo unmixing, where the controls have been pre-gated
# extract the spectral reference matrix
# apply this to the full single-stained control to achieve comparable plots

# if loading locally:
asp.dir <- "/Users/oliverburton/Bioinformatics/AutoSpectral"
devtools::load_all(asp.dir)
asp <- get.autospectral.param()

control.file <- file.path(
  fig8.omip95.dir, "AutoSpectral", "fcs_control_file.csv"
)
control.dir <- file.path(
  fig8.omip95.dir, "AutoSpectral", "RAW Single stain Spleen"
)
fc <- reload.flow.control(control.dir, control.file, asp)

unstained <- file.path(control.dir, "Spleen_Unstained_Raw.fcs")
apc.e780.ssc <- file.path(control.dir, "APCeF780_FceR1a_SS_Raw.fcs")
apc.e780.fmo <- file.path(
  fig8.omip95.dir, "AutoSpectral", "RAW Full FMO", "FMO", "APCeF780_FceR1a_FMO_Raw.fcs"
)

fig8.omip95.dir <- "../Unmixing comparison/OMIP95 comparison"

fig8.omip95.folders <- c(
  Operator1 = file.path( fig8.omip95.dir, "Jay/OMIP95 NO AF EXTRACTION/OMIP95 NO AF EXTRACTION.Expt" ),
  `Operator1 AF`  = file.path( fig8.omip95.dir, "Jay AF/OMIP95 AF EXTRACTION/OMIP95 AF EXTRACTION.Expt" ),
  `Operator3`  = file.path(
    fig8.omip95.dir,
    "Oliver/OMIP95 Oliver pregated controls/OMIP95 Oliver pregated controls.Expt"
  ),
  `Operator3 AF`  = file.path(
    fig8.omip95.dir,
    "Oliver AF/OMIP95 Oliver pregated controls with AF/OMIP95 Oliver pregated controls with AF.Expt"
  ),
  `Operator4`      = file.path( fig8.omip95.dir, "Sameen/40C_No af extraction/40C_No af extraction.Expt" ),
  `Operator4 AF`      = file.path( fig8.omip95.dir, "Sameen AF/40C_ multiple af extraction_imported efluor 660 positive population/40C_ multiple af extraction_imported efluor 660 positive population.Expt" )
)


# load in spectroflo matrices
output.dir <- "spectroflo_matrices"
spectral.matrices <- lapply(names(fig8.omip95.folders), function(nm) {
  file.name <- paste0(nm, "_spectra.csv")
  read.spectroflo.expt(
    fig8.omip95.folders[[nm]],
    output.dir,
    output.filename = file.name,
    fcs.file = apc.e780.ssc
  )
})
names(spectral.matrices) <- names(fig8.omip95.folders)


for (operator in names(spectral.matrices)) {
  mat <- spectral.matrices[[operator]]
  unmix.fcs(
    unstained,
    mat,
    asp,
    fc,
    method = "OLS",
    file.suffix = operator,
    output.dir = output.dir
  )
  unmix.fcs(
    apc.e780.ssc,
    mat,
    asp,
    fc,
    method = "OLS",
    file.suffix = operator,
    output.dir = output.dir
  )
  unmix.fcs(
    apc.e780.fmo,
    mat,
    asp,
    fc,
    method = "OLS",
    file.suffix = operator,
    output.dir = output.dir
  )
}


