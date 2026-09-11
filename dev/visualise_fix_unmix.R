cell.variants <- readRDS( "./figure_spectral_variants/Spectral_variants_cells.rds" )
result <- fix.my.unmix(
  spectra               = bd.spectra,
  unstained.sample      = file.path( cell.dir, "A10 unstained_010_Cells.fcs" ) ,
  fully.stained.sample  = file.path("./concatenated_fcs", "Concatenated_cells.fcs"),
  flow.control          = fc,
  asp                   = asp,
  variants              = cell.variants,
  max.truncated.events  = 5000,
  figures               = FALSE,
  save                  = FALSE,
  verbose               = TRUE
)


fmuviz.pair.gate.scorecard( result$coefficient.log )
fmuviz.signature.gate.scorecard( result$signature.log )
fmuviz.rejection.summary( result$signature.log )

# To drill into one pair's stratification and mask evolution, supply the
# same unmixed abundances and thresholds fix.my.unmix() would have used
# (e.g. from unmix.ols.fast(raw.data, spectra) and get.spread.thresholds()):
unmixed.comp <- unmix.ols.fast(
  raw.data = readFCS(file.path("./concatenated_fcs", "Concatenated_cells.fcs"))[,colnames(bd.spectra)],
  spectra = bd.spectra
)
threshold.matrix <- get.spread.thresholds(
  unmixed.comp[,1:16],
  thresholds = cell.variants$thresholds,
  spillover.spread = cell.variants$spillover.spread
)

strat <- fmuviz.pair.stratification(
  x.source          = unmixed.comp[ , "PerCP-eFluor 710" ],
  x.target          = unmixed.comp[ , "BV711" ],
  threshold.source  = threshold.matrix[ , "PerCP-eFluor 710" ],
  threshold.target  = threshold.matrix[ , "BV711" ],
  spread.var        = cell.variants$spillover.spread[ "PerCP-eFluor 710", "BV711" ],
  neg.var           = 125^2,
  source.name       = "PerCP-eFluor 710", target.name = "BV711" )
strat$plot.mask
strat$plot.strata

fmuviz.spread.boundary(
  unmixed = unmixed.comp[,1:16], thresholds = cell.variants$thresholds,
  spillover.spread = cell.variants$spillover.spread, spread.kappa = 2,
  source = "PerCP-eFluor 710", target = "BV711" )


