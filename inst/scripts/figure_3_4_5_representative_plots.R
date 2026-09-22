# figure_3_4_5_representative_plots.R
#
# Representative biplot panels for manuscript Figures 3, 4 and 5.
#
# Figure 3 (BB700, RLM vs positive-negative gating):
#   3A, 3B - create.biplot() of unstained / BB700-stained splenocytes,
#            B10-A vs V7-A, with a two-box positive/negative gate on 3B.
#   3C, 3D - linear-scale regression.biplot() of the same data, with an
#            LM (3C) or RLM (3D) fit line and R^2 annotation.
#   3G     - unmixed A3 Spleen data, PD-1 BV711 vs PDCA-1 BB700, comparing
#            positive-negative-derived vs RLM-derived BB700 spectra.
#
# Figure 4 (BV650, scatter-matching vs positive-negative gating):
#   4A, 4B - produced automatically by clean.controls() via
#            scatter.match.plot() for BUV395 and BV650; see the note below,
#            no code is needed here.
#   4E     - unmixed A3 Spleen data, CD25 BV480 vs XCR1 BV650, comparing
#            unmatched vs scatter-matched BV650 spectra.
#
# Figure 5 (BV605, intrusive-event removal vs positive-negative gating):
#   5A     - create.biplot() of the unstained/AF control's raw data across
#            three detector pairs, illustrating intrusive AF spikes.
#   5G     - unmixed A3 Spleen data, CD24 Alexa Fluor 532 vs IgE BV605,
#            comparing positive-negative-derived vs intrusive-event-removed
#            BV605 spectra.
#
# Panels 3G, 4E and 5G are pre-gated with gate.large.sample() (large.gate =
# TRUE) before unmixing, matching the large gate used elsewhere for this
# sample.

library( AutoSpectral )


# ---- setup --------------------------------------------------------------

asp <- get.autospectral.param()

control.file <- "control_file_cells.csv"
control.dir <- "./Single stained controls"

flow.control.cells <- define.flow.control(
  control.file = control.file,
  control.dir = control.dir,
  asp = asp
)

spleen.file <- "./Fully stained/A3 Spleen_GFP_003_Samples.fcs"

output.dir <- "./figure_3_4_5"
if ( !dir.exists( output.dir ) ) dir.create( output.dir, recursive = TRUE )


# ---- Figure 3A-B: BB700 unstained vs stained, biexponential -------------

bb700.matching.unstained <- flow.control.cells$universal.negative[ "BB700" ]

unstained.data <- flow.control.cells$expr.data[
  which( flow.control.cells$event.sample == bb700.matching.unstained ), ]
bb700.data <- flow.control.cells$expr.data[
  which( flow.control.cells$event.sample == "BB700" ), ]

fig.3a <- create.biplot(
  plot.data = unstained.data,
  x.dim = "B10-A",
  y.dim = "V7-A",
  asp = asp,
  x.lab = "BB700",
  y.lab = "Autofluorescence",
  save = FALSE
)

# Positive/negative gate boxes for 3B. Coordinates are function arguments,
# adjust as needed.
#
# NOTE: the lower bound of the positive gate was supplied as "e34", which is
# not a valid number - 1e4 is used here as a placeholder pending
# confirmation of the intended value.
bb700.gates <- list(
  list( x.range = c( -2e3, 1e3 ), y.range = c( -2e3, 1e6 ), color = "blue" ),
  list( x.range = c( 1e4, 3e5 ), y.range = c( -2e3, 1e6 ), color = "red" )
)

fig.3b <- create.biplot(
  plot.data = bb700.data,
  x.dim = "B10-A",
  y.dim = "V7-A",
  asp = asp,
  x.lab = "BB700",
  y.lab = "Autofluorescence",
  save = FALSE
)

fig.3b <- add.biplot.gate.boxes(
  plot = fig.3b,
  gates = bb700.gates,
  asp = asp
)

ggplot2::ggsave(
  file.path( output.dir, "Figure_3A_BB700_unstained.jpg" ),
  plot = fig.3a, device = ragg::agg_jpeg, width = 5, height = 5
)
ggplot2::ggsave(
  file.path( output.dir, "Figure_3B_BB700_stained_gated.jpg" ),
  plot = fig.3b, device = ragg::agg_jpeg, width = 5, height = 5
)


# ---- Figure 3C-D: BB700 linear regression, LM vs RLM ---------------------

fig.3c <- regression.biplot(
  data = bb700.data,
  x.dim = "B10-A",
  y.dim = "V7-A",
  asp = asp,
  method = "lm",
  x.lab = "BB700",
  y.lab = "Autofluorescence",
  title = "Figure_3C_BB700_lm",
  output.dir = output.dir
)

fig.3d <- regression.biplot(
  data = bb700.data,
  x.dim = "B10-A",
  y.dim = "V7-A",
  asp = asp,
  method = "rlm",
  x.lab = "BB700",
  y.lab = "Autofluorescence",
  title = "Figure_3D_BB700_rlm",
  output.dir = output.dir
)


# ---- Figure 3G: unmixed A3 Spleen, PD-1 BV711 vs PDCA-1 BB700 -----------

raw.spleen.data <- readFCS( spleen.file )

gated.spleen.data <- gate.large.sample(
  flow.data = raw.spleen.data,
  asp = asp,
  samp = "A3_Spleen_GFP_003",
  output.dir = file.path( output.dir, "gate_definitions" )
)

spectra.bb700.posneg <- read.spectra( "Cell pos neg spectra FlowJo renamed.csv" )
spectra.bb700.rlm <- read.spectra( "Cells RLM autospectral spectra.csv" )

unmixed.spleen.posneg <- unmix.ols( gated.spleen.data, spectra.bb700.posneg )
unmixed.spleen.rlm <- unmix.ols( gated.spleen.data, spectra.bb700.rlm )

fig.3g.posneg <- create.biplot(
  plot.data = unmixed.spleen.posneg,
  x.dim = "BV711",
  y.dim = "BB700",
  asp = asp,
  x.lab = "PD-1 BV711",
  y.lab = "PDCA-1 BB700",
  title = "Figure_3G_posneg",
  output.dir = output.dir
)

fig.3g.rlm <- create.biplot(
  plot.data = unmixed.spleen.rlm,
  x.dim = "BV711",
  y.dim = "BB700",
  asp = asp,
  x.lab = "PD-1 BV711",
  y.lab = "PDCA-1 BB700",
  title = "Figure_3G_rlm",
  output.dir = output.dir
)


# ---- Figure 4A-B: produced by clean.controls() --------------------------
#
# These panels are the scatter-matching plots clean.controls() generates
# internally via scatter.match.plot() for BUV395 (4A) and BV650 (4B), using
# the unstained control named in the control file. No separate plotting
# code is needed here - the JPEGs land in asp$figure.scatter.dir.base as:
#
#   BUV395_<asp$scatter.match.plot.filename>
#   BV650_<asp$scatter.match.plot.filename>
#
# generated the next time clean.controls() is run on flow.control.cells.


# ---- Figure 4E: unmixed A3 Spleen, CD25 BV480 vs XCR1 BV650 -------------

spectra.bv650.unmatched <- read.spectra( "Cell pos neg spectra FlowJo renamed.csv" )
spectra.bv650.matched <- read.spectra( "Cells Matching Negative autospectral spectra.csv" )

unmixed.spleen.unmatched <- unmix.ols( gated.spleen.data, spectra.bv650.unmatched )
unmixed.spleen.matched <- unmix.ols( gated.spleen.data, spectra.bv650.matched )

fig.4e.unmatched <- create.biplot(
  plot.data = unmixed.spleen.unmatched,
  x.dim = "BV480",
  y.dim = "BV650",
  asp = asp,
  x.lab = "CD25 BV480",
  y.lab = "XCR1 BV650",
  title = "Figure_4E_unmatched",
  output.dir = output.dir
)

fig.4e.matched <- create.biplot(
  plot.data = unmixed.spleen.matched,
  x.dim = "BV480",
  y.dim = "BV650",
  asp = asp,
  x.lab = "CD25 BV480",
  y.lab = "XCR1 BV650",
  title = "Figure_4E_matched",
  output.dir = output.dir
)


# ---- Figure 5A: intrusive AF spikes, unstained splenocytes --------------

af.unstained.data <- flow.control.cells$expr.data[
  which( flow.control.cells$event.sample == bb700.matching.unstained ), ]

fig.5a.detector.pairs <- list(
  c( "UV6-A", "UV10-A" ),
  c( "B9-A", "V2-A" ),
  c( "R1-A", "YG7-A" )
)

fig.5a.plots <- lapply( fig.5a.detector.pairs, function( pair ) {
  create.biplot(
    plot.data = af.unstained.data,
    x.dim = pair[ 1 ],
    y.dim = pair[ 2 ],
    asp = asp,
    title = paste0( "Figure_5A_", pair[ 1 ], "_vs_", pair[ 2 ] ),
    output.dir = output.dir
  )
} )


# ---- Figure 5G: unmixed A3 Spleen, CD24 Alexa Fluor 532 vs IgE BV605 ----

spectra.bv605.posneg <- read.spectra( "Cell pos neg spectra FlowJo renamed.csv" )
spectra.bv605.cleaned <- read.spectra( "Cells AF Removed autospectral spectra.csv" )

unmixed.spleen.bv605.posneg <- unmix.ols( gated.spleen.data, spectra.bv605.posneg )
unmixed.spleen.bv605.cleaned <- unmix.ols( gated.spleen.data, spectra.bv605.cleaned )

fig.5g.posneg <- create.biplot(
  plot.data = unmixed.spleen.bv605.posneg,
  x.dim = "Alexa Fluor 532",
  y.dim = "BV605",
  asp = asp,
  x.lab = "CD24 Alexa Fluor 532",
  y.lab = "IgE BV605",
  title = "Figure_5G_posneg",
  output.dir = output.dir
)

fig.5g.cleaned <- create.biplot(
  plot.data = unmixed.spleen.bv605.cleaned,
  x.dim = "Alexa Fluor 532",
  y.dim = "BV605",
  asp = asp,
  x.lab = "CD24 Alexa Fluor 532",
  y.lab = "IgE BV605",
  title = "Figure_5G_cleaned",
  output.dir = output.dir
)
