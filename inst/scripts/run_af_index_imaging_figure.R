# run_af_index_imaging_figure.R
#
# Manuscript figure: correlation between the AutoSpectral AF Index and
# BD FACSDiscover S8 imaging parameters.
#
# Stage A: legacy spectra extraction, AF extraction, AutoSpectral unmixing,
#          retaining imaging + scatter parameters alongside the unmixed data.
# Stage B: FSC/SSC density gate (landmark gating) -> cells only.
# Stage C: FlowSOM + consensus metaclustering and tSNE on imaging parameters.
# Stage D: figure panels.
# Stage E: statistical association tests.


asp.dir      <- "/Users/oliverburton/Bioinformatics/AutoSpectral"
asp.rcpp.dir <- "/Users/oliverburton/Bioinformatics/AutoSpectralRcpp"

if ( requireNamespace( "devtools", quietly = TRUE ) ) {
  try( devtools::load_all( asp.rcpp.dir ), silent = TRUE )
  devtools::load_all( asp.dir )
}
library( ggplot2 )
source( "af_index_imaging_figure.R" )

# ---------------------------------------------------------------------------
# Stage A: legacy spectra extraction, AF extraction, unmixing
# ---------------------------------------------------------------------------

asp <- get.autospectral.param( cytometer = "s8", figures = TRUE )

control.dir      <- "./SSC"
create.control.file(control.dir, asp)
control.def.file <- "fcs_control_file.csv"

check.control.file( control.dir, control.def.file, asp, strict = TRUE )

flow.control <- define.flow.control(
  control.dir, control.def.file, asp,
  gate = TRUE, gating.system = "density"
)
flow.control <- clean.controls( flow.control, asp )
spectra      <- get.fluorophore.spectra( flow.control, asp )

unstained.sample <- "./SSC/Unstained control_brain.fcs"
af.spectra <- get.af.spectra( unstained.sample, asp, spectra, som.dim = 10 )

stained.fcs.file <- "./20240305_AH_S8/Sample_brain 1.1/Sample_brain 1.1.fcs"
raw.data.full <- readFCS( stained.fcs.file )

scatter.params      <- asp$default.scatter.parameter
imaging.params.raw  <- get.imaging.parameter.columns( colnames( raw.data.full ) )

raw.data.spectral <- raw.data.full[ , colnames( spectra ) ]

unmixed <- unmix.autospectral.rcpp(
  raw.data   = raw.data.spectral,
  spectra    = spectra,
  af.spectra = af.spectra,
  verbose    = TRUE,
  parallel   = TRUE,
  threads    = 4L
)

full.data <- cbind(
  unmixed,
  raw.data.full[ , c( scatter.params, imaging.params.raw ), drop = FALSE ]
)

# ---------------------------------------------------------------------------
# Stage B: FSC/SSC density gate (large.gate) -> cells only
# ---------------------------------------------------------------------------

if ( !dir.exists( asp$figure.gate.dir ) ) dir.create( asp$figure.gate.dir, recursive = TRUE )

gate.boundary <- define.gate.landmarks(
  control.def.file, control.dir, asp, gate.name = "Myeloid"
)

gated.data <- apply.gate(
  full.data, gate.boundary, scatter.param = scatter.params, asp = asp
)

# ---------------------------------------------------------------------------
# Stage C: FlowSOM + consensus metaclustering and UMAP on imaging parameters
# ---------------------------------------------------------------------------

imaging.params <- get.imaging.parameter.columns( colnames( gated.data ) )

cluster.result <- cluster.imaging.parameters(
  imaging.data      = gated.data[ , imaging.params ],
  asp               = asp,
  cluster.on        = "embedding",   # "embedding" (default) or "imaging" -- toggle here
  z.score           = TRUE,
  use.pca           = TRUE,
  pca.var.explained = 0.9,
  som.xdim          = 6,
  som.ydim          = 6,
  n.metaclusters    = 10,
  max.embed.events  = 4e4,
  umap.min.dist     = 0.1,
  threads           = 4L
)

# absolute row indices into gated.data for whichever events cluster.id /
# metacluster.id cover -- works regardless of cluster.on
cluster.event.idx <- cluster.event.index( cluster.result )

imaging.for.cluster <- gated.data[ cluster.event.idx, imaging.params ]

af.reorder <- reorder.af.index.by.similarity( af.spectra )
af.index.for.cluster.raw <- gated.data[ cluster.event.idx, "AF Index" ]
af.index.for.cluster     <- af.index.to.continuous( af.index.for.cluster.raw, af.reorder )

# ---------------------------------------------------------------------------
# Stage D: figure panels
# ---------------------------------------------------------------------------

# optional: rename metaclusters for the legend, e.g.
# cluster.labels <- c( "1" = "Granular / high SSC-Imaging", "2" = "Round, low texture" )
cluster.labels <- NULL

# canonical, full set of cluster ids -- pass to every cluster-colored plot
# below so colors are guaranteed concordant across panels
cluster.levels <- sort( unique( cluster.result$metacluster.id ) )

p.embed.clusters <- plot.embedding.clusters(
  embed.coords   = cluster.result$embed.coords,
  cluster.id     = align.to.embedding( cluster.result, cluster.result$metacluster.id ),
  cluster.labels = cluster.labels,
  cluster.levels = cluster.levels,
  asp            = asp
)

star.markers <- if ( length( imaging.params ) > 8 ) {
  param.sd <- apply( imaging.for.cluster, 2, stats::sd, na.rm = TRUE )
  names( sort( param.sd, decreasing = TRUE ) )[ 1:8 ]
} else imaging.params

p.flowsom.stars <- plot.flowsom.stars(
  cluster.result$fsom,
  markers        = star.markers,
  cluster.levels = cluster.levels
)

p.heatmap <- plot.cluster.imaging.heatmap(
  imaging.data   = imaging.for.cluster,
  cluster.id     = cluster.result$metacluster.id,
  imaging.params = imaging.params,
  af.index       = af.index.for.cluster,
  cluster.labels = cluster.labels,
  asp            = asp
)

p.embed.af.index <- plot.embedding.af.index(
  embed.coords = cluster.result$embed.coords,
  af.index     = align.to.embedding( cluster.result, af.index.for.cluster ),
  asp          = asp
)

af.violin.result <- plot.af.index.violin(
  af.index       = af.index.for.cluster,
  cluster.id     = cluster.result$metacluster.id,
  cluster.labels = cluster.labels,
  cluster.levels = cluster.levels,
  asp            = asp
)
p.af.violin <- af.violin.result$plot

# ---------------------------------------------------------------------------
# Stage E: statistics
# ---------------------------------------------------------------------------

chisq.result <- test.af.index.cluster.chisq(
  cluster.id = cluster.result$metacluster.id,
  af.index   = af.index.for.cluster.raw
)

message( sprintf(
  "Cluster x AF Index: chi-sq = %.1f, Cramer's V = %.3f, p (Monte Carlo) = %s",
  chisq.result$statistic, chisq.result$cramers.v,
  format.pval( chisq.result$p.value, eps = 1 / 10000 )
) )

embed.assoc.result <- test.af.index.embedding.location(
  embed.coords = cluster.result$embed.coords,
  af.index     = align.to.embedding( cluster.result, af.index.for.cluster.raw )
)

message( sprintf(
  "Embedding location x AF Index: Pillai = %.3f, permutation p = %s (n.perm = %d)",
  embed.assoc.result$pillai.observed,
  format.pval( embed.assoc.result$p.value, eps = 1 / ( embed.assoc.result$n.perm + 1 ) ),
  embed.assoc.result$n.perm
) )

message( sprintf(
  "AF Index by cluster: %s, p = %s",
  af.violin.result$anova$method,
  format.pval( af.violin.result$anova$p.value, eps = 1e-4 )
) )

print( head( embed.assoc.result$mahalanobis.table, 10 ) )

