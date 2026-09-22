# Get Autofluorescence Spectra

Extracts autofluorescence spectra from an unstained sample. Intended for
use with `unmix.autospectral`. Uses SOM clustering for rapid
identification of cells with similar AF profiles.

Optionally deduplicates the resulting spectra by cosine similarity
(`deduplicate = TRUE`, default) to remove near-identical profiles that
cause spurious over-correction of near-zero events in fully stained
samples. When `refine = TRUE`, a second round of discovery targets cells
that remain far from zero after the first-pass correction: candidate
spectra are built from density-boosted neighbourhoods of these problem
cells and accepted only when the real per-cell AF solver demonstrably
prefers them.

## Usage

``` r
get.af.spectra(
  unstained.sample,
  asp,
  spectra,
  unstained.exprs = NULL,
  som.dim = 10,
  dist = 2L,
  figures = TRUE,
  save = TRUE,
  plot.dir = NULL,
  table.dir = NULL,
  title = "Autofluorescence spectra",
  verbose = TRUE,
  af.assign.method = c("l1", "l2"),
  deduplicate = FALSE,
  duplication.threshold = 0.99,
  use.unmixed = TRUE,
  af.basis.components = NULL,
  raw.pca.components = NULL,
  refine = FALSE,
  k.neighbors = 15L,
  refine.improvement.threshold = 0.02,
  refine.min.shift.n = 8L,
  problem.quantile = 0.99,
  plot.unmixed = FALSE,
  plot.unmixed.n = 30000,
  remove.contaminants = TRUE,
  contaminant.threshold = 0.99,
  parallel = TRUE,
  threads = if (parallel) 0 else 1,
  return.model = FALSE,
  model.rank = 6L,
  model.var.explained = 0.95,
  model.min.events = 50L,
  model.shrinkage = 0.1,
  heatmap.color.palette = "viridis",
  spectral.trace.color.palette = NULL,
  af.fill.color = "red",
  af.line.color = "black"
)
```

## Arguments

- unstained.sample:

  Path and file name for an unstained sample FCS file. The sample type
  and processing (protocol) method should match the fully stained
  samples to which the AF will be applied, ideally.

- asp:

  The AutoSpectral parameter list. Prepare using
  `get.autospectral.param`.

- spectra:

  Spectral signatures of fluorophores, normalized between 0 and 1, with
  fluorophores in rows and detectors in columns.

- unstained.exprs:

  Optional matrix of data from the unstained sample. The columns must
  match `spectra`. Default is `NULL`.

- som.dim:

  Number of x and y dimensions for the SOM. Default is `10`.

- dist:

  Integer 1:4, distance function (1 manhattan, 2 euclidean, 3 chebyshev,
  4 cosine). Default `2`.

- figures:

  Logical, whether to plot the spectral traces and heatmap for the AF
  signatures. Default is `TRUE`.

- save:

  Logical, whether to save the CSV file for the AF signatures. Default
  is `TRUE`.

- plot.dir:

  Directory (folder) where the plots will be saved. Default is `NULL`,
  which inherits from `asp$figure.af.dir`.

- table.dir:

  Directory (folder) where the spectra csv file will be saved. Default
  is `NULL`, which inherits from `asp$table.af.dir`.

- title:

  Title for the output spectral plots and csv file. Default is
  `"Autofluorescence spectra"`.

- verbose:

  Logical, controls messaging. Default is `TRUE`.

- af.assign.method:

  Character, one of `"l1"` (default) or `"l2"`. Controls which per-cell
  AF assignment solver is used everywhere inside this function: the
  first-pass assignment, the refinement loop's candidate-preference
  check, and the second-pass diagnostic unmixing
  (`plot.unmixed = TRUE`). Both options score the joint covariance-
  weighted fluorophore error x raw-space residual error criterion,
  differing only in whether the fluorophore term is L1 (abs, `"l1"`,
  `assign.af.joint.cov`) or L2 (squared, `"l2"`,
  `assign.af.joint.cov.l2`). Each prefers its compiled Rcpp fast path
  ([`AutoSpectralRcpp::assign.af.joint.cov.fast`](https://rdrr.io/pkg/AutoSpectralRcpp/man/assign.af.joint.cov.fast.html)
  / `assign.af.joint.cov.l2.fast`) when available and falls back to pure
  R otherwise.

- deduplicate:

  Logical, default `FALSE`. Whether to deduplicate AF spectra by cosine
  similarity after the base clustering stage and again after the
  refinement stage. Deduplication removes near-identical spectral
  profiles that can cause overzealous matching of near-zero events in
  fully stained samples, reducing apparent "squishing". Deduplication is
  slightly less accurate. Set to `TRUE` to us it.

- duplication.threshold:

  Numeric, default `0.99`. The cosine similarity threshold used for
  deduplication. A spectrum is dropped if its cosine similarity to any
  already-retained spectrum meets or exceeds this value. Only used when
  `deduplicate = TRUE`.

- use.unmixed:

  Logical, default `TRUE`. Whether to include an OLS unmixed-space
  projection (`unmix.ols.fast(unstained.exprs, spectra)`) alongside the
  raw detector data as SOM clustering input. Set to `FALSE` to cluster
  on raw detector space only, which is appropriate when `spectra`
  contains several similar or collinear fluorophores (e.g. a bead-cell
  comparison panel), where an OLS unmix is itself unstable and would
  corrupt the clustering features rather than enrich them.
  `use.unmixed = FALSE` also forces `refine = FALSE` and
  `plot.unmixed = FALSE`, since both rely on the same per-cell OLS
  unmixing residuals and are subject to the same instability.

- af.basis.components:

  Integer, default `NULL`. When supplied, appends this many components
  of the panel-oblique autofluorescence structure of the unstained
  sample to the SOM training features, alongside the raw detector data
  and (when `use.unmixed = TRUE`) the unmixed coefficients. Computed
  directly from the panel residual of the unstained events (trimming the
  top 1% by residual norm, subsampling to 50000 events for the
  decomposition, then projecting every event), reusing the unmixing
  matrix already computed for `unmixed.no.af` rather than re-deriving
  it. Training-time use carries none of the collinearity risk a per-cell
  regression design would, since clustering is not a regression; this
  only gives the SOM extra shape-discriminating features to place nodes
  with. Forced to `NULL` when `use.unmixed = FALSE`, since it depends on
  the same unmixing matrix. Default `NULL` disables it, matching prior
  behaviour.

- raw.pca.components:

  Integer, default `NULL`. When supplied, appends this many raw-event
  principal components of the unstained sample (ordinary singular
  vectors of the raw detector data, not restricted to the part the panel
  cannot explain, subsampled to 50000 events for the decomposition) to
  the SOM training features. Unlike `af.basis.components`, these
  components are free to cross into the panel's span; that carries no
  collinearity risk here for the same reason as above. Not affected by
  `use.unmixed`, since it does not use the unmixing matrix at all.
  Default `NULL` disables it.

- refine:

  Logical, default `FALSE`. Controls whether to perform a second round
  of autofluorescence discovery targeting "problem cells": those with
  the highest residual fluorophore signal after the first-pass per-cell
  AF extraction, as defined by `problem.quantile`. Problem cells are
  grouped by the systematic pattern of their error, and each group's
  cells act as seeds for a nearest-neighbour expansion across the full
  unstained population (`k.neighbors`), so a candidate spectrum is built
  from a locally density-boosted set rather than the (typically sparse)
  seeds alone. A candidate is appended only if, once added to the
  library, the real per-cell AF solver reassigns a sufficient number of
  its own seed cells onto it (`refine.min.shift.n`) and those cells show
  a genuine, paired improvement in fit rather than one that only looks
  better because a candidate was added (`refine.improvement.threshold`).
  When `FALSE`, behavior is identical to versions of AutoSpectral prior
  to 1.0.0. If you are working with samples containing complex
  autofluorescence, e.g. tissues or tumors, using `refine = TRUE` will
  improve autofluorescence extraction at the cost of an increase in
  unmixing time.

- k.neighbors:

  Integer, default `15L`. Used only when `refine = TRUE`. Each error
  cluster's problem cells act as seeds into a nearest-neighbour search
  across the full unstained population, so a candidate spectrum is built
  from a locally density-boosted set rather than from the (typically too
  sparse) problem cells alone. Higher values recruit a larger, more
  stable candidate at the cost of reaching further from the seeds and
  risking dilution by unrelated bulk events.

- refine.improvement.threshold:

  Numeric, default `0.02`. Minimum median gain in cosine similarity (raw
  event to assigned AF spectrum), among the seed cells that shift their
  assignment onto a candidate under the real per-cell solver, before
  that candidate is accepted. The comparison is paired per cell, so
  adding more candidates cannot inflate it. Also requires the 25th
  percentile of that per-cell gain to be positive.

- refine.min.shift.n:

  Integer, default `8L`. Minimum number of a cluster's seed cells that
  must shift their AF assignment onto a candidate spectrum before that
  candidate is evaluated at all. Below this, a median gain is too
  volatile to trust.

- problem.quantile:

  Numeric, default `0.99`. The quantile for determining which cells are
  "problematic" after first-pass per-cell AF extraction. Cells at or
  above this quantile with respect to the L2 norm of their unmixed
  fluorophore channels (i.e. still furthest from zero) are selected for
  the second-round refinement.

- plot.unmixed:

  Logical, default `FALSE`. Whether to unmix the unstained sample before
  and after AF extraction and plot the comparison as a single
  side-by-side biplot (`unmixed.no.af`, `unmixed`, and, when
  `refine = TRUE` and modulation succeeds, `unmixed.second`). This runs
  a full per-cell AF unmixing pass purely for diagnostic plotting, so it
  defaults off. When `refine = FALSE`, the comparison plot is drawn
  immediately from the first-pass unmixing. When `refine = TRUE`,
  plotting is deferred until the refinement loop finishes, so the plot
  always reflects the final (possibly modulated) AF spectra rather than
  an intermediate state.

- plot.unmixed.n:

  Integer, default `30000`. When `plot.unmixed = TRUE` and
  `refine = FALSE`, the unstained sample is subsampled to this many
  events before the diagnostic unmixing pass, since the full refinement
  population isn't otherwise needed. Ignored when `refine = TRUE`, since
  the full population is already required for problem-cell
  identification.

- remove.contaminants:

  Logical, default `TRUE`. A QC check is performed to exclude any
  autofluorescence spectrum that is nearly identical to a fluorophore
  signature in `spectra`. This guards against low-level contamination of
  the unstained sample by single-stained controls.

- contaminant.threshold:

  Numeric, default `0.99`. When `remove.contaminants = TRUE`, events in
  the unstained sample whose cosine similarity to any fluorophore
  spectrum in `spectra` meets or exceeds this value are removed
  **before** SOM construction. This per-event filter is more sensitive
  than the post-SOM centroid check because contaminating events are
  unlikely to dominate an entire SOM node. Lower values are more
  aggressive; the practical range is roughly 0.98–0.999.

- parallel:

  Logical, default `TRUE`, which enables parallel processing for
  per-cell AF identification. Used when `refine = TRUE`.

- threads:

  Numeric, defaults to a single thread for sequential processing
  (`parallel = FALSE`) or all available cores if `parallel = TRUE`. Used
  when `refine = TRUE`.

- return.model:

  Logical. When `TRUE`, attaches an `"af.model"` attribute to the
  returned spectra containing per-node covariance, occupancy priors,
  abundance priors and scatter statistics, for use by
  [`unmix.af.gls()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.af.gls.md).
  The return value is still a matrix, so existing callers are
  unaffected. Default `FALSE`. Requires refine = FALSE for
  well-populated per-node covariances.

- model.rank:

  Integer, maximum rank retained for each node's spectral covariance.
  Default `6`.

- model.var.explained:

  Numeric, fraction of within-node variance to retain. Default `0.95`.

- model.min.events:

  Integer, minimum events for a node to receive a covariance estimate.
  Nodes below this get the pooled covariance. Default `50`.

- model.shrinkage:

  Numeric in `[0, 1]`, shrinkage of each node covariance toward the
  pooled covariance. Guards nodes with few events. Default `0.10`.

- heatmap.color.palette:

  Optional character string defining the viridis color palette for the
  fluorophore heatmap. Default is `"viridis"`. Options: `"magma"`,
  `"inferno"`, `"plasma"`, `"viridis"`, `"cividis"`, `"rocket"`,
  `"mako"`, `"turbo"`.

- spectral.trace.color.palette:

  Optional character string defining the color palette for the AF
  traces. Default is `NULL` (default R Brewer colors). Options: same as
  `heatmap.color.palette`.

- af.fill.color:

  Color for the shaded region indicating the range of autofluorescence
  variation in the variant plot. Default is `"red"`.

- af.line.color:

  Color for the median autofluorescence line in the variant plot.
  Default is `"black"`.

## Value

A matrix of autofluorescence spectra (spectra in rows, detectors in
columns). Row 1 is the population mean of the base spectra; subsequent
rows are the deduplicated base spectra and, if `refine = TRUE`, any
solver- validated spectra discovered from under-represented problem
cells.

## References

Van Gassen S et al. (2015). "FlowSOM: Using self-organizing maps for
87(7), 636-645.
[doi:10.1002/cyto.a.22625](https://doi.org/10.1002/cyto.a.22625) Wehrens
R, Kruisselbrink J (2018). "Flexible Self-Organizing Maps in kohonen
3.0." *Journal of Statistical Software*, *87*(7), 1-18.
[doi:10.18637/jss.v087.i07](https://doi.org/10.18637/jss.v087.i07)
