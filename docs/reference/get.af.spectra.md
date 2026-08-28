# Get Autofluorescence Spectra

Extracts autofluorescence spectra from an unstained sample. Intended for
use with `unmix.autospectral`. Uses SOM clustering for rapid
identification of cells with similar AF profiles.

Optionally deduplicates the resulting spectra by cosine similarity
(`deduplicate = TRUE`, default) to remove near-identical profiles that
cause spurious over-correction of near-zero events in fully stained
samples. When `refine = TRUE`, a second round of targeted modulation is
performed on cells that remain far from zero after the first-pass
correction; modulated spectra are screened for redundancy against each
other and against the base library before being appended.

## Usage

``` r
get.af.spectra(
  unstained.sample,
  asp,
  spectra,
  som.dim = 10,
  figures = TRUE,
  save = TRUE,
  plot.dir = NULL,
  table.dir = NULL,
  title = "Autofluorescence spectra",
  verbose = TRUE,
  deduplicate = FALSE,
  duplication.threshold = 0.99,
  use.unmixed = TRUE,
  refine = TRUE,
  problem.quantile = 0.99,
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

- som.dim:

  Number of x and y dimensions for the SOM. Default is `10`.

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
  `use.unmixed = FALSE` also forces `refine = FALSE`, since the
  second-pass refinement identifies "problem cells" from per-cell
  unmixing residuals and is subject to the same instability.

- refine:

  Logical, default `FALSE`. Controls whether to perform a second round
  of autofluorescence measurement on "problem cells": those with the
  highest residual fluorophore signal after the first-pass per-cell AF
  extraction, as defined by `problem.quantile`. When `FALSE`, behavior
  is identical to versions of AutoSpectral prior to 1.0.0. If you are
  working with samples containing complex autofluorescence, e.g. tissues
  or tumors, using `refine = TRUE` will improve autofluorescence
  extraction at the cost of an increase in unmixing time.

- problem.quantile:

  Numeric, default `0.99`. The quantile for determining which cells are
  "problematic" after first-pass per-cell AF extraction. Cells at or
  above this quantile with respect to the L2 norm of their unmixed
  fluorophore channels (i.e. still furthest from zero) are selected for
  the second-round modulation. A value of `0.99` means the top 1% of
  cells.

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
rows are the deduplicated base spectra and, if `refine = TRUE`,
modulated spectra for problem cells.

## References

Van Gassen S et al. (2015). "FlowSOM: Using self-organizing maps for
87(7), 636-645.
[doi:10.1002/cyto.a.22625](https://doi.org/10.1002/cyto.a.22625) Wehrens
R, Kruisselbrink J (2018). "Flexible Self-Organizing Maps in kohonen
3.0." *Journal of Statistical Software*, *87*(7), 1-18.
[doi:10.18637/jss.v087.i07](https://doi.org/10.18637/jss.v087.i07)
