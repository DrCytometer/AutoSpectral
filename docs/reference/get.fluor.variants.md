# Get Fluorophore Variants

Assesses variation in the spectral signature of a single-stained flow
cytometry control sample using SOM clustering on scatter-matched,
background-corrected positive events.

Autofluorescence is characterised **in situ** from the paired
universal-negative file specified in the control table (or from the
lower 25\\ available). The AF mean vector is unit-normalised and
projected out of each event to identify the empirical peak detector,
mirroring the approach in `get.spectra.automated`. All positive events
above the raw threshold in the (empirical) peak channel are selected, up
to `n.cells` (randomly downsampled when more are present). For each
selected event, \\k\\ scatter-space nearest neighbours are found in the
unstained pool and their spectral values are averaged to form a
per-event background estimate, which is then subtracted. SOM clustering
on the resulting background-corrected matrix recovers the
population-level distribution of spectral shapes. A cosine-similarity QC
step retains only SOM centroids sufficiently similar to the reference
spectrum, followed by off-peak smoothing.

The Spillover Spreading Matrix inputs are now computed differently:
instead of re-unmixing this control's positive population against the
full panel and taking an empirical MAD, this fits a well-conditioned,
low-rank model (autofluorescence components plus this fluorophore's own
reference row – or this fluorophore alone for beads, or a cell control
too collinear with AF to separate safely) and derives spread from the
residual of that fit.

The low-rank fit is exactly the AF-projection step
`get.fluor.variants()` already performs; the difference is that its
residual and this control's own recovered abundance are kept rather than
discarded. Because the residual is, by construction, orthogonal to every
column of this small design, its magnitude is not amplified by
collinearity elsewhere in the panel the way a full `n`-parameter OLS
solve of the raw positive population is. Projecting that residual
through the full-panel pseudoinverse
([`unmix.ols()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.ols.md))
and robustly regressing its squared value against this control's own
recovered abundance – via `.fix.huber.slope()` (falling back to
[`fit.robust.linear.model()`](https://drcytometer.github.io/AutoSpectral/reference/fit.robust.linear.model.md)
if `AutoSpectralRcpp` is unavailable) – gives, for every target
fluorophore at once, the same intercept-plus-slope relationship the
Residual Model of Cai et al. (bioRxiv 2026, USERM) derives analytically:
`Var(unmixed_k) = intercept_k + slope_k * b_h`, where `b_h` is this
control's own on-channel abundance. `slope_k` is exactly `SS(fluor, k)`
in the units used throughout this package (added variance per unit of
on-channel abundance), with no separate unstained-population baseline
required – the regression's own intercept already is that baseline,
estimated across this control's own dim-to- bright range rather than
from a different sample and a different estimator entirely.

This uses every event that clears the raw peak-channel threshold
(`pos.idx`), not just the narrower `keep.idx`/cosine-QC'd subset used
for SOM clustering, so the regression sees the widest available range of
on-channel brightness and does not depend on cosine QC succeeding – a
control that fails cosine QC (and would otherwise fall back to a cruder
spread estimate) gets the identical Residual Model treatment here.

## Usage

``` r
get.fluor.variants(
  fluor,
  file.name,
  control.dir,
  asp,
  spectra,
  figures,
  output.dir,
  verbose,
  spectral.channel,
  scatter.channel,
  universal.negative,
  control.type,
  raw.thresholds,
  unmixed.thresholds,
  flow.channel,
  af.pcs,
  neg.cache = NULL,
  use.unmixed = TRUE,
  n.cells = 10000L,
  som.dim = 10L,
  k.neighbors = 3L,
  sim.threshold = 0.985,
  sim.threshold.floor = 0.9,
  af.collinear.threshold = 0.95,
  noise.floor.tail.fraction = 0.2,
  noise.mask.threshold = 0.05,
  noise.n.cells = 2000L,
  huber.k = 1.345,
  huber.max.iter = 100L,
  variant.fill.color = "red",
  variant.fill.alpha = 0.7,
  median.line.color = "black",
  median.linewidth = 1,
  parallel = TRUE,
  threads = NULL
)
```

## Arguments

- fluor:

  Character. Name of the fluorophore.

- file.name:

  Named character vector of control FCS filenames, named by fluorophore.

- control.dir:

  Character. Directory containing the control FCS files.

- asp:

  The AutoSpectral parameter list from
  [`get.autospectral.param()`](https://drcytometer.github.io/AutoSpectral/reference/get.autospectral.param.md).

- spectra:

  Numeric matrix. Reference spectra; fluorophores in rows, detectors in
  columns.

- figures:

  Logical. Whether to save a spectral-variant plot. Default `TRUE`.

- output.dir:

  Character. Directory for figures.

- verbose:

  Logical. Whether to print progress messages. Default `TRUE`.

- spectral.channel:

  Character vector of spectral detector channel names.

- scatter.channel:

  Character vector of scatter parameter names (e.g. `"FSC-A"`,
  `"SSC-A"`) used for KNN scatter matching against the unstained pool.

- universal.negative:

  Named character vector mapping fluorophore names to their paired
  unstained FCS filename, or `"FALSE"` / `NA` when none is available.

- control.type:

  Character, either "beads" or "cells". Determines the type of control
  sample being used and the subsequent processing steps.

- raw.thresholds:

  Named numeric vector of per-channel positivity thresholds (typically
  the 99.5th percentile of the unstained sample).

- unmixed.thresholds:

  A named vector of numerical values corresponding to the threshold for
  positivity in each unmixed channel. Determined by the 99.5th
  percentile on the unstained sample, typically after single-cell AF
  unmixing.

- flow.channel:

  Named character vector of expected peak raw channels, one per
  fluorophore.

- af.pcs:

  Named list of autofluorescence-defining principal component matrices,
  one per unique unstained FCS file. Names are FCS filenames matching
  entries in `universal.negative`.

- neg.cache:

  Optional named list, one entry per unique unstained FCS file (names
  are filenames), each holding `$spectral` and `$scatter` matrices
  already read from disk. Built once in
  [`get.spectral.variants()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectral.variants.md)
  and reused across all fluorophores sharing the same negative, avoiding
  a repeat disk read/decode per fluorophore. Default `NULL` falls back
  to reading `universal.negative[fluor]` directly.

- use.unmixed:

  Logical, default `TRUE`. Whether to unmix background-corrected
  positive events against the full `spectra` matrix for positivity
  selection and SOM clustering features, and whether the Residual
  Model's residual-projection step runs at all. Set to `FALSE` when
  `spectra` contains several similar or collinear fluorophores (e.g. a
  bead-cell comparison panel), where the full-spectra unmix is itself
  unstable or unsolvable. When `FALSE`, clustering falls back to raw
  detector space only, and the `"spillover.spread"` family of attributes
  is not computed (see Value).

- n.cells:

  Integer, default `10000`. Maximum number of positive events used for
  SOM clustering. Files with more events above threshold are randomly
  downsampled to this number. Does not limit the event set used for the
  Residual Model regression, which uses every event in `pos.idx`.

- som.dim:

  Integer, default `10`. Side length of the square SOM grid. Produces up
  to `som.dim^2` candidate variant spectra before cosine QC.

- k.neighbors:

  Integer, default `3`. Number of scatter-space nearest neighbours from
  the unstained pool used to form the per-event background estimate.

- sim.threshold:

  Numeric, default `0.985`. Minimum cosine similarity between a SOM
  centroid and the reference spectrum for the centroid to be retained as
  a variant.

- sim.threshold.floor:

  Numeric, default `0.90`. Lower bound for adaptive relaxation of
  `sim.threshold` when the initial cutoff retains fewer than 20 events.
  Relaxation is logged via
  [`warning()`](https://rdrr.io/r/base/warning.html) and the threshold
  actually used is returned as the `"cosine.threshold.used"` attribute.

- af.collinear.threshold:

  Numeric, default `0.95`. Minimum cosine similarity between `fluor`'s
  reference spectrum and any of its paired unstained file's AF principal
  directions (`af.pcs`) at or above which the AF-component projection
  step – and the low-rank fit feeding the Residual Model – drops the
  AF-PC terms and uses this fluorophore's reference row alone, since a
  joint OLS fit against near-collinear AF and fluorophore directions can
  push real fluorophore signal into the AF term. Recorded as the
  `"af.collinear"` attribute.

- noise.floor.tail.fraction:

  Numeric in (0, 1), default `0.20`. Per-detector noise floor is the MAD
  (scaled to a Gaussian-equivalent variance) of the lowest fraction of
  raw values in that detector's column, using every event in the control
  file. Lower values isolate a purer background tail but with fewer
  events to estimate from; higher values are more stable but risk
  pulling in dim positive events.

- noise.mask.threshold:

  Numeric in (0, 1), default `0.05`. A detector is excluded from this
  control's contribution to the pooled noise-model regression (see
  `"noise.mask"` below) when this fluorophore's own reference spectrum
  exceeds this fraction of its own peak there – the channels where
  unmodelled spectral-variant wobble, not photon noise, would otherwise
  dominate the residual.

- noise.n.cells:

  Integer, default `2000L`. Maximum events used for the noise-model
  residual (see `"noise.resid"` below), sampled from `keep.idx` – the
  same unambiguously-positive events already selected for SOM input. The
  negative/dim majority of a single-stained control carries no
  information the unstained sample doesn't already supply, and pools to
  a very large, AF-dominated mass across many controls; restricting to
  the positive gate avoids re-fitting that problem at high cost.

- huber.k:

  Numeric, default `1.345`. Huber tuning constant passed to
  `.fix.huber.slope()` for the spillover-spread regression.

- huber.max.iter:

  Integer, default `100L`. Maximum IRLS iterations passed to
  `.fix.huber.slope()`.

- variant.fill.color:

  Color for the shaded ribbon in the variant plot. Default `"red"`.

- variant.fill.alpha:

  Alpha for `variant.fill.color`. Default `0.7`.

- median.line.color:

  Color for the reference-spectrum line. Default `"black"`.

- median.linewidth:

  Width of the reference-spectrum line. Default `1`.

- parallel:

  Logical, default `TRUE`. Enable OpenMP multi-threading for this call's
  batch SOM
  ([`get.som.codes()`](https://drcytometer.github.io/AutoSpectral/reference/get.som.codes.md)).

- threads:

  Numeric or `NULL`. OpenMP threads for the batch SOM. `NULL` defaults
  to `0` (all available cores) when `parallel = TRUE`.

## Value

A numeric matrix; variants in rows, detectors in columns, values
normalised to \\\[0, 1\]\\. Row 1 is always the library reference
spectrum for `fluor`; subsequent rows are SOM-derived variants that
passed cosine QC. When too few positive events are available, or no
centroids survive cosine QC, the single reference spectrum is returned
(one row), still carrying the full set of attributes below. Carries:
`"noise.floor"` (per-detector background SD, described above),
`"noise.events"` (up to `noise.n.cells` events x detectors matrix,
sampled from `keep.idx`, background-corrected – and AF-corrected where
the projection above ran – but not otherwise fit; pooled and fit jointly
against the full panel in
[`get.spectral.variants()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectral.variants.md)),
`"noise.mask"` (logical vector, length `ncol(spectra)`, `TRUE` at
detectors this fluorophore's own spectrum dominates),
`"spillover.spread"` (named numeric vector, one Huber-robust slope per
target fluorophore – `SS(fluor, .)` – or `NULL` if fewer than 20 events
were available or `use.unmixed = FALSE`), `"spillover.spread.intercept"`
(named numeric vector, the matching intercept – this control's own
estimate of each channel's baseline unmixed variance at zero abundance),
`"on.channel.mfi"` (this control's own median recovered abundance from
the low-rank fit, `NA` under the same conditions as
`"spillover.spread"`), `"spillover.spread.n"` (integer, the number of
events behind the regression), `"spillover.spread.range"` (numeric, the
range of this control's own recovered abundance used in the regression –
a control whose events barely vary in brightness gives the slope little
to fit against, regardless of event count), and
`"spillover.spread.source"` (character, `"residual"` when computed, `NA`
otherwise).

## References

Van Gassen S et al. (2015). FlowSOM. *Cytometry Part A*, 87(7), 636-645.
[doi:10.1002/cyto.a.22625](https://doi.org/10.1002/cyto.a.22625)

Cai X et al. (2026). Residual Model for unmixed spread prediction in
spectral flow cytometry. *bioRxiv* 2026.01.27.701929.
