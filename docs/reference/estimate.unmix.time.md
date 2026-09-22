# Estimate Unmixing Time

Benchmarks a small, representative sample of events from an FCS file
through the exact unmixing dispatch
[`unmix.fcs()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.fcs.md)
would use, then scales the measured per-event rate up to the full file
to give a fast, approximate estimate of total unmixing time. Most useful
for the `AutoSpectral` method with `spectra.variants` supplied, where
per-cell variant selection makes unmixing the slowest step and its cost
is hard to predict analytically: it depends on the number of events,
detectors, fluorophores, the size of `spectra.variants`, and how many
cells in the file are above positivity thresholds and therefore trigger
variant search. Because the probe runs the real dispatch code on real
data, all of these factors are captured automatically without being
modelled explicitly.

## Usage

``` r
estimate.unmix.time(
  fcs.file,
  spectra,
  asp,
  method = c("AutoSpectral", "OLS", "WLS", "Poisson", "FastPoisson"),
  af.spectra = NULL,
  spectra.variants = NULL,
  use.dist0 = TRUE,
  speed = c("fast", "medium", "slow"),
  parallel = TRUE,
  threads = if (parallel) 0 else 1,
  n.variants = NULL,
  pipeline = c("joint", "legacy"),
  n.passes = 1L,
  n.af.passes = 1L,
  cell.weight = if (asp$cytometer == "ID7000") TRUE else FALSE,
  noise.floor = 125,
  alpha = 0.5,
  collinear.threshold = 0.5,
  joint.pair.resolution = TRUE,
  refine.af.quantile = 0.5,
  weights = NULL,
  divergence.threshold = 10000,
  divergence.handling = "Balance",
  balance.weight = 0.5,
  sample.events = 5000,
  chunk.size = 2e+06,
  verbose = TRUE
)
```

## Arguments

- fcs.file:

  A character string specifying the path to the FCS file.

- spectra:

  A matrix containing the spectral data, as passed to
  [`unmix.fcs()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.fcs.md).

- asp:

  The AutoSpectral parameter list, as passed to
  [`unmix.fcs()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.fcs.md).

- method:

  A character string specifying the unmixing method, as in
  [`unmix.fcs()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.fcs.md).
  Default `"AutoSpectral"`.

- af.spectra, spectra.variants, use.dist0, speed:

  AutoSpectral-specific arguments, passed through unchanged to the same
  unmixing dispatch
  [`unmix.fcs()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.fcs.md)
  uses; see
  [`?unmix.fcs`](https://drcytometer.github.io/AutoSpectral/reference/unmix.fcs.md)
  for details on each.

- parallel, threads, n.variants, pipeline:

  Execution and pipeline selection arguments, passed through unchanged;
  see
  [`?unmix.fcs`](https://drcytometer.github.io/AutoSpectral/reference/unmix.fcs.md).

- n.passes, n.af.passes, cell.weight, noise.floor:

  Joint-pipeline tuning arguments, passed through unchanged; see
  [`?unmix.fcs`](https://drcytometer.github.io/AutoSpectral/reference/unmix.fcs.md).

- alpha, collinear.threshold, joint.pair.resolution, refine.af.quantile:

  Further joint-pipeline tuning arguments, passed through unchanged; see
  [`?unmix.fcs`](https://drcytometer.github.io/AutoSpectral/reference/unmix.fcs.md).

- weights, divergence.threshold, divergence.handling, balance.weight:

  Weighting and IRLS-divergence arguments (`WLS`/`Poisson`/
  `FastPoisson`), passed through unchanged; see
  [`?unmix.fcs`](https://drcytometer.github.io/AutoSpectral/reference/unmix.fcs.md).

- sample.events:

  Numeric, number of events to draw for the timing probe. Default
  `5000`. Larger samples give a more stable estimate at the cost of a
  longer probe run.

- chunk.size:

  Numeric, as passed to
  [`unmix.fcs()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.fcs.md);
  used only to report how many chunks the full run will use. Default
  `2e6`.

- verbose:

  Logical, whether to print the estimate. Default `TRUE`.

## Value

Invisibly, a list with `total.events`, `sample.events`,
`events.per.second`, `estimated.unmix.seconds`,
`estimated.read.seconds`, `estimated.total.seconds`, and `chunk.n`; or
`NULL` (invisibly) if `method` is not `"AutoSpectral"` or
`spectra.variants` is `NULL`, in which case the timing probe is skipped
as unnecessary and, if `verbose`, a message explains why.

## Details

For the range of event counts typical of flow cytometry files (tens of
thousands to tens of millions), per-event unmixing cost dominates any
fixed, event-count-independent overhead (thread pool setup, small object
allocation), so total dispatch time is well approximated as directly
proportional to the number of events. A per-event rate is measured once
from a small probe of `sample.events` drawn from the middle of the file
(not the first events, to avoid acquisition start-up artefacts) and then
used to extrapolate to the full event count. This is expected to be
accurate to within roughly 4-5x under normal conditions. It will be less
accurate if the fraction of multi-positive events, which trigger more
per-cell variant search, varies substantially over the course of
acquisition, since the probe's local complexity may then not represent
the file as a whole.

File reading is timed and extrapolated the same way; file writing
([`writeFCS()`](https://drcytometer.github.io/AutoSpectral/reference/writeFCS.md))
is not estimated and is excluded from the total, since it is typically
small compared to `AutoSpectral` per-cell unmixing when
`spectra.variants` is supplied.
