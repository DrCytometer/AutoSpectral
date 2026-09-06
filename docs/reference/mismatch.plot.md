# Plot Spectral Mismatch, Angle, Variability, Alignment, and Brightness by Dye Class

Produces a standard set of diagnostic plots comparing a particle type's
spectral mismatch, spectral angle, variability, variability/mismatch
alignment, cosine similarity, and brightness against fluorophore dye
class, plus every pairwise correlation between those six metrics with
simple linear-fit statistics. Each individual plot is saved as a JPEG in
`output.dir`. A curated subset of these, centered on spectral angle as
the primary accuracy metric, is additionally combined into a single
multi-panel PDF report.

## Usage

``` r
mismatch.plot(
  cosine.data,
  angle.data,
  mismatch.data,
  variability.data,
  alignment.data,
  brightness.data,
  fluor.df,
  particle.name = "UltraComp",
  output.dir = "./results/aurora",
  mismatch.limits = c(0, 3.5),
  sim.limits = c(1, 0.93),
  angle.limits = c(0, 25),
  alignment.limits = c(0, 1),
  cytometer
)
```

## Arguments

- cosine.data:

  One-column numeric matrix of cosine similarity values, rownames are
  fluorophore names, as returned by
  [`assess.mismatch()`](https://drcytometer.github.io/AutoSpectral/reference/assess.mismatch.md).

- angle.data:

  One-column numeric matrix of spectral angle values (degrees), rownames
  are fluorophore names, as returned by
  [`assess.mismatch.angle()`](https://drcytometer.github.io/AutoSpectral/reference/assess.mismatch.angle.md).

- mismatch.data:

  Named numeric vector of per-fluorophore mismatch magnitudes, e.g.
  `rowSums(abs(bead.cell.dist(...)))`.

- variability.data:

  Named numeric vector of per-fluorophore variability magnitudes, e.g.
  `rowSums(abs(assess.variability.mad(...)))`.

- alignment.data:

  One-column numeric matrix of per-fluorophore variability/mismatch
  alignment values, rownames are fluorophore names, as returned by
  [`assess.variability.alignment()`](https://drcytometer.github.io/AutoSpectral/reference/assess.variability.alignment.md).

- brightness.data:

  One-column numeric matrix of per-fluorophore brightness (MFI) values,
  rownames are fluorophore names, as returned by
  [`get.brightness.automated()`](https://drcytometer.github.io/AutoSpectral/reference/get.brightness.automated.md).

- fluor.df:

  Data frame with at least `Fluorophore` and `Class` columns, used to
  annotate each fluorophore with a dye class for the violin plots.

- particle.name:

  Character. Particle type label used in plot titles and output
  filenames. Default `"UltraComp"`.

- output.dir:

  Character. Directory for the saved JPEGs and PDF report. Default
  `"./results/aurora"`.

- mismatch.limits:

  Numeric vector of length 2, y-axis limits for mismatch plots. Default
  `c(0, 3.5)`.

- sim.limits:

  Numeric vector of length 2, y/x-axis limits for cosine similarity
  plots (reversed axis). Default `c(1, 0.93)`.

- angle.limits:

  Numeric vector of length 2, y/x-axis limits for spectral angle plots.
  Default `c(0, 25)`.

- alignment.limits:

  Numeric vector of length 2, y/x-axis limits for variability/mismatch
  alignment plots. Default `c(0, 1)`.

- cytometer:

  Character. Cytometer label used in plot titles and output filenames.

## Value

A data frame of pairwise linear-fit statistics (R-squared and p-value)
for each pair of metrics, one row per comparison. Returns
`invisible(NULL)` if fewer than 5 fluorophores have complete data.

## Details

Six metrics are compared: `Mismatch`, `Cosine`, `Angle`, `Variability`,
`Alignment`, and `Brightness`. For each metric, a violin/jitter plot by
dye class is produced. For every pair of metrics, a scatter plot with an
[`lm()`](https://rdrr.io/r/stats/lm.html) trendline is produced,
annotated with that pair's R-squared and p-value; since simple linear
regression's R-squared and F-test p-value are symmetric in the two
variables, these statistics are unaffected by which metric in a pair
ends up on the plot's x- versus y-axis. `Cosine` is the only metric
plotted on a reversed axis (cosine similarity decreases with divergence,
unlike the other five metrics, which all increase with divergence).

The JPEGs cover all six metrics and every pairwise combination. The
consolidated PDF report is narrower and fixes `Angle` as the reference
metric throughout: four violin plots (`Angle`, `Variability`,
`Alignment`, `Brightness`) followed by three pairwise scatter plots
(`Variability` vs `Angle`, `Alignment` vs `Angle`, `Brightness` vs
`Angle`), each with spectral angle fixed on the y-axis. `Mismatch` and
`Cosine` are excluded from the PDF report; they remain available as
individual JPEGs.
