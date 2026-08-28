# Plot Spectral Mismatch, Variability, and Brightness by Dye Class

Produces a standard set of diagnostic plots comparing a particle type's
spectral mismatch, variability, cosine similarity, and brightness
against fluorophore dye class, plus pairwise correlations between all
four metrics with simple linear-fit statistics. Each individual plot is
saved as a JPEG in `output.dir`, and all plots are additionally combined
into a single multi-panel PDF report.

## Usage

``` r
mismatch.plot(
  cosine.data,
  mismatch.data,
  variability.data,
  brightness.data,
  fluor.df,
  particle.name = "UltraComp",
  output.dir = "./results/aurora",
  mismatch.limits = c(0, 3.5),
  sim.limits = c(1, 0.93),
  cytometer
)
```

## Arguments

- cosine.data:

  One-column numeric matrix of cosine similarity values, rownames are
  fluorophore names, as returned by
  [`assess.mismatch()`](https://drcytometer.github.io/AutoSpectral/reference/assess.mismatch.md).

- mismatch.data:

  Named numeric vector of per-fluorophore mismatch magnitudes, e.g.
  `rowSums(abs(bead.cell.dist(...)))`.

- variability.data:

  Named numeric vector of per-fluorophore variability magnitudes, e.g.
  `rowSums(abs(assess.variability.mad(...)))`.

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

- cytometer:

  Character. Cytometer label used in plot titles and output filenames.

## Value

A data frame of pairwise linear-fit statistics (R-squared and p-value)
for each pair of metrics, one row per comparison. Returns
`invisible(NULL)` if fewer than 5 fluorophores have complete data.
