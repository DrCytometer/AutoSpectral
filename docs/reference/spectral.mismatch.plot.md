# Spectral Mismatch Plot

Produces a two-panel diagnostic figure comparing two fluorophore
spectra. The upper panel overlays the normalised spectral traces for a
reference and a test spectrum across all detectors. The lower panel
shows the per-detector spectral difference (reference minus test) as a
bar chart, making it easy to identify channels where the spectra
diverge.

The combined figure is saved as a JPEG to `plot.dir` and returned
invisibly as a `ggplot` / cowplot object so it can be embedded in
downstream reports or further modified.

## Usage

``` r
spectral.mismatch.plot(
  ref.spectrum,
  test.spectrum,
  ref.label = NULL,
  test.label = NULL,
  title = "Spectral comparison",
  bar.color = "lightblue",
  fluorophore.colors = c("black", "red"),
  plot.dir = NULL,
  line.width = 1,
  point.size = 1,
  color.palette = NULL,
  show.legend = TRUE,
  plot.width = NULL,
  plot.height = NULL
)
```

## Arguments

- ref.spectrum:

  A named numeric vector (or single-row matrix) of normalised spectral
  intensities for the reference spectrum. Names must correspond to
  detector channel names.

- test.spectrum:

  A named numeric vector (or single-row matrix) of normalised spectral
  intensities for the test spectrum. Must have the same names as
  `ref.spectrum`.

- ref.label:

  Character string used as the legend label for the reference spectrum.
  If `NULL` (default), the
  [`rownames()`](https://rdrr.io/r/base/colnames.html) of `ref.spectrum`
  is used.

- test.label:

  Character string used as the legend label for the test spectrum. If
  `NULL` (default), the
  [`rownames()`](https://rdrr.io/r/base/colnames.html) of
  `test.spectrum` is used.

- title:

  Character string used as both the plot title and the JPEG filename
  stem. Default is `"Spectral comparison"`.

- bar.color:

  Fill color for the spectral-difference bar chart. Default is
  `"lightblue"`.

- fluorophore.colors:

  Character vector of length 2 giving the line colors for the reference
  and test spectra respectively. Overrides `color.palette` when
  non-`NULL`. Default is `c("black", "red")`.

- plot.dir:

  Directory to save the JPEG file. If `NULL` (default), the current
  working directory is used. Created automatically if absent.

- line.width:

  Numeric. Line width for the spectral traces. Default `1`.

- point.size:

  Numeric. Point size for the spectral traces. Default `1`.

- color.palette:

  Optional viridis palette name (e.g., `"viridis"`, `"plasma"`) used
  when `fluorophore.colors` is `NULL`. Ignored when `fluorophore.colors`
  is set. Default is `NULL`.

- show.legend:

  Logical. Whether to display the color legend on the upper panel.
  Default is `TRUE`.

- plot.width:

  Numeric. Width of the saved figure in inches. If `NULL` (default), an
  automatic width is derived from the number of detectors (approximately
  12 inches per 64 channels, minimum 3 inches).

- plot.height:

  Numeric. Height of the saved figure in inches. If `NULL` (default), an
  automatic height is derived from the number of spectra rows.

## Value

The combined cowplot object is returned invisibly. The figure is always
saved to disk as `<plot.dir>/<title>.jpg`.
