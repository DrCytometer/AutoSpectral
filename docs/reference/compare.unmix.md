# Compare Unmixing Quality Across Two Spectral References

A diagnostic utility for comparing unmixing performance when two
different spectral references are available for the same fluorophore —
for example, a cell-derived reference versus a bead-derived reference,
or two different lot measurements. Both references are used to unmix the
same single-stained FCS file and the resulting Secondary Stain Index
(SSI) and secondary spillover values are compared side-by-side.

The function:

1.  Reads the single-stained and unstained FCS files.

2.  Applies a FSC/SSC scatter gate (auto-detected or user-supplied).

3.  Unmixes both gated datasets with each spectral reference matrix
    using OLS (or WLS for the ID7000).

4.  Identifies the fluorophore channel with the largest absolute SSI
    across both unmixings and uses it as the y-axis for biplot
    visualisation.

5.  Annotates each biplot with the median ± rSD of the unstained
    distribution and the SSI / spillover value.

6.  Saves a side-by-side JPEG to `plot.dir` and returns a summary table.

**Note:** This function calls
[`calculate.ssi()`](https://drcytometer.github.io/AutoSpectral/reference/calculate.ssi.md),
which must be available in the package namespace.

## Usage

``` r
compare.unmix(
  single.stained.fcs,
  unstained.fcs,
  fluorophore,
  spectra,
  ref.spectrum,
  test.spectrum,
  cytometer,
  ref.spectra = NULL,
  test.spectra = NULL,
  x.channel = NULL,
  y.channel = NULL,
  x.min = -1000,
  y.min = -1000,
  x.width.basis = -1000,
  y.width.basis = -1000,
  gate = TRUE,
  gate.bound = NULL,
  ref.label = "Cells",
  test.label = "Beads",
  biplot.width = 5,
  biplot.height = 5,
  title = paste(ref.label, "vs", test.label, fluorophore),
  plot.dir = "./figure_compare_unmix"
)
```

## Arguments

- single.stained.fcs:

  Character string. Path to the single-stained FCS file for the
  fluorophore of interest.

- unstained.fcs:

  Character string. Path to the matched unstained (negative control) FCS
  file.

- fluorophore:

  Character string. Name of the target fluorophore (must match a row
  name in `spectra`, or will be appended automatically).

- spectra:

  Numeric matrix of spectral references (fluorophores × detectors,
  values normalised 0–1) used as the starting point for constructing
  `ref.spectra` and `test.spectra` when those are not supplied directly.

- ref.spectrum:

  Named numeric vector or single-row matrix. The reference spectrum for
  `fluorophore` used in the first unmixing.

- test.spectrum:

  Named numeric vector or single-row matrix. The test spectrum for
  `fluorophore` used in the second unmixing.

- cytometer:

  Character string identifying the cytometer model, passed to
  [`get.autospectral.param()`](https://drcytometer.github.io/AutoSpectral/reference/get.autospectral.param.md).
  Also determines the unmixing algorithm: `"id7000"` uses WLS, all
  others use OLS.

- ref.spectra:

  Optional numeric matrix. A pre-built full spectral reference matrix
  (fluorophores × detectors) for the reference condition. When `NULL`
  (default), this is constructed from `spectra` and `ref.spectrum`.

- test.spectra:

  Optional numeric matrix. A pre-built full spectral reference matrix
  for the test condition. When `NULL` (default), this is constructed
  from `spectra` and `test.spectrum`.

- x.channel:

  Unused. Reserved for future use.

- y.channel:

  Unused. Reserved for future use.

- x.min:

  Numeric. Minimum x-axis value (in data units) passed to
  [`create.biplot()`](https://drcytometer.github.io/AutoSpectral/reference/create.biplot.md).
  Default `-1000`.

- y.min:

  Numeric. Minimum y-axis value (in data units) passed to
  [`create.biplot()`](https://drcytometer.github.io/AutoSpectral/reference/create.biplot.md).
  Default `-1000`.

- x.width.basis:

  Numeric. Width basis for the biexponential x-axis transform passed to
  [`create.biplot()`](https://drcytometer.github.io/AutoSpectral/reference/create.biplot.md).
  Default `-1000`.

- y.width.basis:

  Numeric. Width basis for the biexponential y-axis transform passed to
  [`create.biplot()`](https://drcytometer.github.io/AutoSpectral/reference/create.biplot.md).
  Default `-1000`.

- gate:

  Logical. Reserved for future gating control; currently the scatter
  gate is always applied. Default `TRUE`.

- gate.bound:

  Optional list with elements `$x` and `$y` defining a polygon gate
  boundary in FSC/SSC space. When `NULL` (default), the gate is detected
  automatically via
  [`do.gate()`](https://drcytometer.github.io/AutoSpectral/reference/do.gate.md).

- ref.label:

  Character string. Display label for the reference condition (e.g.,
  `"Cells"`). Also passed as `control.type` to
  [`do.gate()`](https://drcytometer.github.io/AutoSpectral/reference/do.gate.md).
  Default `"Cells"`.

- test.label:

  Character string. Display label for the test condition (e.g.,
  `"Beads"`). Default `"Beads"`.

- biplot.width:

  Numeric. Width of each individual biplot in inches. The saved figure
  is twice this wide (two panels). Default `5`.

- biplot.height:

  Numeric. Height of the saved figure in inches. Default `5`.

- title:

  Character string. Title displayed on the figure and used as the JPEG
  filename stem. Default is constructed from `ref.label`, `test.label`,
  and `fluorophore`.

- plot.dir:

  Character string. Directory for saving the output JPEG. Created
  automatically if absent. Default `"./figure_compare_unmix"`.

## Value

A data frame with one row per fluorophore in the reference spectra
matrix and four columns:

- `Fluorophore`:

  Fluorophore name.

- `Reference.SSI`:

  SSI under the reference spectrum.

- `Test.SSI`:

  SSI under the test spectrum.

- `Reference.Spill`:

  Fractional spillover under the reference spectrum.

- `Test.Spill`:

  Fractional spillover under the test spectrum.

The figure is saved to `plot.dir` and the combined plot is also printed
to the active graphics device.
