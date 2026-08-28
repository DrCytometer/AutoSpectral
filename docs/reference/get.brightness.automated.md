# Get Brightness (Automated)

Fast, non-interactive estimate of fluorophore brightness (background-
subtracted MFI) for every single-stained control listed in
`control.def.file`. Unlike `get.brightness()`, no scatter gate is drawn:
brightness is measured directly from the whole acquired file using the
top `n.cells` events in each fluorophore's peak detector, with
background subtracted using the matched universal-negative control.

The peak detector for each fluorophore is taken from `spectra` (the
output of
[`get.spectra.automated()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectra.automated.md))
rather than the nominal `channel` column in the control table, so
brightness is always measured on the same channel the extracted spectrum
is actually built from, including any legacy refinement.

## Usage

``` r
get.brightness.automated(
  control.dir,
  control.def.file,
  asp,
  spectra,
  n.cells = 500L,
  verbose = TRUE
)
```

## Arguments

- control.dir:

  Character. Path to the single-stained control FCS files.

- control.def.file:

  Character. Path to the control definition CSV.

- asp:

  The AutoSpectral parameter list from
  [`get.autospectral.param()`](https://drcytometer.github.io/AutoSpectral/reference/get.autospectral.param.md).
  Currently unused directly but accepted for interface consistency with
  [`get.spectra.automated()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectra.automated.md)
  /
  [`get.spectral.variants()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectral.variants.md),
  and to allow future use (e.g. saturation thresholds).

- spectra:

  Numeric matrix as returned by
  [`get.spectra.automated()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectra.automated.md);
  fluorophores in rows (may include `"AF"`, which is ignored), detectors
  in columns.

- n.cells:

  Integer, default `500L`. Number of brightest events (in the
  fluorophore's peak channel) averaged to estimate the positive MFI.

- verbose:

  Logical, default `TRUE`.

## Value

A one-column numeric matrix of background-subtracted MFI values;
rownames are fluorophore names, column name `"MFI"`.
