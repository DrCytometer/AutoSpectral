# Get Fluorophore Spectra

This function retrieves the fluorophore spectra for flow cytometry data,
optionally using cleaned expression data. It also plots and saves the
spectra, and performs cosine similarity QC for controls.

## Usage

``` r
get.fluorophore.spectra(
  flow.control,
  asp,
  use.clean.expr = TRUE,
  af.spectra = NULL,
  title = NULL
)
```

## Arguments

- flow.control:

  A list containing flow cytometry control data.

- asp:

  The AutoSpectral parameter list. Prepare using
  `get.autospectral.param`

- use.clean.expr:

  Logical indicating whether to use cleaned expression data, default is
  `TRUE`

- af.spectra:

  Optional autofluorescence spectra to include.

- title:

  Optional prefix for plot titles, default is `NULL`, which gives
  "Initial" when `use.clean.expr` is `FALSE` and "Clean" when
  `use.clean.expr` is `TRUE`.

## Value

A matrix with the fluorophore spectra.
