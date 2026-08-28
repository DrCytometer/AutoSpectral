# Get AutoSpectral Parameters

Retrieves autospectral parameters for a specified cytometer.

## Usage

``` r
get.autospectral.param(cytometer = "aurora", figures = TRUE)
```

## Arguments

- cytometer:

  The type of cytometer, default is `aurora`. Supported options include
  `aurora`, `auroraNL` for Northern Lights, `id7000`, `discover` (BD
  FACSDiscover family), `a8` and `s8` (specific FACSDiscover models),
  `a5se`, `opteon`, `mosaic`, `xenith` and `cytostellar`. Matching is
  case-insensitive and supports unambiguous partial matches (e.g.
  `"aur"` matches `"aurora"`).

- figures:

  Logical indicating whether to set up directory parameters for figures
  and tables, default is `TRUE`

## Value

A list of AutoSpectral parameters.
