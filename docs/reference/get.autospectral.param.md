# Get AutoSpectral Parameters

Retrieves autospectral parameters for a specified cytometer.

## Usage

``` r
get.autospectral.param(cytometer = "aurora", figures = TRUE)
```

## Arguments

- cytometer:

  The type of cytometer, default is `aurora`. Supported options include
  `aurora`, `auroraNL` for Northern Lights, `id7000`, `a8`, `s8`,
  `a5se`, `opteon`, `mosaic` and `xenith`.

- figures:

  Logical indicating whether to set up directory parameters for figures
  and tables, default is `TRUE`

## Value

A list of AutoSpectral parameters.
