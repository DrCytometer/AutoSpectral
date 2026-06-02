# Benchmark Scatter-Match Against Existing AF Assignment Methods

Convenience wrapper that runs `assign.af.scatter.match` alongside the
three existing methods (`assign.af.fluorophores`, `assign.af.residuals`,
`assign.af.joint.cov`) on the same test and reference unstained data.
For each existing method, the assigned AF variant spectrum is looked up
and its cosine similarity to each test cell is computed, allowing direct
apples-to-apples comparison with the scatter-match approach.

## Usage

``` r
benchmark.af.scatter.match(
  test.data,
  ref.data,
  scatter.param,
  spectra,
  af.spectra,
  k = 5L,
  verbose = TRUE
)
```

## Arguments

- test.data:

  Numeric matrix or FCS file path. Test unstained data (cells x
  channels).

- ref.data:

  Numeric matrix or FCS file path. Reference unstained data (cells x
  channels).

- scatter.param:

  Character vector of scatter channel names.

- spectra:

  Numeric matrix. Fluorophore spectra (fluorophores x detectors), as
  used by the existing assign.af.\* functions.

- af.spectra:

  Numeric matrix. AF variant spectra (variants x detectors).

- k:

  Integer. Neighbours for scatter-matching. Default `5`.

- verbose:

  Logical. Default `TRUE`.

## Value

A list with:

- `scatter.match`:

  Full output of `assign.af.scatter.match`.

- `comparison`:

  Data frame with one row per method and columns: `method`,
  `mean.cosine`, `median.cosine`, `sd.cosine`, `pct.above.0.9`,
  `pct.above.0.95`.

- `per.cell`:

  Data frame with one row per test cell containing cosine similarities
  from all four methods, for cell-level analysis.
