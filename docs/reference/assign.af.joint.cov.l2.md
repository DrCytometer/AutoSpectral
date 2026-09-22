# Assign AF Spectrum By Joint Covariance-Weighted Squared Error

Assigns each cell to the best-fitting autofluorescence spectral variant
using a joint scoring criterion that multiplies two proportional error
terms: a covariance-weighted fluorophore error and a raw-space residual
error, both measured as squared (L2) deviations. The covariance of the
AF spectra library is propagated into fluorophore space via the unmixing
matrix to derive per-channel error weights, giving channels where AF
variation matters most a proportionally larger influence on the
assignment decision. Multiplying the two terms rewards variants that
achieve large improvements on either axis, without requiring an explicit
mixing parameter.

Because both error terms are quadratic in the per-cell, per-variant AF
abundance, each can be expanded into a baseline term, a cross term, and
a curvature term. All three are computed as single matrix products
across every cell and every variant simultaneously, so no per-variant
loop is required. This makes the function substantially faster than the
L1 (`abs`-based) formulation in `assign.af.joint.cov`, at the cost of
being somewhat less robust to outlier channels, since squared error
weights large deviations more heavily than L1.

## Usage

``` r
assign.af.joint.cov.l2(raw.data, spectra, af.spectra, return.scores = FALSE)
```

## Arguments

- raw.data:

  Expression data from raw FCS files. Cells in rows and detectors in
  columns. Columns should be fluorescent data only and must match the
  columns in `spectra`.

- spectra:

  Spectral signatures of fluorophores, normalized between 0 and 1, with
  fluorophores in rows and detectors in columns.

- af.spectra:

  Spectral signatures of autofluorescences, normalized between 0 and 1,
  with AF variants in rows and detectors in columns. Prepare using
  `get.af.spectra`.

- return.scores:

  Logical, default `FALSE`. If `\code{TRUE}`, also returns the unmixed
  data and scores for each AF variant per cell.

## Value

Integer vector of length `nrow(raw.data)` giving the row index (into
`af.spectra`) of the best-fitting AF variant for each cell.
