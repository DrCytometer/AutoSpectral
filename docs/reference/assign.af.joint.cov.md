# Assign AF Spectrum By Joint Covariance-Weighted Error

Assigns each cell to the best-fitting autofluorescence spectral variant
using a joint scoring criterion that multiplies two proportional error
terms: a covariance-weighted fluorophore error and a raw-space residual
error. The covariance of the AF spectra library is propagated into
fluorophore space via the unmixing matrix to derive per-channel error
weights, giving channels where AF variation matters most a
proportionally larger influence on the assignment decision. Multiplying
the two terms rewards variants that achieve large improvements on either
axis, without requiring an explicit mixing parameter.

More principled than `assign.af.fluorophores` (plain L1) or
`assign.af.residuals` (simple residual dot-product) when the AF library
contains spectrally diverse variants, because the covariance weights
naturally downweight channels where all AF variants look similar and
upweight channels where they diverge.

## Usage

``` r
assign.af.joint.cov(raw.data, spectra, af.spectra)
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

## Value

Integer vector of length `nrow(raw.data)` giving the row index (into
`af.spectra`) of the best-fitting AF variant for each cell.
