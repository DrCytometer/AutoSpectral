# Assign AF Spectrum By Residual Alignment

Aligns the residuals with the autofluorescence spectral variants to
determine which best fits each cell (event). A fast approximation for
brute force sequential unmixing method in early versions of
AutoSpectral. Provides similar results to residual minimization, and
when combined with per-cell optimization, works well. Substantially
faster.

## Usage

``` r
assign.af.residuals(raw.data, spectra, af.spectra)
```

## Arguments

- raw.data:

  Expression data from raw FCS files. Cells in rows and detectors in
  columns. Columns should be fluorescent data only and must match the
  columns in spectra.

- spectra:

  Spectral signatures of fluorophores, normalized between 0 and 1, with
  fluorophores in rows and detectors in columns.

- af.spectra:

  Spectral signatures of autofluorescences, normalized between 0 and 1,
  with fluorophores in rows and detectors in columns.

## Value

Row indices for best-fitting AF spectra (from `af.spectra`)
