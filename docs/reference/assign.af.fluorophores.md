# Assign AF Spectrum By Fluorophore Projection

Projects the autofluorescence spectral variants into fluorophore unmixed
space to determine which best fits each cell (event). A fast
approximation for brute force sequential unmixing method in early
versions of AutoSpectral. Provides essentially identical results to
minimization of fluorophore signal (dist0 method). Substantially faster.
Use L1 (absolute value) minimization, which works better than the
standard L2 (squared error).

## Usage

``` r
assign.af.fluorophores(raw.data, spectra, af.spectra)
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
