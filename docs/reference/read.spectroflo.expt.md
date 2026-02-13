# Extract Spectra From SpectroFlo Expt File

Reads an Experiment (.Expt) file from SpectroFlo and extracts the
spillover matrix that was used for the unmixing. Note that the
information stored in the SpectroFlo .Expt files for spillover is
sometimes perfect for unmixing and is at other times very odd. I do not
currently know why this is, but it may relate to the way the software
updates the spectral profiles using the QC information. If so, spectra
extracted using `read.spectroflo.expt()` are only likely to be accurate
when the .Expt file has been last opened on the instrument used to
acquire the sames, not an offline copy.

## Usage

``` r
read.spectroflo.expt(
  expt.file,
  output.dir,
  output.filename = "SpectroFlo_spillover_matrix.csv",
  fcs.file = NULL
)
```

## Arguments

- expt.file:

  File name and path to the .Expt file to be read.

- output.dir:

  Directory where the spillover .csv file will be written.

- output.filename:

  Name for the output spillover file. Default is
  `"SpectroFlo_spillover_matrix.csv"`.

- fcs.file:

  File name and path to an FCS file to be read. Optional, but inclusion
  of an FCS file will allow for the detector names (channels) to be
  added to the returned spillover matrix.

## Value

Spillover matrix (fluorophores x detectors). Fluorophore names will be
automatically extracted from the control .FCS files linked to the .Expt
file, with an attempt to clean them up. Detector names (columns) are not
extracted.
