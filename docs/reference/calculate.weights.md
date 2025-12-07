# Calculate Weights

This function calculates weights for generating a weighted least-squares
(WLS) unmixing matrix. It does so using mean expression levels in the
provided FCS file. The provided FCS file should thus be representative
of the data to which the unmixing matrix will be applied. One option
would be to concatenate samples of all FCS files to be unmixed and use
that as input to `fcs.file`.

## Usage

``` r
calculate.weights(
  fcs.file,
  spectral.channels,
  save = FALSE,
  output.dir = "./table_spectra",
  filename = "weights.csv"
)
```

## Arguments

- fcs.file:

  A character string specifying the path to the FCS file.

- spectral.channels:

  Character string specifying the names of the spectral detectors in the
  FCS file. Can be obtained from `flow.control$spectral.channel` or
  `colnames(spectra)`.

- save:

  Logical, if `TRUE`, save the weights to a CSV file in `output.dir`.

- output.dir:

  Character string specifying the directory to save the CSV file, if
  generated.

- filename:

  Character string specifying the filename for the CSV file. Default is
  `weights.csv`.

## Value

A named numeric vector with weighting for each detector.
