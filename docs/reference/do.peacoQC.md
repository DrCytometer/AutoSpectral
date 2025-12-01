# Perform Quality Control on Flow Cytometry Data using PeacoQC

This function performs quality control on flow cytometry data using the
PeacoQC method. It transforms the data, removes margins, and identifies
good cells. The function also optionally plots the results and saves the
cleaned data.

## Usage

``` r
do.peacoQC(
  dirty.expr,
  sample.name,
  spectral.channel,
  biexp.transform,
  transform.inv,
  output.dir,
  time.param,
  all.channels,
  method = "MAD",
  figures = TRUE,
  verbose = TRUE
)
```

## Arguments

- dirty.expr:

  A matrix containing the raw expression data, pre-cleaning.

- sample.name:

  The name of the sample.

- spectral.channel:

  The spectral channels to be used.

- biexp.transform:

  The biexponential transformation function.

- transform.inv:

  The inverse transformation function.

- output.dir:

  The directory to save output files.

- time.param:

  The time channel parameter.

- all.channels:

  A vector of all channels to be included in the final cleaned data.

- method:

  The PeacoQC method to use. Inherited from `PeacoQC`. Options are
  `all`, `MAD` or `IT`. Default is `MAD`.

- figures:

  Logical. Controls whether PeacoQC plots are created. Default is
  `TRUE`.

- verbose:

  Logical, default is `TRUE`. Set to `FALSE` to suppress messages.

## Value

A matrix with the cleaned expression data.
