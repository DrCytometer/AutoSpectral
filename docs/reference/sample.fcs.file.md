# Sample FCS File

This function samples events from an FCS file based on specified
parameters and downsampling criteria.

## Usage

``` r
sample.fcs.file(file.name, control.dir, downsample.n, asp)
```

## Arguments

- file.name:

  A character string specifying the name of the FCS file.

- control.dir:

  A character string specifying the directory containing the control FCS
  file.

- downsample.n:

  A numeric value indicating the number of events to downsample.

- asp:

  The AutoSpectral parameter list. Prepare using
  `get.autospectral.param`

## Value

A matrix containing the sampled events from the FCS file.
