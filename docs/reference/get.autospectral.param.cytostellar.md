# Get Autospectral Parameters for CytoStellar Cytometer

Returns parameters for running a calculation of unmixing with
AutoSpectral for the CytoStellar, without creating any figures or
tables.

The CytoStellar allows users to acquire spectral parameters as `-A`
(area), `-H` (height), or both. `AutoSpectral` prefers `-A` when both
are present; the specific suffix actually present in a given FCS file is
resolved at read time by `.resolve.cytostellar.suffix()`, not fixed
here. Detector names in `cytometer_database.csv` and
`fluorophore_database.csv` are therefore stored WITHOUT a suffix for
this cytometer.

## Usage

``` r
get.autospectral.param.cytostellar(autosp.param)
```

## Arguments

- autosp.param:

  A list of initial AutoSpectral parameters.

## Value

A list of AutoSpectral parameters specific to the CytoStellar cytometer.
