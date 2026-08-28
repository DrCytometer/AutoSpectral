# Get Autospectral Parameters for BD FACSDiscover Cytometers

Returns parameters for running a calculation of unmixing with
AutoSpectral for the BD FACSDiscover family (A8, S8), without creating
any figures or tables. A8 and S8 share identical acquisition parameters;
the caller
([`get.autospectral.param()`](https://drcytometer.github.io/AutoSpectral/reference/get.autospectral.param.md))
sets the specific `cytometer` label ("FACSDiscover A8" / "FACSDiscover
S8" / generic "FACSDiscover") after this function returns, based on
which alias was requested.

## Usage

``` r
get.autospectral.param.discover(autosp.param)
```

## Arguments

- autosp.param:

  A list of initial AutoSpectral parameters.

## Value

A list of AutoSpectral parameters specific to the FACSDiscover cytometer
family.
