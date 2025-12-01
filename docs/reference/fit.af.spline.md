# Fit Spline to Autofluorescence Data

This function fits a spline to autofluorescence data, removing extreme
events and defining bounds equally far from zero. It uses robust linear
modeling and identifies events within a specified number of standard
deviations from the spline.

## Usage

``` r
fit.af.spline(af.cells, non.af.cells, asp)
```

## Arguments

- af.cells:

  A matrix containing the autofluorescence data.

- non.af.cells:

  A matrix containing the low autofluorescence data.

- asp:

  The AutoSpectral parameter list. Prepare using
  `get.autospectral.param`

## Value

The boundary of the events within the specified number of standard
deviations from the spline.
