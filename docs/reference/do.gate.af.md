# Perform Gating on Autofluorescence Parameters

This function returns a vector with the indexes of events inside the
initial gate on autofluorescence parameters. It proceeds through several
steps to define the gate boundaries and identify density maxima using
numerical search and Voronoi tessellations.

## Usage

``` r
do.gate.af(gate.data, samp, asp, intermediate.figures = FALSE)
```

## Arguments

- gate.data:

  A data frame containing the gate data.

- samp:

  A sample identifier.

- asp:

  The AutoSpectral parameter list. Prepare using
  `get.autospectral.param`

- intermediate.figures:

  Logical, if `TRUE` returns additional figures to show the inner
  workings of the cleaning, including definition of low-AF cell gates on
  the PCA-unmixed unstained and spectral ribbon plots of the AF
  exclusion from the unstained.

## Value

A vector with the indexes of events inside the initial gate.
