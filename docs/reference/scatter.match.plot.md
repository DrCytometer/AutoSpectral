# Plot Scatter-Matching of Universal Negative

This function generates scatter plots for matching positive and negative
expression data, selected based on scatter parameters gates.

## Usage

``` r
scatter.match.plot(
  pos.expr.data,
  neg.expr.data,
  fluor.name,
  scatter.param,
  asp,
  color.palette = "rainbow"
)
```

## Arguments

- pos.expr.data:

  A matrix containing the positive expression data.

- neg.expr.data:

  A matrix containing the negative expression data.

- fluor.name:

  A character string specifying the fluorophore name.

- scatter.param:

  A character vector specifying the scatter parameters.

- asp:

  The AutoSpectral parameter list. Prepare using
  `get.autospectral.param`

- color.palette:

  Optional character string defining the viridis color palette to be
  used for the fluorophore traces. Default is `rainbow`, which will be
  similar to FlowJo or SpectroFlo. Other pptions are the viridis color
  options: `magma`, `inferno`, `plasma`, `viridis`, `cividis`, `rocket`,
  `mako` and `turbo`.

## Value

None. The function saves the generated scatter plot to a file.
