# Apply Gate

Applies a previously-defined scatter gate boundary to a set of flow
cytometry expression data and returns only the events falling inside it.
This is a lightweight, standalone counterpart to the gating step used
internally by
[`get.gated.flow.expression.data()`](https://drcytometer.github.io/AutoSpectral/reference/get.gated.flow.expression.data.md),
intended for interactive use and for testing other `AutoSpectral`
functions against a gated subset of data without going through the full
control-file/FCS-reading pipeline.

## Usage

``` r
apply.gate(
  flow.data,
  gate.boundary,
  scatter.param = asp$default.scatter.parameter,
  asp = NULL,
  min.fraction = 0.01
)
```

## Arguments

- flow.data:

  A matrix or data frame of flow cytometry data (for example, unmixed or
  raw expression data) with named columns, including the two scatter
  parameters named in `scatter.param`.

- gate.boundary:

  A gate boundary, as returned by
  [`define.gate.landmarks()`](https://drcytometer.github.io/AutoSpectral/reference/define.gate.landmarks.md),
  [`define.gate.density()`](https://drcytometer.github.io/AutoSpectral/reference/define.gate.density.md),
  or
  [`do.gate()`](https://drcytometer.github.io/AutoSpectral/reference/do.gate.md)
  — a list containing at least numeric `x` and `y` components describing
  the polygon vertices.

- scatter.param:

  Character vector of length 2 giving the names of the two scatter
  columns in `flow.data` to gate on. Default is
  `asp$default.scatter.parameter`.

- asp:

  The AutoSpectral parameter list, prepared using
  [`get.autospectral.param()`](https://drcytometer.github.io/AutoSpectral/reference/get.autospectral.param.md).
  Only used to supply the default for `scatter.param`; not required if
  `scatter.param` is supplied directly.

- min.fraction:

  Numeric between `0` and `1`, default `0.01`. If the fraction of events
  retained by the gate falls below this value, a warning is issued (the
  function still returns the gated data). Set to `0` to disable this
  check.

## Value

`flow.data`, subset to only those events falling inside `gate.boundary`.

## See also

- [`define.gate.landmarks()`](https://drcytometer.github.io/AutoSpectral/reference/define.gate.landmarks.md)

- [`define.gate.density()`](https://drcytometer.github.io/AutoSpectral/reference/define.gate.density.md)

- [`do.gate()`](https://drcytometer.github.io/AutoSpectral/reference/do.gate.md)

- [`check.gates()`](https://drcytometer.github.io/AutoSpectral/reference/check.gates.md)

- [`get.gated.flow.expression.data()`](https://drcytometer.github.io/AutoSpectral/reference/get.gated.flow.expression.data.md)
