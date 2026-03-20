# Check Gates For Errors

A small helper function to check for non-standard or crash-inducing
input when users supply gates to
[`define.flow.control()`](https://drcytometer.github.io/AutoSpectral/reference/define.flow.control.md).

## Usage

``` r
check.gates(gate.list, control.table, asp)
```

## Arguments

- gate.list:

  Named list of gates. To use this, pre-define the gates using
  [`define.gate.landmarks()`](https://drcytometer.github.io/AutoSpectral/reference/define.gate.landmarks.md)
  and/or
  [`define.gate.density()`](https://drcytometer.github.io/AutoSpectral/reference/define.gate.density.md),
  ensure that the names of the gates correspond to the names in the
  `control.def.file`, and ensure that the `gate.name` column has been
  filled in for the `control.def.file`.

- control.table:

  Dataframe or table of the control file, read in via
  [`define.flow.control()`](https://drcytometer.github.io/AutoSpectral/reference/define.flow.control.md)
  and cleaned up by that function.

- asp:

  The AutoSpectral parameter list defined using
  `get.autospectral.param`.

## Value

Silently returns `TRUE` if all checks pass. If any check fails, the
pipeline halts.

## See also

- [`tune.gate()`](https://drcytometer.github.io/AutoSpectral/reference/tune.gate.md)

- [`define.gate.landmarks()`](https://drcytometer.github.io/AutoSpectral/reference/define.gate.landmarks.md)

- [`define.gate.density()`](https://drcytometer.github.io/AutoSpectral/reference/define.gate.density.md)

- [`do.gate()`](https://drcytometer.github.io/AutoSpectral/reference/do.gate.md)

- [`define.flow.control()`](https://drcytometer.github.io/AutoSpectral/reference/define.flow.control.md)
