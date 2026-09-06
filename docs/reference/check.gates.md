# Check Gates For Errors

A small helper function to check for non-standard or crash-inducing
input when users supply gates to
[`define.flow.control()`](https://drcytometer.github.io/AutoSpectral/reference/define.flow.control.md).

## Usage

``` r
check.gates(gate.list, control.table, asp, bound.tolerance = 0.05)
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

- bound.tolerance:

  Numeric, default `0.05`. Fractional margin, relative to each axis's
  configured range (`scatter.data.max.* - scatter.data.min.*`), allowed
  outside `scatter.data.min.*`/`scatter.data.max.*` before a gate's
  coordinates are rejected.
  [`do.gate()`](https://drcytometer.github.io/AutoSpectral/reference/do.gate.md)'s
  own auto-gating search region is clamped to these same bounds rather
  than erroring when real data extends slightly past them (raw scatter
  from baseline-corrected digital detectors commonly dips a little
  negative near the origin), so a pre-defined gate is held to the same
  loose ceiling here rather than an exact boundary. Set to `0` to
  restore a hard boundary.

## Value

Silently returns `TRUE` if all checks pass. If any check fails, the
pipeline halts.

## See also

- [`tune.gate()`](https://drcytometer.github.io/AutoSpectral/reference/tune.gate.md)

- [`define.gate.landmarks()`](https://drcytometer.github.io/AutoSpectral/reference/define.gate.landmarks.md)

- [`define.gate.density()`](https://drcytometer.github.io/AutoSpectral/reference/define.gate.density.md)

- [`do.gate()`](https://drcytometer.github.io/AutoSpectral/reference/do.gate.md)

- [`define.flow.control()`](https://drcytometer.github.io/AutoSpectral/reference/define.flow.control.md)
