# Check Gate Consistency

A small helper function to check for non-unique input factors when
setting gate definitions.

## Usage

``` r
check.consistency(vec, var.name)
```

## Arguments

- vec:

  Vector of a factor, such as `large.gate` settings.

- var.name:

  Name of the factor (variable) being checked for multiplicity

## Value

The unique element in the factor `vec`, if all elements of `vec` are
identical.

## See also

- [`tune.gate()`](https://drcytometer.github.io/AutoSpectral/reference/tune.gate.md)

- [`define.gate.landmarks()`](https://drcytometer.github.io/AutoSpectral/reference/define.gate.landmarks.md)

- [`define.gate.density()`](https://drcytometer.github.io/AutoSpectral/reference/define.gate.density.md)

- [`do.gate()`](https://drcytometer.github.io/AutoSpectral/reference/do.gate.md)

- [`define.flow.control()`](https://drcytometer.github.io/AutoSpectral/reference/define.flow.control.md)
