# Fast, uncapped biexponential (logicle) transform

A self-contained implementation of the logicle transform (Parks,
Roederer & Moore, Cytometry A 2006), used as a replacement for
[`flowWorkspace::flowjo_biexp()`](https://rdrr.io/pkg/flowWorkspace/man/flowjo_biexp.html)
with the same `channelRange`/`maxValue`/`pos`/`neg`/`widthBasis`
argument names, but with no hard-coded `-1000` floor on `widthBasis`.
Uses only base R
([`stats::uniroot()`](https://rdrr.io/r/stats/uniroot.html) for a
one-time scalar solve).

## Usage

``` r
biexp.transform(
  channelRange = 4096,
  maxValue = 262144,
  pos = 4.5,
  neg = 0,
  widthBasis = -10,
  inverse = FALSE,
  newton.iter = 5L
)
```

## Arguments

- channelRange:

  numeric. Maximum value of the transformed (display) scale. Default
  `4096`.

- maxValue:

  numeric. Maximum value of the input (raw) scale. Default `262144`.

- pos:

  numeric. Number of positive decades spanned by the transform (`M` in
  the logicle parameterisation). Default `4.5`.

- neg:

  numeric. Must be `0`; see Details.

- widthBasis:

  numeric, negative. FlowJo-style width parameter; see Details for its
  relationship to the logicle `W`. Default `-10`.

- inverse:

  logical. If `TRUE`, returns the inverse transform (display scale back
  to raw values) instead of the forward transform.

- newton.iter:

  integer. Fixed number of vectorised Newton iterations used to solve
  the forward transform. Default `5` – the smallest value that reaches
  double-precision round-trip error across the ranges swept in testing.
  Has no effect on the inverse direction, which is closed-form. If a
  different `pos`/`widthBasis` combination needs more, the built-in
  residual check will warn rather than silently under-converge.

## Value

A function mapping a numeric vector on the input scale to the
transformed display scale (or the reverse, if `inverse = TRUE`).

## Details

FlowJo's `widthBasis` is a legacy re-parameterisation of the logicle
transform's width parameter `W` (in asymptotic decades): \$\$W =
\log\_{10}(-\mathrm{widthBasis}) / 2\$\$ The logicle transform requires
\\0 \< W \le M/2\\, where `M` is `pos`. FlowJo's UI, and
[`flowWorkspace::flowjo_biexp()`](https://rdrr.io/pkg/flowWorkspace/man/flowjo_biexp.html),
impose a hard floor of `widthBasis = -1000` (`W = 1.5`) regardless of
`pos` – a historical GUI choice, not a mathematical one. This function
checks the real constraint against whatever `pos` you supply and errors
with the exact numbers if it's violated.

Only `neg = 0` is supported (the value used by every cytometer profile
in this package). Extending to `neg > 0` requires the shifted-origin
form of the transform (GatingML 2.0's `A` parameter) and is deliberately
left unimplemented until there's a concrete need for it and data to
verify it against.

The forward direction (raw intensity -\> display decade) has no closed
form and is solved by Newton-Raphson, applied to the whole input vector
at once for a fixed number of iterations (not a per-point loop),
starting from `asinh(x / (2 * a.scale))` – correct in the large-`|x|`
limit and reasonable near zero. The inverse direction (display decade
-\> raw intensity) is closed-form.

## References

Parks DR, Roederer M, Moore WA (2006). A new "Logicle" display method
avoids deceptive effects of logarithmic scaling for low signals and
compensated data. Cytometry A, 69(6):541-551.
