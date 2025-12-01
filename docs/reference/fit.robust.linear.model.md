# Fit Robust Linear Model

Returns a matrix by rows, with the intercept and coefficient of a robust
linear model fitted to the input data. Reverts to a standard linear
model in case of no convergence.

## Usage

``` r
fit.robust.linear.model(
  x.data,
  y.data,
  x.name,
  y.name,
  max.iter = 100,
  fix.unmix = FALSE
)
```

## Arguments

- x.data:

  A vector containing the predictor variable data.

- y.data:

  A vector containing the response variable data.

- x.name:

  The name of the predictor variable.

- y.name:

  The name of the response variable.

- max.iter:

  Numeric. Maximum number of iterations for the robust linear model.
  Default is `100`.

- fix.unmix:

  Logical, default is `FALSE`. If `TRUE`, sets coefficient to zero in
  case of failed convergence. Used for `fix.my.unmix.`

## Value

A vector with the intercept and coefficient.
