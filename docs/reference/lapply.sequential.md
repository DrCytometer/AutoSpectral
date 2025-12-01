# Sequential lapply

This function allows switching between `lapply` and `future_lapply`
without `future.seed` issues by using `lapply` sequentially.

## Usage

``` r
lapply.sequential(X, FUN, ..., future.seed = NULL)
```

## Arguments

- X:

  A vector (atomic or list) or an expression object. Other objects
  (including classed objects) will be coerced by `as.list`.

- FUN:

  The function to be applied to each element of `X`.

- ...:

  Optional arguments to `FUN`.

- future.seed:

  Ignored in this function, included for compatibility with
  `future_lapply`.

## Value

A list of the same length as `X` and named by `X`.
