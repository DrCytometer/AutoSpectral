# Perform Critical Check

Checks condition, if not true prints error message and stops execution.

## Usage

``` r
check.critical(condition, error.msg)
```

## Arguments

- condition:

  Essential condition(s) evaluated as a logical TRUE/FALSE.

- error.msg:

  Message to be returned to the user if condition is FALSE.

## Value

Returns error.msg and breaks if FALSE, no return if TRUE.
