# Handle and Plot Gating Failures

A helper function that dictates what occurs when a gating error is
detected.

## Usage

``` r
handle.gating.error(
  e,
  gate.id,
  files.to.gate,
  scatter.coords,
  samp,
  viability.gate,
  control.type,
  flow.scatter.and.channel.label,
  asp
)
```

## Arguments

- e:

  The error object from tryCatch

- gate.id:

  The ID of the gate that failed

- files.to.gate:

  Character vector of filenames used

- scatter.coords:

  The data matrix used for gating

- samp:

  String identifier for the sample/gate type

- viability.gate:

  Logical, is this a viability gate

- control.type:

  "beads" or "cells"

- flow.scatter.and.channel.label:

  Mapping vector for axis labels

- asp:

  The parameter list

## Value

None. Stops flow.
