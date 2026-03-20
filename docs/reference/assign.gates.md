# Assign Gates

Assigns gate labels to samples based on whether they share
characteristics, such as being the same for `viability.gate`,
`large.gate` or `control.type` in the control file and whether the
universal negative used is the same. Updates the supplied control table
and returns it with gate assignments.

## Usage

``` r
assign.gates(control.table, gating.system, gate, verbose = TRUE)
```

## Arguments

- control.table:

  Dataframe or table of the control file, read in via
  [`define.flow.control()`](https://drcytometer.github.io/AutoSpectral/reference/define.flow.control.md)
  and cleaned up by that function.

- gating.system:

  Character string selecting the automated gating system to employ in
  defining the initial scatter gates for identifying cells in the FCS
  files. Options are `landmarks` and `density`. The `density` option
  uses the original gating system from AutoSpill, picking cell
  populations based on dense regions on FSC/SSC. The default `landmarks`
  system picks out the cell regions using the brightest events in the
  peak channel for the single-stained control(s). This approach is
  generally more robust. A good way to use this is to utilize
  `landmarks` in combination with specifying which single-stained
  control files should be used for defining the gates. For example,
  abundant, bright cell markers such as CD4 will reliably identify the
  lymphocyte region, as CD14 will identify the monocyte region, etc. For
  more instructions, see the help pages on GitHub.

- gate:

  Logical, default is `TRUE`, in which case, automated gating will be
  performed. If `FALSE`, the FCS files will be imported without
  automatically generated gates applied. That is, all data in the files
  will be used. This is intended to allow the user to pre-gate the files
  in commercial software.

- verbose:

  Logical, default is `TRUE`. Set to `FALSE` to suppress messages.

## Value

A modified `control.table`, complete with gate.name and gate.define

## See also

- [`tune.gate()`](https://drcytometer.github.io/AutoSpectral/reference/tune.gate.md)

- [`define.flow.control()`](https://drcytometer.github.io/AutoSpectral/reference/define.flow.control.md)
