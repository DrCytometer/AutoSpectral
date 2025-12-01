# Check Control File

Attempts to find potential failure points and input errors in the
control file `control.def.file` used to define the single-stained
control setup for AutoSpectral.

## Usage

``` r
check.control.file(control.dir, control.def.file, asp, strict = FALSE)
```

## Arguments

- control.dir:

  File path to the single-stained control FCS files.

- control.def.file:

  CSV file defining the single-color control file names, fluorophores
  they represent, marker names, peak channels, and gating requirements.

- asp:

  The AutoSpectral parameter list defined using
  `get.autospectral.param`.

- strict:

  Logical. Controls whether the function triggers a break or continues
  and outputs a list of errors. Default is `FALSE`.

## Value

A named list of errors and warnings intended to help the user fix
problems with the `control.def.file`.
