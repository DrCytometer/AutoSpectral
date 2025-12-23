# Validate Control File

Checks whether the control file used to define the single-stained
control setup for AutoSpectral conforms to expectations and follows
rules required for successful running of AutoSpectral.

## Usage

``` r
validate.control.file(control.dir, control.def.file, asp)
```

## Arguments

- control.dir:

  File path to the single-stained control FCS files.

- control.def.file:

  CSV file defining the single-color control file names, fluorophores
  they represent, marker names, peak channels, and gating requirements.

- asp:

  The AutoSpectral parameter list.

## Value

A dataframe of errors and warnings intended to help the user fix
problems with the `control.def.file`.
