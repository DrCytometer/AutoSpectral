# Reload Flow Control Information

This function reloads essential information from control files to permit
rapid unmixing at a later date, without recalculating spectra or gates.

## Usage

``` r
reload.flow.control(control.dir, control.def.file, asp)
```

## Arguments

- control.dir:

  file path to the single stained control fcs files

- control.def.file:

  csv file defining the single color control file names, fluorophores
  they represent, marker names, peak channels and gating requirements.

- asp:

  The AutoSpectral parameter list. Prepare using
  `get.autospectral.param`

## Value

A list containing the reloaded flow control information.
