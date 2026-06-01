# Validate Control File

Checks whether the control file used to define the single-stained
control setup for AutoSpectral conforms to expectations and follows
rules required for successful running of AutoSpectral.

## Usage

``` r
validate.control.file(
  control.dir,
  control.def.file,
  asp,
  min.event.warning,
  min.event.error,
  legacy = FALSE
)
```

## Arguments

- control.dir:

  File path to the single-stained control FCS files.

- control.def.file:

  CSV file defining the single-color control file names, fluorophores
  they represent, marker names, peak channels, and gating requirements.

- asp:

  The AutoSpectral parameter list.

- min.event.warning:

  The number of events in the entire FCS file that will trigger a
  warning if not met.

- min.event.error:

  The number of events in the entire FCS file that will trigger an error
  if not met.

- legacy:

  Logical. If `FALSE`, gating-related columns will not be created and
  the control file will be suitable only for the new automated spectral
  extraction pipeline using
  [`get.spectra.automated()`](https://drcytometer.github.io/AutoSpectral/reference/get.spectra.automated.md).
  To use the version 1 "legacy" pipeline for the extraction of
  fluorophore spectra, using gating and
  [`define.flow.control()`](https://drcytometer.github.io/AutoSpectral/reference/define.flow.control.md),
  set `legacy=TRUE`.

## Value

A dataframe of errors and warnings intended to help the user fix
problems with the `control.def.file`.
