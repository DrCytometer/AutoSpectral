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
  min.event.error
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

## Value

A dataframe of errors and warnings intended to help the user fix
problems with the `control.def.file`.
