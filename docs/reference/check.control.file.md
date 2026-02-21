# Check Control File

Attempts to find potential failure points and input errors in the
control file `control.def.file` used to define the single-stained
control setup for AutoSpectral.

## Usage

``` r
check.control.file(
  control.dir,
  control.def.file,
  asp,
  strict = FALSE,
  min.event.warning = 5000,
  min.event.error = 1000
)
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

- min.event.warning:

  The number of events in the entire FCS file that will trigger a
  warning if not met. Default is `5000`.

- min.event.error:

  The number of events in the entire FCS file that will trigger an error
  if not met. Default is `1000`.

## Value

A dataframe of errors and warnings intended to help the user fix
problems with the `control.def.file`.
