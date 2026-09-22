# Select Retained Channels

Internal helper for
[`unmix.fcs()`](https://drcytometer.github.io/AutoSpectral/reference/unmix.fcs.md).
Selects the parameters of a raw FCS file that are carried through
unchanged into the unmixed file, alongside the newly unmixed fluorophore
parameters.

Some instruments (e.g. the BD FACSDiscover family and the FACSymphony A5
SE) write "raw" FCS files that also hold the instrument's own unmixed
fluorophore parameters. These are never carried through: they would
duplicate, and could be mistaken for, the AutoSpectral output. A
parameter is retained only if it matches one of the named parameter
families in `asp$non.spectral.channel` (time, scatter, imaging and other
acquisition parameters). Bare suffix entries in that vector (`-H`, `-W`,
`-T`) exist to keep detector companions out of the spectral channel set
and are not used to retain anything, since fluorophore names such as
`APC-H7` would otherwise match them.

Raw detector data are written only when `include.raw = TRUE`. For each
detector used in unmixing this is the `-A` (or `-H`) parameter that was
unmixed, plus its `-T` companion where the file has one. The `-H` and
`-W` companions of the detectors are never written.

## Usage

``` r
.select.retained.channels(
  original.param,
  spectral.channel,
  asp,
  include.raw = FALSE,
  include.imaging = TRUE,
  verbose = TRUE
)
```

## Arguments

- original.param:

  Character vector of the parameter names (`$PnN`) of the raw FCS file,
  in file order.

- spectral.channel:

  Character vector of the detector channels used for unmixing (the
  column names of the spectra matrix).

- asp:

  The AutoSpectral parameter list. Prepare using
  `get.autospectral.param`.

- include.raw:

  Logical, whether to retain the raw detector data. Default is `FALSE`.

- include.imaging:

  Logical, whether to retain imaging parameters. Only applies to the BD
  FACSDiscover family, where `FALSE` restricts the retained parameters
  to time and scatter. Default is `TRUE`.

- verbose:

  Logical, whether to report the parameters that are not carried
  through. Default is `TRUE`.

## Value

Character vector of the parameter names to carry through, in the order
they are to be written.
