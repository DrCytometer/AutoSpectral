# Read Channel Information

This function reads channel information from control files and corrects
channel names based on specified forbidden characters.

## Usage

``` r
read.channel(control.dir, control.def.file, asp)
```

## Arguments

- control.dir:

  Directory containing control files.

- control.def.file:

  File containing control definitions.

- asp:

  The AutoSpectral parameter list. Prepare using
  `get.autospectral.param`

## Value

A data frame containing the original and corrected channel names.
