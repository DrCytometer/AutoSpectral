# Spectral Voltage/Gain Keyword Suffix

Internal helper. Determines which FCS keyword suffix (\$PnV, \$PnG, or
\$PnR) should be compared as the "voltage" for a given cytometer.

- Mosaic uses \$PnG (gain) instead of \$PnV.

- ID7000 has no per-channel PMT voltage.

- All other supported cytometers use \$PnV.

## Usage

``` r
.spectral.voltage.suffix(asp, header)
```

## Arguments

- asp:

  The AutoSpectral parameter list; `asp$cytometer` is used.

- header:

  Parsed FCS header/keyword list for the file being checked.
