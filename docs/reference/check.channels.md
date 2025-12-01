# Check Channels

A helper function to reorganize the spectral channels in a nice order
for plotting. Puts them in excitation/emission order.

## Usage

``` r
check.channels(spectral.channels, asp)
```

## Arguments

- spectral.channels:

  Vector of initial spectral channel names.

- asp:

  The AutoSpectral parameter list. Generate using
  `get.autospectral.param`

## Value

Returns the vector of spectral channels re-organized in excitation-
emission order. That is, narrowest to longest excitation laser, wtih
narrowest to longest emission wavelength inside each laser group.
