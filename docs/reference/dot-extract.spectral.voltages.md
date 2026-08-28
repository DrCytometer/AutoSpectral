# Extract Spectral Voltages/Gains

Internal helper. Given an FCS header and a keyword suffix from
[`.spectral.voltage.suffix()`](https://drcytometer.github.io/AutoSpectral/reference/dot-spectral.voltage.suffix.md),
extracts the voltage/gain value for each named spectral channel (matched
by channel name, not position).

## Usage

``` r
.extract.spectral.voltages(header, spectral.channel, suffix)
```

## Arguments

- header:

  Parsed FCS header/keyword list for the file being checked.

- spectral.channel:

  Character vector of spectral detector/channel names.

- suffix:

  Character keyword suffix (`"V"`, `"G"`, or `"R"`) from
  [`.spectral.voltage.suffix()`](https://drcytometer.github.io/AutoSpectral/reference/dot-spectral.voltage.suffix.md).
