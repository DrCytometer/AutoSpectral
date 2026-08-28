# Check Control Voltage Consistency

Internal helper. Compares detector voltage/gain settings for each
spectral channel across a set of single-stained control FCS files. The
first file in `filenames` is the reference; every other file is compared
against it, channel by channel. Returns a data frame of mismatches (zero
rows if everything is consistent, or if the cytometer/file exposes no
usable voltage/gain keyword).

## Usage

``` r
.check.control.voltages(control.dir, filenames, spectral.channel, asp)
```

## Arguments

- control.dir:

  File path to the single-stained control FCS files.

- filenames:

  Character vector of filenames for the FCS files.

- spectral.channel:

  Character vector of spectral detector/channel names.

- asp:

  The AutoSpectral parameter list
