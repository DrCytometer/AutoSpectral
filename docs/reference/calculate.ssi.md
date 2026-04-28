# Calculate Secondary Stain Index and Spillover

Calculates the Secondary Stain Index (SSI) and secondary spillover for
every fluorophore channel in an unmixed dataset, given a single-stained
control and a matched unstained control. SSI is defined as:

\$\$SSI = \frac{median(pos) - median(neg)}{2 \times rSD(neg)}\$\$

where *rSD* is the robust standard deviation (median absolute deviation,
MAD). Secondary spillover is defined as:

\$\$Spillover = \frac{\Delta MFI\_{channel}}{\Delta
MFI\_{on\text{-}channel}}\$\$

Positive events are identified as single-stained cells whose intensity
in the on-target channel exceeds the 99th percentile of the unstained
control. Up to 500 of the brightest positive events and up to 500
randomly sampled unstained events are used for robustness and speed.

## Usage

``` r
calculate.ssi(unstained.unmixed, single.stained.unmixed, fluor.name)
```

## Arguments

- unstained.unmixed:

  A numeric matrix of unmixed fluorescence values for the unstained
  control. Rows are cells; columns are fluorophore channels.

- single.stained.unmixed:

  A numeric matrix of unmixed fluorescence values for the single-stained
  control. Must have the same column names as `unstained.unmixed`.

- fluor.name:

  Character string. The column name of the target fluorophore (the
  "on-channel") used to define positive events and calculate \\\Delta
  MFI\_{on\text{-}channel}\\.

## Value

A data frame with one row per fluorophore channel (matching the column
names of `single.stained.unmixed`) and two columns:

- `Spillover`:

  Fractional secondary spillover into each channel.

- `SSI`:

  Secondary Stain Index for each channel. The SSI for `fluor.name`
  itself is set to `0` (on-channel reference).
