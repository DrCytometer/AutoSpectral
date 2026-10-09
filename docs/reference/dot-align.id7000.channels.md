# Align ID7000 Database Channel Names to an FCS File

Internal helper. ID7000 detector names are built from a laser label,
"CH", a detector number and a suffix (e.g. "320CH12-A"). Some
instruments write a different laser label into the FCS file (e.g.
"DUVCH12-A"). The `cytometer_database.csv` therefore holds a second name
for each ID7000 detector in the `ID7000_alt` column.

For each database detector, the alternate name is used only when the
primary name (`ID7000` column) is not among `file.channels` and the
alternate name is. Otherwise the primary name is kept. The choice is
made one detector at a time, so an instrument that relabels only some
lasers is handled, and an alternate name that never occurs in the data
has no effect.

Names in `x` that are not ID7000 database names (scatter, time, blanks,
`NA`) are returned unchanged.

## Usage

``` r
.align.id7000.channels(x, file.channels, database = NULL)
```

## Arguments

- x:

  Character vector of ID7000 channel names as written in the `ID7000`
  column of `cytometer_database.csv` (or in the ID7000 spectral
  reference library).

- file.channels:

  Character vector of channel names found in the acquired data, or the
  detector names of a spectra matrix built from it.

- database:

  Data frame read from `cytometer_database.csv`. Default `NULL`, which
  reads the copy installed with AutoSpectral.

## Value

Character vector the same length and order as `x`.
