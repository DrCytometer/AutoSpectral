# Get Top (Brightest) Events

Retrieves the brightest `n` events in the specified FCS file in the
specified `channel`.

## Usage

``` r
get.top.events(fcs.path, channel, n = 2000, scatter.param)
```

## Arguments

- fcs.path:

  File path and name for the FCS file to be read.

- channel:

  Channel (detector) name for the selection of the brightest events.

- n:

  Number of brightest events to retrieve, default is `2000`.

- scatter.param:

  The names of the scatter channels to define which data are returned.

## Value

Returns a matrix of scatter (FSC, SSC) data
