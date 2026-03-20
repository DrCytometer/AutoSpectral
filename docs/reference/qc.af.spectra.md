# QC Autofluorescence Spectra

This function performs a comparison between the user's fluorophore
spectra, from `spectra`, and identified autofluorescence spectra to
check for any contamination. Contamination may occur if the unstained
sample used to generate the AF spectra was not actually unstained,
contained fluorescent protein reporter constructs, or was contaminated
during acquisition, for example by crossover from nearby wells on a
plate. Quality control is performed using cosine similarity.

## Usage

``` r
qc.af.spectra(
  af.spectra,
  spectra,
  output.dir = "./figure_autofluorescence",
  remove = TRUE,
  pass = 1
)
```

## Arguments

- af.spectra:

  Matrix or dataframe containing autofluorescence spectral signatures.
  This should be in format AFn x detectors. Row names will be used as
  the names.

- spectra:

  Matrix or dataframe containing spectral data. This should be in format
  fluorophores x detectors. Row names will be used as the fluorophore
  names. Column names will be used as the detectors (channels).

- output.dir:

  Directory where the files will be saved. Default is
  `./figure_autofluorescence`.

- remove:

  Logical. When `TRUE`, identified highly similar AF spectra will be
  removed, returning only AF spectra that at least 0.995 or less from
  any fluorophore spectrum.

- pass:

  Numeric, default `1`. Counter to separate multiple passes of AF
  extraction, such as occur when `refine=TRUE` in
  [`get.af.spectra()`](https://drcytometer.github.io/AutoSpectral/reference/get.af.spectra.md).

## Value

None. Plots are saved to a PDF in `plot.dir`.

## See also

- [`af.qc.plot()`](https://drcytometer.github.io/AutoSpectral/reference/af.qc.plot.md)

- [`get.af.spectra()`](https://drcytometer.github.io/AutoSpectral/reference/get.af.spectra.md)
