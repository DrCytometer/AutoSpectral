# Autofluorescence QC Plot

This function plots a comparison between the user's fluorophore
spectrum, from `spectra`, and an identified autofluorescence spectrum
that may represent contamination. Quality control is performed using
cosine similarity between the two spectral profiles.

## Usage

``` r
af.qc.plot(
  af.spectra,
  spectra,
  qc.table,
  af.color = "black",
  fluor.color = "blue",
  af.line.type = "solid",
  fluor.line.type = "dotted",
  linewidth = 1,
  plot.dir = "./figure_autofluorescence",
  filename = "autofluorescence_qc_report.pdf"
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

- qc.table:

  Dataframe containing the pre-screened problematic AF spectra and the
  fluorophores to which they are similar. Format: columes named `AF`,
  `Fluorophore` and `Similarity` containing, respectively, the names of
  the AF spectra, the names of the similar fluorophores and the cosine
  similarity values.

- af.color:

  Color for the line representing the user's autofluorescence spectrum
  from the unstained reference control. Default is `black`.

- fluor.color:

  Color for the line representing the similar fluorophore spectrum.
  Default is `blue`.

- af.line.type:

  Line style for the line representing the user's AF spectrum. Default
  is `solid`.

- fluor.line.type:

  Line style for the line representing thefluorophore spectrum. Default
  is `dotted`.

- linewidth:

  Width of the line for the spectral traces. Default is `1`.

- plot.dir:

  Directory where the files will be saved. Default is
  `./figure_autofluorescence`.

- filename:

  Name for the output PDF file. Default is
  `autofluorescence_qc_report.pdf`.

## Value

None. Plots are saved to a PDF in `plot.dir`.

## See also

- [`qc.af.spectra()`](https://drcytometer.github.io/AutoSpectral/reference/qc.af.spectra.md)

- [`get.af.spectra()`](https://drcytometer.github.io/AutoSpectral/reference/get.af.spectra.md)
