# Spectral Reference Plot

This function plots a comparison between the user's fluorophore
spectrum, from `spectra`, and a known reference spectrum (if available)
for the same fluorophore. Quality control is performed using cosine
similarity between the two spectral profiles.

## Usage

``` r
spectral.reference.plot(
  spectra,
  asp,
  fluorophore = rownames(spectra),
  qc.threshold.warn = 0.98,
  qc.threshold.fail = 0.9,
  experiment.control.color = "black",
  library.reference.color = "blue",
  comparison.color = "red",
  experiment.line.type = "solid",
  library.line.type = "dotted",
  comparison.line.type = "dashed",
  pass.color = "darkgreen",
  warn.color = "darkorange",
  fail.color = "red",
  linewidth = 1,
  plot.dir = "./figure_spectra",
  filename = "spectral_qc_report.pdf",
  comparison.spectra = NULL,
  comparison.label = "Pre-refinement",
  highlight.fluors = character(0L)
)
```

## Arguments

- spectra:

  Matrix or dataframe containing spectral data. This should be in format
  fluorophores x detectors. Row names are used to label each trace (e.g.
  `sample` identifiers, which may disambiguate multiple controls for the
  same fluorophore) and column names will be used as the detectors
  (channels).

- asp:

  The AutoSpectral parameter list. Used to determine which cytometer
  produced the data.

- fluorophore:

  Optional character vector, the same length and order as
  `nrow(spectra)`, giving the true fluorophore identity of each row.
  Used to look up the library reference spectrum and any
  `comparison.spectra`. Defaults to `rownames(spectra)`, matching the
  previous behavior when each row already corresponds to a single,
  uniquely-named fluorophore.

- qc.threshold.warn:

  Numeric, default `0.98`. The similarity value to trigger a QC warning.

- qc.threshold.fail:

  Numeric, default `0.90`. The similarity value to trigger a QC failure.

- experiment.control.color:

  Color for the line representing the user's fluorophore spectrum from
  the single-stained reference control. Default is `black`.

- library.reference.color:

  Color for the line representing the library reference standard
  fluorophore spectrum. Default is `"blue"`.

- comparison.color:

  Color for the optional comparison spectrum. Default is `"red"`.

- experiment.line.type:

  Line style for the line representing the user's fluorophore spectrum
  from the single-stained reference control. Default is `"solid"`.

- library.line.type:

  Line style for the line representing the library reference standard
  fluorophore spectrum. Default is `"dotted"`.

- comparison.line.type:

  Line style for the line representing the optional comparison spectrum.
  Default is `"dashed"`.

- pass.color:

  Color to label similarity values above the `qc.threshold.warn`, i.e.,
  fluorophores passing QC. Default is `darkgreen`.

- warn.color:

  Color to label similarity values above the `qc.threshold.fail` but
  below the `qc.threshold.warn`, i.e., fluorophores in a potentially
  problematic zone. Default is `darkorange`.

- fail.color:

  Color to label similarity values below the `qc.threshold.fail`, i.e.,
  fluorophores failing QC. Default is `red`.

- linewidth:

  Width of the line for the spectral traces. Default is `1`.

- plot.dir:

  Directory where the files will be saved. Default is
  `./figure_spectra`.

- filename:

  Name for the output PDF file. Default is `spectral_qc_report.pdf`.

- comparison.spectra:

  Optional matrix of additional spectral profiles to be plotted as
  comparisons. Used for `get.spectra.automated` for legacy fallback.

- comparison.label:

  Optional character string, default `"Pre-refinement"`, for the
  `comparison.spectra`.

- highlight.fluors:

  Optional numeric for differentiating comparison fluors in the plots
  with the "refined" tag.

## Value

None. Plots are saved to a PDF in `plot.dir`.
