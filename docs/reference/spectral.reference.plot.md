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
  qc.threshold.warn = 0.98,
  qc.threshold.fail = 0.9,
  experiment.control.color = "black",
  library.reference.color = "blue",
  experiment.line.type = "solid",
  library.line.type = "dotted",
  pass.color = "darkgreen",
  warn.color = "darkorange",
  fail.color = "red",
  linewidth = 1,
  plot.dir = "./figure_spectra",
  filename = "spectral_qc_report.pdf"
)
```

## Arguments

- spectra:

  Matrix or dataframe containing spectral data. This should be in format
  fluorophores x detectors. Row names will be used as the fluorophore
  names. Column names will be used as the detectors (channels).

- asp:

  The AutoSpectral parameter list. Used to determine which cytometer
  produced the data.

- qc.threshold.warn:

  Numeric, default `0.98`. The similarity value to trigger a QC warning.

- qc.threshold.fail:

  Numeric, default `0.90`. The similarity value to trigger a QC failure.

- experiment.control.color:

  Color for the line representing the user's fluorophore spectrum from
  the single-stained reference control. Default is `black`.

- library.reference.color:

  Color for the line representing the library reference standard
  fluorophore spectrum. Default is `blue`.

- experiment.line.type:

  Line style for the line representing the user's fluorophore spectrum
  from the single-stained reference control. Default is `solid`.

- library.line.type:

  Line style for the line representing the library reference standard
  fluorophore spectrum. Default is `dotted`.

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

## Value

None. Plots are saved to a PDF in `plot.dir`.
