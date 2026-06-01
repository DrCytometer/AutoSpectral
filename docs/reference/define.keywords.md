# Define Keywords

Updates FCS file keywords after unmixing to define the new parameters.
Tracks existing keywords from the input FCS file for metadata
compatibility.

Per [Data File Standard for Flow Cytometry, version FCS
3.1](https://doi.org/10.1002/cyto.a.20825):

- \$LAST_MODIFIED/dd-mmm-yyyy hh:mm:ss.cc

- \$LAST_MODIFIER/string/

- \$ORIGINALITY/string

are added as additional keywords to denote modification of the FCS file.

## Usage

``` r
define.keywords(
  fcs.keywords,
  final.matrix,
  spectra,
  af.spectra,
  flow.control,
  asp,
  method,
  file.name,
  weights = NULL,
  spectral.channel = NULL,
  include.imaging = FALSE
)
```

## Arguments

- fcs.keywords:

  The input keywords obtained from the read FCS file.

- final.matrix:

  The expression data, containing both unmixed data and any retained
  parameters such as scatter and time.

- spectra:

  A matrix containing the spectral data.

- af.spectra:

  Spectral signatures of autofluorescences, normalized between 0 and 1,
  with fluorophores in rows and detectors in columns.

- flow.control:

  A list containing flow cytometry control parameters.

- asp:

  The AutoSpectral parameter list.

- method:

  A character string specifying the unmixing method used.

- file.name:

  The name of the FCS file to be written.

- weights:

  Optional numeric vector of weights (one per fluorescent detector).

- spectral.channel:

  Optional character vector of the channels used for unmixing. Should
  match `weights` in length (one weight per channel).

- include.imaging:

  A logical value indicating whether to include imaging parameters in
  the written FCS file. Default is `FALSE` to provide smaller output
  files.

## Value

The updated keyword list for writing the FCS file.
