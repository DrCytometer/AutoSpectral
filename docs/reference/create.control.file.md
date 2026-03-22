# Create Control File

A helper function to draft a description of your single stained control
files such that AutoSpectral can understand and process them correctly.
Given a set of single stained control fcs files, `create.control.file`
will produce a csv file listing the matching peak detector channels for
your fluorophores (if known). If your files contain bead or cell tags in
the filename, it will assign your controls as cells or beads. You will
need to fill in any "No Match" results manually. You will need to set
universal negatives manually. You will need to add marker names
manually.

## Usage

``` r
create.control.file(
  control.dir,
  asp,
  fill.gate.name = TRUE,
  filename = "fcs_control_file"
)
```

## Arguments

- control.dir:

  file path to the single stained control fcs files

- asp:

  The AutoSpectral parameter list. Generate using
  `get.autospectral.param`

- fill.gate.name:

  Logical, default is `TRUE`. Will attempt to automatically assign gate
  names for the `gate.name` column if `TRUE`.

- filename:

  Character string defining the output filename. Default is
  "fcs_control_file", to which .csv will be appended.

## Value

No returns. Outputs a csv file called fcs_control_file.csv
