# Sanitize Optimization Inputs

Checks the delta list and delta norms for empty, NULL or redundant
content, removes fluorophores from the list of ones to be optimized if
they fail checks.

## Usage

``` r
sanitize.optimization.inputs(spectra, optimize.fluors, variants, delta.norms)
```

## Arguments

- spectra:

  Numeric matrix (fluors x detectors)

- optimize.fluors:

  Character vector of fluorophores present in variants

- variants:

  List of variant matrices per fluorophore

- delta.norms:

  List of delta norms per fluorophore

## Value

Updated list of fluorophores to be optimized
