# New Issue

Helper function to log issues in control file

## Usage

``` r
.new_issue(severity, rule, filename = NA, column = NA, message)
```

## Arguments

- severity:

  Character, severity of the issue, e.g., "error", "warning".

- rule:

  Character, rule being broken.

- filename:

  Character, name of the FCS file linked to the issue.

- column:

  Character, column of the control file linked to the issue.

- message:

  Character, output message to guide the user in fixing the error.

## Value

A dataframe of issues
