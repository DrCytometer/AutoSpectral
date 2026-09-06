# Check GitHub for a newer tagged release of a package

Queries the GitHub tags API for `repo`, compares the newest semantic-
version tag against `installed.version`, and returns that version if it
is newer. Silently returns NULL on any failure (no network, rate limit,
unparseable response, non-interactive session) since this is a
convenience check and must never interrupt or slow down a script.
Results are cached on disk for `cache.days` so a session does not
re-query GitHub every time the package is attached.

## Usage

``` r
.AS.check.github.update(pkg, repo, installed.version, cache.days = 1)
```

## Arguments

- pkg:

  Package name, used as the cache subdirectory.

- repo:

  GitHub "owner/repo" string.

- installed.version:

  A `package_version` to compare against.

- cache.days:

  Minimum days between live checks. Default 1.

## Value

A `package_version` if a newer tag exists, otherwise NULL (invisibly).
