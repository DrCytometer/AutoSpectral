# plot_unmix_comparison_biplots.R
#
# Per-operator representative biplots for the unmixing comparison figure:
# for a chosen on-target fluorophore's single-stained control, reads the
# already-unmixed FCS file for each operator/method folder and plots it
# against one off-target fluorophore, using the channel mapping already
# resolved by `setup.unmix.comparison()` (which itself uses
# `match.fluorophores()`), rather than re-matching fluorophore names here.
# This is the same channel-resolution problem `resolve.fluorophore.channels()`
# solves for the Figure 2-style scripts, already solved for this workflow by
# `setup.unmix.comparison()`'s own per-folder setup CSV, so no further
# fluorophore matching is done in this file.

#' @title Build One Operator's Representative Biplot
#'
#' @description
#' Reads a single operator/method's already-unmixed single-stained control
#' for `on.target.fluor` (using the filename and channel mapping already
#' resolved in that folder's `setup.unmix.comparison()` CSV), and plots it
#' against the matched channel for `off.target.fluor` - normally an
#' off-target, unstained channel included to show spillover or its absence.
#'
#' @param setup.csv.path Path to the folder's setup CSV, as written by
#' `setup.unmix.comparison()`.
#' @param folder.path Path to the folder of already-unmixed FCS files, as
#' passed to `setup.unmix.comparison()`/`compare.unmix.folders()`.
#' @param on.target.fluor Character, the fluorophore whose single-stained
#' control file is read (plotted on the y-axis).
#' @param off.target.fluor Character, the fluorophore whose matched channel
#' is read from that same control file (plotted on the x-axis).
#' @param asp The AutoSpectral parameter list.
#' @param x.lab,y.lab Axis labels, passed to `create.biplot()`.
#' @param subtitle Character, added to the plot via `labs(subtitle = )`.
#' Typically the operator/method's display label.
#' @param subtitle.text.size Numeric, font size for `subtitle`. Default
#' `NULL` uses `asp$figure.axis.title.size`.
#' @param x.min,x.max,y.min,y.max Numeric, the untransformed axis bounds
#' passed to `create.biplot()`'s biexponential transform. Defaults match
#' `create.biplot()`'s own (`x.min`/`y.min = -5000`,
#' `x.max`/`y.max = asp$expr.data.max`).
#' @param x.width.basis,y.width.basis Width basis for the x/y-axis
#' biexponential transform, passed to `create.biplot()`. Default `-1000`
#' for both, matching the default used throughout
#' `plot_biplot_annotations.R`.
#'
#' @seealso
#' * [setup.unmix.comparison()]
#' * [match.fluorophores()]
#' * [build.operator.biplot.row()]
#'
#' @return A ggplot object.
#'
#' @export

build.operator.biplot <- function(
    setup.csv.path,
    folder.path,
    on.target.fluor,
    off.target.fluor,
    asp,
    x.lab,
    y.lab,
    subtitle = NULL,
    subtitle.text.size = NULL,
    x.min = -5000,
    x.max = asp$expr.data.max,
    y.min = -5000,
    y.max = asp$expr.data.max,
    x.width.basis = -1000,
    y.width.basis = -1000
) {

  if ( !file.exists( setup.csv.path ) )
    stop( paste0( "Setup CSV not found: '", setup.csv.path, "'." ), call. = FALSE )

  setup.table <- utils::read.csv(
    setup.csv.path, stringsAsFactors = FALSE, strip.white = TRUE
  )

  on.row <- setup.table[ setup.table$fluorophore == on.target.fluor, ]
  off.row <- setup.table[ setup.table$fluorophore == off.target.fluor, ]

  if ( nrow( on.row ) != 1 || is.na( on.row$filename[ 1 ] ) ||
       on.row$filename[ 1 ] == "No match" )
    stop(
      paste0(
        "No matched control file for '", on.target.fluor, "' in '",
        setup.csv.path, "'."
      ),
      call. = FALSE
    )

  if ( nrow( off.row ) != 1 || is.na( off.row$channel[ 1 ] ) ||
       off.row$channel[ 1 ] == "No match" )
    stop(
      paste0(
        "No matched channel for '", off.target.fluor, "' in '",
        setup.csv.path, "'."
      ),
      call. = FALSE
    )

  control.data <- readFCS( file.path( folder.path, on.row$filename[ 1 ] ) )

  biplot <- create.biplot(
    control.data,
    x.dim = off.row$channel[ 1 ],
    y.dim = on.row$channel[ 1 ],
    asp,
    x.lab = x.lab,
    y.lab = y.lab,
    x.min = x.min, x.max = x.max, y.min = y.min, y.max = y.max,
    x.width.basis = x.width.basis, y.width.basis = y.width.basis,
    save = FALSE
  )

  # Keeps every operator panel square regardless of how many columns Figure
  # 8 ends up giving this row, and stops subtitle text (added next) from
  # being clipped at the panel edge when the row renders small - the same
  # pattern used for the Figure 2 panels.
  biplot <- biplot + ggplot2::theme(
    aspect.ratio = 1,
    plot.margin = ggplot2::margin( 12, 12, 12, 12 )
  )
  biplot$coordinates$clip <- "off"

  if ( !is.null( subtitle ) ) {
    subtitle.size <- if ( is.null( subtitle.text.size ) )
      asp$figure.axis.title.size else subtitle.text.size
    biplot <- biplot +
      ggplot2::labs( subtitle = subtitle ) +
      ggplot2::theme( plot.subtitle = ggplot2::element_text( size = subtitle.size ) )
  }

  biplot
}


#' @title Build a Row of Operator Biplots With Shared Axis Labels
#'
#' @description
#' Calls `build.operator.biplot()` once per entry in `panel.order`, arranges
#' the results in a single row, and wraps the row with a larger shared axis
#' label (with an arrow) on each side via `add.biplot.shared.axis.labels()`.
#'
#' @param setup.dir Directory containing the `<label>_unmix_comparison_setup.csv`
#' files written by `setup.unmix.comparison()`.
#' @param folders Named character vector of folder paths, as passed to
#' `setup.unmix.comparison()`/`compare.unmix.folders()`. Names must include
#' every entry of `panel.order`.
#' @param panel.order Character vector, the folder/operator keys to include,
#' in display order. Must match both `names(folders)` and the setup CSV
#' filenames' `<label>` stem.
#' @param panel.labels Optional named character vector, keyed by entries of
#' `panel.order`, giving a display subtitle different from the key itself
#' (for example, relabeling a folder as "Operator 3 (beads)"). Entries not
#' present use the key as-is.
#' @param on.target.fluor,off.target.fluor As in `build.operator.biplot()`.
#' @param asp The AutoSpectral parameter list.
#' @param x.lab,y.lab Axis labels, used both on the small per-panel axes
#' (via `create.biplot()`) and on the larger shared axis-label strip.
#' @param subtitle.text.size Numeric, font size for each panel's operator
#' subtitle. Default `NULL` uses `asp$figure.axis.title.size`. Passed
#' through to `build.operator.biplot()`.
#' @param axis.label.text.size Numeric, font size for the larger shared
#' axis-label strip. Default `16`. Passed through to
#' `add.biplot.shared.axis.labels()`.
#' @param x.min,x.max,y.min,y.max Numeric, the untransformed axis bounds
#' passed to `create.biplot()`'s biexponential transform, applied
#' identically to every panel in the row. Passed through to
#' `build.operator.biplot()`.
#' @param x.width.basis,y.width.basis Width basis for the x/y-axis
#' biexponential transform, applied identically to every panel in the row.
#' Passed through to `build.operator.biplot()`.
#'
#' @seealso
#' * [build.operator.biplot()]
#' * [add.biplot.shared.axis.labels()]
#'
#' @return A cowplot/ggplot object: the row of biplots with shared axis
#' label strips added.
#'
#' @export

build.operator.biplot.row <- function(
    setup.dir,
    folders,
    panel.order,
    panel.labels = NULL,
    on.target.fluor,
    off.target.fluor,
    asp,
    x.lab,
    y.lab,
    subtitle.text.size = NULL,
    axis.label.text.size = 16,
    x.min = -5000,
    x.max = asp$expr.data.max,
    y.min = -5000,
    y.max = asp$expr.data.max,
    x.width.basis = -1000,
    y.width.basis = -1000
) {

  plots <- lapply( panel.order, function( key ) {

    subtitle <- if ( !is.null( panel.labels ) && key %in% names( panel.labels ) )
      panel.labels[[ key ]] else key

    build.operator.biplot(
      setup.csv.path = file.path(
        setup.dir, paste0( key, "_unmix_comparison_setup.csv" )
      ),
      folder.path = folders[[ key ]],
      on.target.fluor = on.target.fluor,
      off.target.fluor = off.target.fluor,
      asp = asp,
      x.lab = x.lab,
      y.lab = y.lab,
      subtitle = subtitle,
      subtitle.text.size = subtitle.text.size,
      x.min = x.min, x.max = x.max, y.min = y.min, y.max = y.max,
      x.width.basis = x.width.basis, y.width.basis = y.width.basis
    )
  } )

  row.plots <- cowplot::plot_grid( plotlist = plots, nrow = 1 )

  add.biplot.shared.axis.labels(
    row.plots, x.lab = x.lab, y.lab = y.lab, text.size = axis.label.text.size
  )
}
