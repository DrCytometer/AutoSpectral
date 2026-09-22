# plot_biplot_annotations.R

#' @title Position of the Biexponential Positive-Decade Parameter
#'
#' @description
#' Private helper computing the `pos` argument to `biexp.transform()` from a
#' given `maxValue` and `widthBasis`, following the same excess-width-basis
#' work-around used in `create.biplot()`. Shared by every function in this
#' file that needs to build a biexponential transform matching one produced
#' inside `create.biplot()`, so the logic is defined once.
#'
#' @param max.value Maximum value of the input (raw) scale.
#' @param width.basis Width basis for the biexponential transform.
#'
#' @return Numeric, the `pos` argument to pass to `biexp.transform()`.

.biexp.pos.log <- function( max.value, width.basis ) {

  if ( width.basis < -1000 ) {
    excess.width.basis <- width.basis + 1000
    pos.log.delta <- log10( abs( excess.width.basis ) )
    pos.log <- log10( max.value ) - 1 - pos.log.delta
    pmax( pos.log, 2 )
  } else {
    log10( max.value ) - 1
  }
}


#' @title Strip the Acquisition-Statistic Suffix From a Channel Name
#'
#' @description
#' Private helper used as the default display label wherever a metric
#' annotation (`annotate.biplot.metrics()`, `add.biplot.gated.metric()`)
#' would otherwise fall back to a raw FCS channel name (for example
#' "BUV496-A"). Strips a trailing "-A"/"-H"/"-W" acquisition-statistic
#' suffix, so a metric label reads "BUV496 SD: ..." rather than
#' "BUV496-A SD: ..." unless an explicit `channel.label`/`channel.labels`
#' is supplied.
#'
#' @param x Character vector of channel names.
#'
#' @return Character vector, same length as `x`, with any trailing
#' "-A"/"-H"/"-W" removed.

.strip.channel.suffix <- function( x ) {
  sub( "-[AHW]$", "", x )
}


#' @title Add a Gate Box to a Biplot
#'
#' @description
#' Overlays a single rectangular gate boundary onto a plot produced by
#' `create.biplot()`, using the same biexponential transform so the box lines
#' up with the transformed axes.
#'
#' @importFrom ggplot2 geom_path aes
#'
#' @param plot A ggplot object, typically returned by `create.biplot()` with
#' `save = FALSE`.
#' @param x.range Numeric vector of length 2, the untransformed lower and
#' upper x-axis bounds of the gate box.
#' @param y.range Numeric vector of length 2, the untransformed lower and
#' upper y-axis bounds of the gate box.
#' @param asp The AutoSpectral parameter list.
#' @param x.max Maximum value for the x-axis, used to build the biexponential
#' transform. Should match the value passed to `create.biplot()`. Default is
#' `asp$expr.data.max`.
#' @param y.max Maximum value for the y-axis, used to build the biexponential
#' transform. Should match the value passed to `create.biplot()`. Default is
#' `asp$expr.data.max`.
#' @param x.width.basis Width basis for the x-axis biexponential transform.
#' Should match the value passed to `create.biplot()`. Default is `-1000`.
#' @param y.width.basis Width basis for the y-axis biexponential transform.
#' Should match the value passed to `create.biplot()`. Default is `-1000`.
#' @param color Color of the gate box outline. Default is `black`.
#' @param linewidth Line width of the gate box outline. Default is `1`.
#'
#' @return The input `plot`, with the gate box added as an additional layer.
#'
#' @export

add.biplot.gate.box <- function(
    plot,
    x.range,
    y.range,
    asp,
    x.max = asp$expr.data.max,
    y.max = asp$expr.data.max,
    x.width.basis = -1000,
    y.width.basis = -1000,
    color = "black",
    linewidth = 1
) {

  biexp.trans.x <- biexp.transform(
    channelRange = asp$default.transformation.param$length,
    maxValue = x.max,
    pos = .biexp.pos.log( x.max, x.width.basis ),
    neg = asp$default.transformation.param$neg,
    widthBasis = x.width.basis,
    inverse = FALSE
  )

  biexp.trans.y <- biexp.transform(
    channelRange = asp$default.transformation.param$length,
    maxValue = y.max,
    pos = .biexp.pos.log( y.max, y.width.basis ),
    neg = asp$default.transformation.param$neg,
    widthBasis = y.width.basis,
    inverse = FALSE
  )

  box.data <- data.frame(
    x = biexp.trans.x( c(
      x.range[ 1 ], x.range[ 2 ], x.range[ 2 ], x.range[ 1 ], x.range[ 1 ]
    ) ),
    y = biexp.trans.y( c(
      y.range[ 1 ], y.range[ 1 ], y.range[ 2 ], y.range[ 2 ], y.range[ 1 ]
    ) )
  )

  plot + geom_path(
    data = box.data,
    aes( x = x, y = y ),
    color = color,
    linewidth = linewidth,
    inherit.aes = FALSE
  )
}


#' @title Add Multiple Gate Boxes to a Biplot
#'
#' @description
#' Overlays one or more rectangular gate boundaries onto a plot produced by
#' `create.biplot()`, by repeated calls to `add.biplot.gate.box()`. Useful for
#' representative plots showing, for example, separate positive and negative
#' gates on the same biplot.
#'
#' @param plot A ggplot object, typically returned by `create.biplot()` with
#' `save = FALSE`.
#' @param gates A list of gate specifications. Each element must be a list
#' with named entries `x.range` and `y.range` (each a numeric vector of
#' length 2, untransformed axis bounds), and may optionally include `color`
#' (a single color string for that gate's outline; defaults to `black` if
#' omitted).
#' @param asp The AutoSpectral parameter list.
#' @param x.max Maximum value for the x-axis, used to build the biexponential
#' transform. Should match the value passed to `create.biplot()`. Default is
#' `asp$expr.data.max`.
#' @param y.max Maximum value for the y-axis, used to build the biexponential
#' transform. Should match the value passed to `create.biplot()`. Default is
#' `asp$expr.data.max`.
#' @param x.width.basis Width basis for the x-axis biexponential transform.
#' Should match the value passed to `create.biplot()`. Default is `-1000`.
#' @param y.width.basis Width basis for the y-axis biexponential transform.
#' Should match the value passed to `create.biplot()`. Default is `-1000`.
#' @param linewidth Line width of the gate box outlines. Default is `1`.
#'
#' @seealso [add.biplot.gate.box()]
#'
#' @return The input `plot`, with each gate box added as an additional layer.
#'
#' @export

add.biplot.gate.boxes <- function(
    plot,
    gates,
    asp,
    x.max = asp$expr.data.max,
    y.max = asp$expr.data.max,
    x.width.basis = -1000,
    y.width.basis = -1000,
    linewidth = 1
) {

  for ( gate in gates ) {

    gate.color <- if ( is.null( gate$color ) ) "black" else gate$color

    plot <- add.biplot.gate.box(
      plot = plot,
      x.range = gate$x.range,
      y.range = gate$y.range,
      asp = asp,
      x.max = x.max,
      y.max = y.max,
      x.width.basis = x.width.basis,
      y.width.basis = y.width.basis,
      color = gate.color,
      linewidth = linewidth
    )
  }

  plot
}


#' @title Add Positivity Threshold Lines to a Biplot
#'
#' @description
#' Overlays a vertical and/or horizontal positivity threshold line onto a
#' plot produced by `create.biplot()`, using the same biexponential
#' transform so the lines line up with the transformed axes.
#'
#' @importFrom ggplot2 geom_vline geom_hline
#'
#' @param plot A ggplot object, typically returned by `create.biplot()` with
#' `save = FALSE`.
#' @param x.thresh Numeric, untransformed x-axis threshold. `NULL` (default)
#' omits the vertical line.
#' @param y.thresh Numeric, untransformed y-axis threshold. `NULL` (default)
#' omits the horizontal line.
#' @param asp The AutoSpectral parameter list.
#' @param x.max Maximum value for the x-axis, used to build the biexponential
#' transform. Should match the value passed to `create.biplot()`. Default is
#' `asp$expr.data.max`.
#' @param y.max Maximum value for the y-axis, used to build the biexponential
#' transform. Should match the value passed to `create.biplot()`. Default is
#' `asp$expr.data.max`.
#' @param x.width.basis Width basis for the x-axis biexponential transform.
#' Should match the value passed to `create.biplot()`. Default is `-1000`.
#' @param y.width.basis Width basis for the y-axis biexponential transform.
#' Should match the value passed to `create.biplot()`. Default is `-1000`.
#' @param x.color Color of the vertical threshold line. Default `"red"`.
#' @param y.color Color of the horizontal threshold line. Default `"blue"`.
#' @param linewidth Line width of both threshold lines. Default `1`.
#'
#' @return The input `plot`, with the threshold line(s) added as additional
#' layers.
#'
#' @export

add.biplot.threshold.lines <- function(
    plot,
    x.thresh = NULL,
    y.thresh = NULL,
    asp,
    x.max = asp$expr.data.max,
    y.max = asp$expr.data.max,
    x.width.basis = -1000,
    y.width.basis = -1000,
    x.color = "red",
    y.color = "blue",
    linewidth = 1
) {

  if ( !is.null( x.thresh ) ) {

    biexp.trans.x <- biexp.transform(
      channelRange = asp$default.transformation.param$length,
      maxValue = x.max,
      pos = .biexp.pos.log( x.max, x.width.basis ),
      neg = asp$default.transformation.param$neg,
      widthBasis = x.width.basis,
      inverse = FALSE
    )

    plot <- plot + geom_vline(
      xintercept = biexp.trans.x( x.thresh ),
      color = x.color,
      linewidth = linewidth
    )
  }

  if ( !is.null( y.thresh ) ) {

    biexp.trans.y <- biexp.transform(
      channelRange = asp$default.transformation.param$length,
      maxValue = y.max,
      pos = .biexp.pos.log( y.max, y.width.basis ),
      neg = asp$default.transformation.param$neg,
      widthBasis = y.width.basis,
      inverse = FALSE
    )

    plot <- plot + geom_hline(
      yintercept = biexp.trans.y( y.thresh ),
      color = y.color,
      linewidth = linewidth
    )
  }

  plot
}


#' @title Compute Quadrant Percentages for a Two-Channel Gate
#'
#' @description
#' Given a pair of positivity thresholds, computes the percentage of events
#' falling in each of the four quadrants they define. Used to annotate a
#' biplot with, for example, the percentage of events that are BUV661+
#' APC-, BUV661- APC+, and BUV661+ APC+.
#'
#' @param data A matrix or data frame with named columns, including `x.dim`
#' and `y.dim`.
#' @param x.dim String, the column of `data` the x threshold applies to.
#' @param y.dim String, the column of `data` the y threshold applies to.
#' @param x.thresh Numeric, the x-axis positivity threshold.
#' @param y.thresh Numeric, the y-axis positivity threshold.
#'
#' @return A named numeric vector of length 4, giving the percentage of
#' events in each quadrant: `x-y-`, `x+y-`, `x-y+`, `x+y+`.
#'
#' @export

compute.quadrant.percentages <- function(
    data,
    x.dim,
    y.dim,
    x.thresh,
    y.thresh
) {

  x.pos <- data[ , x.dim ] > x.thresh
  y.pos <- data[ , y.dim ] > y.thresh

  n <- nrow( data )

  c(
    "x-y-" = 100 * sum( !x.pos & !y.pos ) / n,
    "x+y-" = 100 * sum(  x.pos & !y.pos ) / n,
    "x-y+" = 100 * sum( !x.pos &  y.pos ) / n,
    "x+y+" = 100 * sum(  x.pos &  y.pos ) / n
  )
}


#' @title Annotate a Biplot with Quadrant Percentages
#'
#' @description
#' Adds the percentage of events falling in each of the four quadrants
#' defined by an x and y positivity threshold (see
#' `compute.quadrant.percentages()`) to the corresponding corner of a plot
#' produced by `create.biplot()`. Threshold lines are added first, via
#' `add.biplot.threshold.lines()`.
#'
#' @importFrom ggplot2 annotate
#'
#' @param plot A ggplot object, typically returned by `create.biplot()` with
#' `save = FALSE`.
#' @param data The data plotted, used to compute the quadrant percentages.
#' Should match the data passed to `create.biplot()` for this plot.
#' @param x.dim String, the column of `data` for the x-axis. Should match
#' the `x.dim` passed to `create.biplot()`.
#' @param y.dim String, the column of `data` for the y-axis. Should match
#' the `y.dim` passed to `create.biplot()`.
#' @param x.thresh Numeric, the x-axis positivity threshold, untransformed.
#' @param y.thresh Numeric, the y-axis positivity threshold, untransformed.
#' @param asp The AutoSpectral parameter list.
#' @param x.max Maximum value for the x-axis biexponential transform.
#' Default is `asp$expr.data.max`.
#' @param y.max Maximum value for the y-axis biexponential transform.
#' Default is `asp$expr.data.max`.
#' @param x.width.basis Width basis for the x-axis biexponential transform.
#' Default is `-1000`.
#' @param y.width.basis Width basis for the y-axis biexponential transform.
#' Default is `-1000`.
#' @param x.color Color of the vertical threshold line. Default `"red"`.
#' @param y.color Color of the horizontal threshold line. Default `"blue"`.
#' @param linewidth Line width of the threshold lines. Default `1`.
#' @param digits Integer, decimal places for the percentage labels. Default
#' `1`.
#' @param text.size Numeric, font size for the percentage labels. Default
#' `NULL` derives it from `asp$figure.axis.text.size`.
#'
#' @seealso
#' * [compute.quadrant.percentages()]
#' * [add.biplot.threshold.lines()]
#'
#' @return The input `plot`, with threshold lines and quadrant percentage
#' labels added.
#'
#' @export

annotate.biplot.quadrants <- function(
    plot,
    data,
    x.dim,
    y.dim,
    x.thresh,
    y.thresh,
    asp,
    x.max = asp$expr.data.max,
    y.max = asp$expr.data.max,
    x.width.basis = -1000,
    y.width.basis = -1000,
    x.color = "red",
    y.color = "blue",
    linewidth = 1,
    digits = 1,
    text.size = NULL
) {

  plot <- add.biplot.threshold.lines(
    plot = plot,
    x.thresh = x.thresh,
    y.thresh = y.thresh,
    asp = asp,
    x.max = x.max,
    y.max = y.max,
    x.width.basis = x.width.basis,
    y.width.basis = y.width.basis,
    x.color = x.color,
    y.color = y.color,
    linewidth = linewidth
  )

  quadrants <- compute.quadrant.percentages( data, x.dim, y.dim, x.thresh, y.thresh )

  if ( is.null( text.size ) ) text.size <- asp$figure.axis.text.size / 2.8

  fmt <- function( value ) sprintf( paste0( "%.", digits, "f%%" ), value )

  plot +
    annotate(
      "text", x = -Inf, y = -Inf, hjust = -0.1, vjust = -0.5,
      label = fmt( quadrants[ "x-y-" ] ), size = text.size
    ) +
    annotate(
      "text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
      label = fmt( quadrants[ "x+y-" ] ), size = text.size
    ) +
    annotate(
      "text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
      label = fmt( quadrants[ "x-y+" ] ), size = text.size
    ) +
    annotate(
      "text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
      label = fmt( quadrants[ "x+y+" ] ), size = text.size
    )
}


#' @title Add Power-of-Ten Gridlines to a Biplot
#'
#' @description
#' Adds light gray reference gridlines at each power-of-ten break on both
#' axes of a plot produced by `create.biplot()`, using the same
#' biexponential transform, and by default the same break values,
#' `create.biplot()` itself uses for its axis labels (`asp$ribbon.breaks`).
#' `create.biplot()`'s own theme blanks `panel.grid.major`/`panel.grid.minor`
#' (`element_blank()`), so gridlines are not drawn by default; this adds
#' them back as an explicit, optional layer. The gridlines are inserted as
#' the first layer of the plot, so they render behind the density and point
#' layers rather than obscuring them.
#'
#' @importFrom ggplot2 geom_vline geom_hline
#'
#' @param plot A ggplot object, typically returned by `create.biplot()` with
#' `save = FALSE`.
#' @param asp The AutoSpectral parameter list.
#' @param x.max Maximum value for the x-axis, used to build the biexponential
#' transform. Should match the value passed to `create.biplot()`. Default is
#' `asp$expr.data.max`.
#' @param y.max Maximum value for the y-axis, used to build the biexponential
#' transform. Should match the value passed to `create.biplot()`. Default is
#' `asp$expr.data.max`.
#' @param x.width.basis Width basis for the x-axis biexponential transform.
#' Should match the value passed to `create.biplot()`. Default is `-1000`.
#' @param y.width.basis Width basis for the y-axis biexponential transform.
#' Should match the value passed to `create.biplot()`. Default is `-1000`.
#' @param breaks Numeric vector of untransformed axis values to draw
#' gridlines at. Default `NULL` uses `asp$ribbon.breaks`, the same set
#' `create.biplot()` labels its axes with.
#' @param color Gridline color. Default `"gray80"`.
#' @param linewidth Gridline width. Default `0.3`.
#'
#' @return The input `plot`, with gridlines inserted as its first layer.
#'
#' @export

add.biplot.gridlines <- function(
    plot,
    asp,
    x.max = asp$expr.data.max,
    y.max = asp$expr.data.max,
    x.width.basis = -1000,
    y.width.basis = -1000,
    breaks = NULL,
    color = "gray80",
    linewidth = 0.3
) {

  if ( is.null( breaks ) ) breaks <- asp$ribbon.breaks

  x.breaks <- breaks[ breaks < x.max ]
  y.breaks <- breaks[ breaks < y.max ]

  biexp.trans.x <- biexp.transform(
    channelRange = asp$default.transformation.param$length,
    maxValue = x.max,
    pos = .biexp.pos.log( x.max, x.width.basis ),
    neg = asp$default.transformation.param$neg,
    widthBasis = x.width.basis,
    inverse = FALSE
  )

  biexp.trans.y <- biexp.transform(
    channelRange = asp$default.transformation.param$length,
    maxValue = y.max,
    pos = .biexp.pos.log( y.max, y.width.basis ),
    neg = asp$default.transformation.param$neg,
    widthBasis = y.width.basis,
    inverse = FALSE
  )

  grid.layers <- list(
    geom_vline(
      xintercept = biexp.trans.x( x.breaks ),
      color = color,
      linewidth = linewidth
    ),
    geom_hline(
      yintercept = biexp.trans.y( y.breaks ),
      color = color,
      linewidth = linewidth
    )
  )

  plot$layers <- c( grid.layers, plot$layers )

  plot
}


#' @title Compute a Summary Metric Per Channel
#'
#' @description
#' Computes one summary statistic per named column of `data`. Three of the
#' four available methods follow conventions already used elsewhere in
#' AutoSpectral: `rsd` is the robust standard deviation (median absolute
#' deviation, as used by `calculate.ssi()` and `assess.variability()`), and
#' `mfi` is the median (as used for positive/negative comparisons in
#' `compare_unmix_folders.R`). `sd` is the plain (non-robust) standard
#' deviation, kept distinct from `rsd` since the two are not interchangeable.
#'
#' @param data A matrix or data frame with named columns.
#' @param channels Character vector of column names in `data` to summarize.
#' @param method One of `rsd` (median absolute deviation), `mfi` (median),
#' `sd` (standard deviation), or `mean`.
#'
#' @return A named numeric vector, one value per element of `channels`.
#'
#' @export

compute.channel.metric <- function(
    data,
    channels,
    method = c( "rsd", "mfi", "sd", "mean" )
) {

  method <- match.arg( method )

  metric.fun <- switch(
    method,
    rsd  = stats::mad,
    mfi  = stats::median,
    sd   = stats::sd,
    mean = base::mean
  )

  values <- vapply( channels, function( ch ) metric.fun( data[ , ch ] ), numeric( 1 ) )
  names( values ) <- channels

  values
}


#' @title Compare a Channel Metric Between Two Datasets
#'
#' @description
#' Computes a summary metric (see `compute.channel.metric()`) for the same
#' set of channels in two datasets, and returns the paired values in one data
#' frame - for example, comparing the robust SD of each fluorophore channel
#' before and after a change in unmixing or autofluorescence extraction.
#' `data.1` and `data.2` need not share a channel-naming convention: when
#' they don't (for example, two independently-unmixed exports), resolve
#' each dataset's own column name for the shared set of fluorophores with
#' `resolve.fluorophore.channels()` first, and pass the two resulting
#' column-name vectors as `channels` and `channels.2` - `channels[i]` and
#' `channels.2[i]` are then treated as the same fluorophore's channel in
#' `data.1` and `data.2` respectively, whatever each is individually named.
#'
#' @param data.1 A matrix or data frame with named columns.
#' @param data.2 A second matrix or data frame with named columns. Shares
#' `channels` with `data.1` when `channels.2` is not supplied.
#' @param channels Character vector of column names in `data.1` to compare.
#' @param channels.2 Character vector of column names in `data.2`, the same
#' length as `channels` and in the same order, giving `data.2`'s own column
#' name for each entry of `channels`. Default `NULL` uses `channels` for
#' both datasets (they must then share a naming convention).
#' @param method One of `rsd`, `mfi`, `sd`, or `mean`. See
#' `compute.channel.metric()`.
#' @param label.1 Column name to use for `data.1`'s values in the result.
#' Default is `x`.
#' @param label.2 Column name to use for `data.2`'s values in the result.
#' Default is `y`.
#' @param channel.labels Character vector, same length as `channels`, used
#' as the result's `Channel` column in place of `channels` itself - useful
#' when `channels`/`channels.2` are two datasets' differing raw column
#' names for the same fluorophores. Default `NULL` uses `channels`.
#'
#' @seealso
#' * [compute.channel.metric()]
#' * [resolve.fluorophore.channels()]
#'
#' @return A data frame with one row per channel and columns `Channel`,
#' `label.1`, and `label.2`.
#'
#' @export

compare.channel.metric <- function(
    data.1,
    data.2,
    channels,
    channels.2 = NULL,
    method = c( "rsd", "mfi", "sd", "mean" ),
    label.1 = "x",
    label.2 = "y",
    channel.labels = NULL
) {

  method <- match.arg( method )

  if ( is.null( channels.2 ) ) channels.2 <- channels
  if ( length( channels.2 ) != length( channels ) )
    stop( "`channels.2` must be the same length as `channels`.", call. = FALSE )

  if ( is.null( channel.labels ) ) channel.labels <- channels

  values.1 <- compute.channel.metric( data.1, channels, method )
  values.2 <- compute.channel.metric( data.2, channels.2, method )

  result <- data.frame( Channel = channel.labels, stringsAsFactors = FALSE )
  result[[ label.1 ]] <- as.numeric( values.1 )
  result[[ label.2 ]] <- as.numeric( values.2 )

  result
}


#' @title Annotate a Biplot with Whole-Plot Channel Metrics
#'
#' @description
#' Adds a text block to a plot produced by `create.biplot()`, reporting a
#' summary metric (see `compute.channel.metric()`) for one or more channels,
#' computed over every event shown on the plot (as opposed to
#' `add.biplot.gated.metric()`, which restricts the calculation to events
#' inside a gate box).
#'
#' @importFrom ggplot2 annotate
#'
#' @param plot A ggplot object, typically returned by `create.biplot()` with
#' `save = FALSE`.
#' @param data The data plotted, used to compute the metric. Should match
#' the data passed to `create.biplot()` for this plot.
#' @param channels Character vector of column names in `data` to report.
#' @param asp The AutoSpectral parameter list.
#' @param method One of `rsd`, `mfi`, `sd`, or `mean`. See
#' `compute.channel.metric()`.
#' @param metric.label Character, the label used for the metric in each
#' printed line (for example `rSD` or `SD`). Default `NULL` derives it from
#' `method`.
#' @param channel.labels Character vector, same length as `channels`, used
#' in place of the raw channel names in each printed line. Default `NULL`
#' strips any trailing "-A"/"-H"/"-W" acquisition-statistic suffix from
#' `channels` (see `.strip.channel.suffix()`) - for example "BUV496-A"
#' becomes "BUV496".
#' @param position One of `topleft`, `topright`, `bottomleft`, `bottomright`,
#' giving the corner of the plot the text block is anchored to. Default is
#' `topleft`.
#' @param text.size Numeric, font size for the annotation. Default `NULL`
#' derives it from `asp$figure.axis.text.size`.
#'
#' @return The input `plot`, with the metric annotation added.
#'
#' @export

annotate.biplot.metrics <- function(
    plot,
    data,
    channels,
    asp,
    method = c( "rsd", "mfi", "sd", "mean" ),
    metric.label = NULL,
    channel.labels = NULL,
    position = c( "topleft", "topright", "bottomleft", "bottomright" ),
    text.size = NULL
) {

  method <- match.arg( method )
  position <- match.arg( position )

  if ( is.null( metric.label ) )
    metric.label <- switch( method, rsd = "rSD", mfi = "MFI", sd = "SD", mean = "Mean" )

  if ( is.null( channel.labels ) ) channel.labels <- .strip.channel.suffix( channels )

  values <- compute.channel.metric( data, channels, method )

  label.text <- paste(
    sprintf( "%s %s: %s", channel.labels, metric.label, signif( values, 3 ) ),
    collapse = "\n"
  )

  is.left <- grepl( "left", position )
  is.top <- grepl( "top", position )

  if ( is.null( text.size ) ) text.size <- asp$figure.axis.text.size / 2.8

  plot + annotate(
    "text",
    x = if ( is.left ) -Inf else Inf,
    y = if ( is.top ) Inf else -Inf,
    hjust = if ( is.left ) -0.05 else 1.05,
    vjust = if ( is.top ) 1.1 else -0.1,
    label = label.text,
    size = text.size,
    lineheight = 0.9
  )
}


#' @title Add a Gate Box with a Within-Gate Metric Annotation
#'
#' @description
#' Draws a single rectangular gate box on a plot produced by
#' `create.biplot()` (via `add.biplot.gate.box()`), computes a summary
#' metric (see `compute.channel.metric()`) for one channel restricted to the
#' events falling inside that box, and labels the box with the result.
#'
#' @importFrom ggplot2 annotate
#'
#' @param plot A ggplot object, typically returned by `create.biplot()` with
#' `save = FALSE`.
#' @param data The data plotted, used both to identify events inside the
#' gate box and to compute the metric. Should match the data passed to
#' `create.biplot()` for this plot.
#' @param x.dim String, the column of `data` the gate box's x-axis bounds
#' apply to. Should match the `x.dim` passed to `create.biplot()`.
#' @param y.dim String, the column of `data` the gate box's y-axis bounds
#' apply to. Should match the `y.dim` passed to `create.biplot()`.
#' @param x.range Numeric vector of length 2, the untransformed lower and
#' upper x-axis bounds of the gate box.
#' @param y.range Numeric vector of length 2, the untransformed lower and
#' upper y-axis bounds of the gate box.
#' @param metric.channel String, the column of `data` the metric is computed
#' on, for events inside the gate box. Often, but not necessarily, the same
#' as `x.dim`.
#' @param asp The AutoSpectral parameter list.
#' @param method One of `rsd`, `mfi`, `sd`, or `mean`. See
#' `compute.channel.metric()`.
#' @param metric.label Character, the label used for the metric in the
#' printed line (for example `rSD` or `MFI`). Default `NULL` derives it from
#' `method`.
#' @param channel.label Character, used in place of the raw channel name in
#' the printed line. Default `NULL` strips any trailing "-A"/"-H"/"-W"
#' acquisition-statistic suffix from `metric.channel` (see
#' `.strip.channel.suffix()`) - for example "BV650-A" becomes "BV650".
#' @param x.max Maximum value for the x-axis, used to build the biexponential
#' transform. Should match the value passed to `create.biplot()`. Default is
#' `asp$expr.data.max`.
#' @param y.max Maximum value for the y-axis, used to build the biexponential
#' transform. Should match the value passed to `create.biplot()`. Default is
#' `asp$expr.data.max`.
#' @param x.width.basis Width basis for the x-axis biexponential transform.
#' Should match the value passed to `create.biplot()`. Default is `-1000`.
#' @param y.width.basis Width basis for the y-axis biexponential transform.
#' Should match the value passed to `create.biplot()`. Default is `-1000`.
#' @param box.color Color of the gate box outline and its metric label.
#' Default is `black`.
#' @param linewidth Line width of the gate box outline. Default is `1`.
#' @param text.size Numeric, font size for the metric label. Default `NULL`
#' derives it from `asp$figure.axis.text.size`.
#'
#' @seealso [add.biplot.gate.box()]
#'
#' @return The input `plot`, with the gate box and its metric label added.
#'
#' @export

add.biplot.gated.metric <- function(
    plot,
    data,
    x.dim,
    y.dim,
    x.range,
    y.range,
    metric.channel,
    asp,
    method = c( "rsd", "mfi", "sd", "mean" ),
    metric.label = NULL,
    channel.label = NULL,
    x.max = asp$expr.data.max,
    y.max = asp$expr.data.max,
    x.width.basis = -1000,
    y.width.basis = -1000,
    box.color = "black",
    linewidth = 1,
    text.size = NULL
) {

  method <- match.arg( method )

  if ( is.null( metric.label ) )
    metric.label <- switch( method, rsd = "rSD", mfi = "MFI", sd = "SD", mean = "Mean" )

  if ( is.null( channel.label ) ) channel.label <- .strip.channel.suffix( metric.channel )

  gated.idx <- which(
    data[ , x.dim ] >= min( x.range ) & data[ , x.dim ] <= max( x.range ) &
      data[ , y.dim ] >= min( y.range ) & data[ , y.dim ] <= max( y.range )
  )

  if ( length( gated.idx ) == 0 ) {
    warning(
      paste0(
        "No events fall within the gate box (x.range = ", paste( x.range, collapse = ", " ),
        ", y.range = ", paste( y.range, collapse = ", " ), ") on ", x.dim, " vs ", y.dim, "."
      ),
      call. = FALSE
    )
    value <- NA_real_
  } else {
    value <- compute.channel.metric( data[ gated.idx, , drop = FALSE ], metric.channel, method )
  }

  plot <- add.biplot.gate.box(
    plot = plot,
    x.range = x.range,
    y.range = y.range,
    asp = asp,
    x.max = x.max,
    y.max = y.max,
    x.width.basis = x.width.basis,
    y.width.basis = y.width.basis,
    color = box.color,
    linewidth = linewidth
  )

  biexp.trans.x <- biexp.transform(
    channelRange = asp$default.transformation.param$length,
    maxValue = x.max,
    pos = .biexp.pos.log( x.max, x.width.basis ),
    neg = asp$default.transformation.param$neg,
    widthBasis = x.width.basis,
    inverse = FALSE
  )

  biexp.trans.y <- biexp.transform(
    channelRange = asp$default.transformation.param$length,
    maxValue = y.max,
    pos = .biexp.pos.log( y.max, y.width.basis ),
    neg = asp$default.transformation.param$neg,
    widthBasis = y.width.basis,
    inverse = FALSE
  )

  label.x <- biexp.trans.x( max( x.range ) )
  label.y <- biexp.trans.y( max( y.range ) )

  if ( is.null( text.size ) ) text.size <- asp$figure.axis.text.size / 2.8

  label.text <- sprintf( "%s %s: %s", channel.label, metric.label, signif( value, 3 ) )

  plot + annotate(
    "text",
    x = label.x, y = label.y,
    hjust = -0.05, vjust = -0.3,
    label = label.text,
    color = box.color,
    size = text.size
  )
}


#' @title Linear-Scale Density Biplot
#'
#' @description
#' Private helper producing a pseudocolour-density biplot on linear axes,
#' following the same visual conventions as `create.biplot()` (scattermore
#' point layer, 2D density contour fill, package figure theme), but without
#' the biexponential transform. Used as the base layer for
#' `regression.biplot()`.
#'
#' @importFrom ggplot2 ggplot aes after_stat stat_density_2d
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous theme_bw theme
#' @importFrom ggplot2 margin element_line element_text element_rect element_blank
#' @importFrom ggplot2 scale_fill_viridis_c scale_fill_gradientn
#' @importFrom scattermore geom_scattermore
#'
#' @param data A matrix or data frame containing the flow cytometry data to
#' be plotted.
#' @param x.dim String specifying the column of `data` for the x-axis.
#' @param y.dim String specifying the column of `data` for the y-axis.
#' @param asp The AutoSpectral parameter list.
#' @param x.lab Label for the x-axis. Defaults to `x.dim`.
#' @param y.lab Label for the y-axis. Defaults to `y.dim`.
#' @param color.palette Viridis palette name, or `rainbow` for the package's
#' default gradient. Default is `rainbow`.
#' @param max.points Number of points to plot. Default is `5e5`.
#'
#' @return A ggplot object.

.linear.density.biplot <- function(
    data,
    x.dim,
    y.dim,
    asp,
    x.lab = NULL,
    y.lab = NULL,
    color.palette = "rainbow",
    max.points = 5e5
) {

  if ( is.null( x.lab ) ) x.lab <- x.dim
  if ( is.null( y.lab ) ) y.lab <- y.dim

  if ( nrow( data ) > max.points ) {
    set.seed( asp$bird.seed )
    data <- data[ sample( seq_len( nrow( data ) ), max.points ), ]
  }

  plot.data <- data.frame(
    x = data[ , x.dim ],
    y = data[ , y.dim ]
  )

  biplot <- ggplot( plot.data, aes( x, y ) ) +
    geom_scattermore(
      pointsize = asp$figure.gate.point.size,
      color = "black",
      alpha = 1,
      na.rm = TRUE
    ) +
    stat_density_2d(
      aes( fill = after_stat( level ) ),
      geom = "polygon",
      na.rm = TRUE
    ) +
    scale_x_continuous( name = x.lab ) +
    scale_y_continuous( name = y.lab ) +
    theme_bw() +
    theme(
      plot.margin = margin(
        asp$figure.margin, asp$figure.margin, asp$figure.margin, asp$figure.margin
      ),
      legend.position = "none",
      axis.ticks = element_line( linewidth = asp$figure.panel.line.size ),
      axis.text = element_text( size = asp$figure.axis.text.size ),
      axis.title = element_text( size = asp$figure.axis.title.size ),
      panel.border = element_rect( fill = NA, linewidth = asp$figure.panel.line.size ),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )

  viridis.colors <- c(
    "magma", "inferno", "plasma", "viridis",
    "cividis", "rocket", "mako", "turbo"
  )

  if ( color.palette %in% viridis.colors ) {
    biplot <- biplot + scale_fill_viridis_c( option = color.palette )
  } else {
    biplot <- biplot + scale_fill_gradientn( colors = asp$density.palette.base.color )
  }

  biplot
}


#' @title Regression Biplot
#'
#' @description
#' Produces a linear-scale density biplot with a fitted regression line
#' overlaid, using either ordinary least squares (`lm`) or a robust linear
#' model (`rlm`). The proportion of variance explained (R-squared) is
#' annotated on the plot; no significance test is reported, since a p-value
#' from either fit is not informative about how well the line summarizes the
#' bulk of biexponentially-distributed spectral flow data.
#'
#' @importFrom ggplot2 geom_abline annotate
#' @importFrom ggplot2 ggsave
#' @importFrom ragg agg_jpeg
#' @importFrom stats lm coef fitted cor
#' @importFrom MASS rlm
#'
#' @param data A matrix or data frame containing the flow cytometry data to
#' be plotted, on a linear scale.
#' @param x.dim String specifying the column of `data` for the x-axis
#' (predictor).
#' @param y.dim String specifying the column of `data` for the y-axis
#' (response).
#' @param asp The AutoSpectral parameter list.
#' @param method Either `lm` (ordinary least squares) or `rlm` (robust linear
#' model, `MASS::rlm`).
#' @param line.color Color for the fitted line. Default is `NULL`, which
#' uses `red` for `lm` and `blue` for `rlm`.
#' @param x.lab Label for the x-axis. Defaults to `x.dim`.
#' @param y.lab Label for the y-axis. Defaults to `y.dim`.
#' @param rlm.maxit Maximum number of iterations for `MASS::rlm`. Default is
#' `100`.
#' @param save Logical, if `TRUE`, saves a JPEG file to `output.dir`.
#' @param title Optional title for the plot filename. Defaults to
#' `x.lab vs y.lab method`.
#' @param output.dir Optional output directory. Default is the current
#' working directory.
#' @param width Numeric, width of the saved plot. Default is `5`.
#' @param height Numeric, height of the saved plot. Default is `5`.
#'
#' @return A ggplot object.
#'
#' @export

regression.biplot <- function(
    data,
    x.dim,
    y.dim,
    asp,
    method = c( "lm", "rlm" ),
    line.color = NULL,
    x.lab = NULL,
    y.lab = NULL,
    rlm.maxit = 100,
    save = TRUE,
    title = NULL,
    output.dir = NULL,
    width = 5,
    height = 5
) {

  method <- match.arg( method )

  if ( is.null( x.lab ) ) x.lab <- x.dim
  if ( is.null( y.lab ) ) y.lab <- y.dim

  if ( is.null( line.color ) )
    line.color <- if ( method == "lm" ) "red" else "blue"

  x.vals <- data[ , x.dim ]
  y.vals <- data[ , y.dim ]

  if ( method == "lm" ) {

    fit <- lm( y.vals ~ x.vals )
    r.squared <- summary( fit )$r.squared
    intercept <- coef( fit )[ 1 ]
    slope <- coef( fit )[ 2 ]

  } else {

    fit <- rlm( y.vals ~ x.vals, maxit = rlm.maxit )
    intercept <- coef( fit )[ 1 ]
    slope <- coef( fit )[ 2 ]
    r.squared <- cor( fitted( fit ), y.vals ) ^ 2

  }

  stats.label <- sprintf( "R^2 == %.3f", r.squared )

  biplot <- .linear.density.biplot(
    data = data,
    x.dim = x.dim,
    y.dim = y.dim,
    asp = asp,
    x.lab = x.lab,
    y.lab = y.lab
  )

  biplot <- biplot +
    geom_abline(
      intercept = intercept,
      slope = slope,
      color = line.color,
      linewidth = 1
    ) +
    annotate(
      "text",
      x = -Inf, y = Inf,
      hjust = -0.1, vjust = 1.5,
      label = stats.label,
      parse = TRUE,
      size = asp$figure.axis.text.size / 2.8
    )

  if ( is.null( title ) )
    title <- paste( x.lab, "vs", y.lab, method )

  if ( is.null( output.dir ) )
    output.dir <- getwd()

  if ( save ) {
    ggsave(
      file.path( output.dir, sprintf( "%s.jpg", title ) ),
      plot = biplot,
      device = ragg::agg_jpeg,
      width = width,
      height = height,
      limitsize = FALSE
    )
  }

  biplot
}


#' @title Compute a Scatter Gate Boundary
#'
#' @description
#' Computes a scatter gate boundary on a flow cytometry sample's own scatter
#' data, without applying it. This separates gate computation from gate
#' application (`apply.gate()`) so that a single boundary, computed once from
#' one representative dataset, can be applied identically to several other
#' datasets (for example, the same sample unmixed by two different
#' pipelines), rather than each being gated independently.
#'
#' @param flow.data A matrix or data frame of flow cytometry data with named
#' columns, including the two scatter parameters named in `scatter.param`.
#' @param asp The AutoSpectral parameter list.
#' @param scatter.param Character vector of length 2 giving the names of the
#' two scatter columns in `flow.data` to gate on. Default is
#' `asp$default.scatter.parameter`.
#' @param large.gate Logical, whether to extend the gate upwards and
#' outwards for larger cells. Default is `FALSE`.
#' @param viability.gate Logical, whether to extend the gate to the left to
#' include dead cells. Default is `FALSE`.
#' @param control.type Either `cells` or `beads`, selecting which set of
#' `asp` gating defaults to use. Default is `cells`.
#' @param samp Character, a label for this sample, used only for the
#' diagnostic gate plot `do.gate()` produces.
#' @param output.dir File path where the diagnostic gate plot is saved.
#' Default is `./figure_gate`.
#' @param color.palette Viridis palette used for the diagnostic gate plot.
#' Default is `plasma`.
#' @param max.points Number of points to show on the diagnostic gate plot.
#' Default is `5e4`.
#' @param gate.color Color of the gate boundary line on the diagnostic gate
#' plot. Default is `darkgoldenrod1`.
#'
#' @seealso
#' * [do.gate()]
#' * [apply.gate()]
#' * [gate.large.sample()]
#'
#' @return A gate boundary, as returned by `do.gate()`.
#'
#' @export

compute.scatter.gate <- function(
    flow.data,
    asp,
    scatter.param = asp$default.scatter.parameter,
    large.gate = FALSE,
    viability.gate = FALSE,
    control.type = "cells",
    samp = "sample",
    output.dir = "./figure_gate",
    color.palette = "plasma",
    max.points = 5e4,
    gate.color = "darkgoldenrod1"
) {

  asp$figure.gate.dir <- output.dir
  if ( !dir.exists( asp$figure.gate.dir ) )
    dir.create( asp$figure.gate.dir, recursive = TRUE )

  gate.data <- flow.data[ , scatter.param ]

  do.gate(
    gate.data = gate.data,
    viability.gate = viability.gate,
    large.gate = large.gate,
    samp = samp,
    scatter.and.channel.label = scatter.param,
    control.type = control.type,
    asp = asp,
    color.palette = color.palette,
    max.points = max.points,
    gate.color = gate.color
  )
}


#' @title Pre-Gate a Sample Using a Large Scatter Gate
#'
#' @description
#' Computes a large scatter gate on a flow cytometry sample and returns only
#' the events falling inside it. This is a standalone counterpart to the
#' `large.gate` option used within the control-file pipeline
#' (`define.gate.density()`, `define.gate.landmarks()`), for use on an
#' arbitrary sample (for example, a fully-stained sample being unmixed for a
#' representative plot) that is not itself part of a control file.
#'
#' Internally, the gate boundary is computed by `compute.scatter.gate()` and
#' applied with `apply.gate()`. To apply the same gate boundary to several
#' datasets (for example, a sample unmixed by two different pipelines), call
#' those two functions directly rather than this convenience wrapper.
#'
#' @param flow.data A matrix or data frame of flow cytometry data (raw or
#' unmixed) with named columns, including the two scatter parameters named
#' in `scatter.param`.
#' @param asp The AutoSpectral parameter list.
#' @param scatter.param Character vector of length 2 giving the names of the
#' two scatter columns in `flow.data` to gate on. Default is
#' `asp$default.scatter.parameter`.
#' @param large.gate Logical, whether to extend the gate upwards and
#' outwards for larger cells. Default is `TRUE`.
#' @param viability.gate Logical, whether to extend the gate to the left to
#' include dead cells. Default is `FALSE`.
#' @param control.type Either `cells` or `beads`, selecting which set of
#' `asp` gating defaults to use. Default is `cells`.
#' @param samp Character, a label for this sample, used only for the
#' diagnostic gate plot `do.gate()` produces.
#' @param output.dir File path where the diagnostic gate plot is saved.
#' Default is `./figure_gate`.
#' @param color.palette Viridis palette used for the diagnostic gate plot.
#' Default is `plasma`.
#' @param max.points Number of points to show on the diagnostic gate plot.
#' Default is `5e4`.
#' @param gate.color Color of the gate boundary line on the diagnostic gate
#' plot. Default is `darkgoldenrod1`.
#' @param min.fraction Numeric between `0` and `1`, passed to `apply.gate()`.
#' Default is `0.01`.
#'
#' @seealso
#' * [compute.scatter.gate()]
#' * [apply.gate()]
#'
#' @return `flow.data`, subset to only those events falling inside the
#' computed large gate.
#'
#' @export

gate.large.sample <- function(
    flow.data,
    asp,
    scatter.param = asp$default.scatter.parameter,
    large.gate = TRUE,
    viability.gate = FALSE,
    control.type = "cells",
    samp = "sample",
    output.dir = "./figure_gate",
    color.palette = "plasma",
    max.points = 5e4,
    gate.color = "darkgoldenrod1",
    min.fraction = 0.01
) {

  gate.boundary <- compute.scatter.gate(
    flow.data = flow.data,
    asp = asp,
    scatter.param = scatter.param,
    large.gate = large.gate,
    viability.gate = viability.gate,
    control.type = control.type,
    samp = samp,
    output.dir = output.dir,
    color.palette = color.palette,
    max.points = max.points,
    gate.color = gate.color
  )

  apply.gate(
    flow.data = flow.data,
    gate.boundary = gate.boundary,
    scatter.param = scatter.param,
    asp = asp,
    min.fraction = min.fraction
  )
}


#' @title Build a Shared Axis Label Strip
#'
#' @description
#' Builds a small, borderless panel containing a single axis label with a
#' single-headed arrow spanning it - arrowhead at the higher-value end (the
#' right end for a horizontal, x-axis strip; the top end for a vertical,
#' y-axis strip) - sized to sit alongside a row or column of biplots as a
#' larger shared axis label (see `add.biplot.shared.axis.labels()`), in
#' place of relying on each individual biplot's own, smaller axis title.
#'
#' @importFrom ggplot2 ggplot coord_cartesian theme_void annotate
#' @importFrom grid arrow unit
#'
#' @param label Character, the axis label text.
#' @param orientation One of `horizontal` (a shared x-axis label, meant to
#' sit below a row of biplots) or `vertical` (a shared y-axis label, meant
#' to sit to the left of a row of biplots; the label text is rotated 90
#' degrees).
#' @param text.size Numeric, font size for the label. Default `16`.
#' @param arrow.length Numeric between `0` and `1`, the fraction of the
#' strip's length spanned by the arrow. Default `0.9`.
#'
#' @seealso [add.biplot.shared.axis.labels()]
#'
#' @return A ggplot object: a borderless panel containing the label and its
#' arrow.
#'
#' @export

build.axis.label.strip <- function(
    label,
    orientation = c( "horizontal", "vertical" ),
    text.size = 16,
    arrow.length = 0.9
) {

  orientation <- match.arg( orientation )

  half.gap <- ( 1 - arrow.length ) / 2

  strip <- ggplot() +
    coord_cartesian( xlim = c( 0, 1 ), ylim = c( 0, 1 ), clip = "off" ) +
    theme_void()

  if ( orientation == "horizontal" ) {

    strip <- strip +
      annotate(
        "segment",
        x = half.gap, xend = 1 - half.gap, y = 0.5, yend = 0.5,
        arrow = arrow( ends = "last", length = unit( 0.1, "inches" ) )
      ) +
      annotate(
        "text", x = 0.5, y = 0.5, label = label, vjust = -0.8, size = text.size / 2.8
      )

  } else {

    strip <- strip +
      annotate(
        "segment",
        x = 0.5, xend = 0.5, y = half.gap, yend = 1 - half.gap,
        arrow = arrow( ends = "last", length = unit( 0.1, "inches" ) )
      ) +
      annotate(
        "text", x = 0.5, y = 0.5, label = label, angle = 90, vjust = -0.8,
        size = text.size / 2.8
      )
  }

  strip
}


#' @title Add Shared Axis Labels to a Row of Biplots
#'
#' @description
#' Wraps a row of biplots (typically a `cowplot::plot_grid()` row) with a
#' larger shared x-axis label strip below it and a larger shared y-axis
#' label strip to its left, each built by `build.axis.label.strip()`, in
#' place of relying on each individual biplot's own, smaller axis title.
#'
#' @importFrom cowplot plot_grid
#'
#' @param plot.row A single plot object (typically a `cowplot::plot_grid()`
#' row of biplots) to wrap.
#' @param x.lab Character, the shared x-axis label.
#' @param y.lab Character, the shared y-axis label.
#' @param x.strip.height Numeric, fraction of the wrapped plot's total
#' height given to the x-axis label strip. Default `0.12`.
#' @param y.strip.width Numeric, fraction of the wrapped plot's total width
#' given to the y-axis label strip. Default `0.08`.
#' @param text.size Numeric, font size for both labels. Default `16`.
#'
#' @seealso [build.axis.label.strip()]
#'
#' @return A cowplot object: `plot.row`, wrapped with the shared axis label
#' strips.
#'
#' @export

add.biplot.shared.axis.labels <- function(
    plot.row,
    x.lab,
    y.lab,
    x.strip.height = 0.12,
    y.strip.width = 0.08,
    text.size = 16
) {

  x.strip <- build.axis.label.strip( x.lab, orientation = "horizontal", text.size = text.size )
  y.strip <- build.axis.label.strip( y.lab, orientation = "vertical", text.size = text.size )

  with.x.label <- plot_grid(
    plot.row, x.strip, ncol = 1, rel_heights = c( 1 - x.strip.height, x.strip.height )
  )

  plot_grid(
    y.strip, with.x.label, nrow = 1, rel_widths = c( y.strip.width, 1 - y.strip.width )
  )
}
