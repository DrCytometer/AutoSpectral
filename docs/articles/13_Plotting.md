# 11 Plotting

## Plotting with AutoSpectral

In this article, let’s look at a couple of the functions that
AutoSpectral offers for plotting your spectra and data. AutoSpectral
isn’t intended to replace flow analysis software like FlowJo or FCS
Express, but it’s essential to have some options for visualizing the
outputs.

``` r
library( AutoSpectral )
```

Let’s load in some spectra, here from OMIP-102.

``` r
s8.fcs <- "~/AutoSpectral_data/S8_data/BP0323502_1.fcs"
chorus.spectra <- read.bd.spectra( s8.fcs )
chorus.spectra[ 1:6, 1:5 ]
#>                UV1 (375)-A UV2 (390)-A UV3 (420)-A UV4 (440)-A UV5 (460)-A
#> CD40 BUV395-A   0.15491615  0.43476247   1.0000000   0.7179330  0.36265173
#> CD3 BUV496-A    0.02629492  0.08357235   0.2908688   0.2602035  0.24447993
#> CD56 BUV563-A   0.04424246  0.10662948   0.2355437   0.1685517  0.08813200
#> CD141 BUV615-A  0.02715625  0.06347369   0.1510007   0.1212349  0.06775691
#> CD303 BUV661-A  0.01522760  0.04873425   0.1771819   0.1613903  0.09411419
#> CD86 BUV737-A   0.03180002  0.10011904   0.3424454   0.3125159  0.18939716
```

The basic tools are
[`spectral.trace()`](https://drcytometer.github.io/AutoSpectral/reference/spectral.trace.md),
which gives line graphs of the spectra for each fluorophore in the
matrix, and
[`spectral.heatmap()`](https://drcytometer.github.io/AutoSpectral/reference/spectral.heatmap.md),
which gives the same information as a heatmap. The function
[`create.heatmap()`](https://drcytometer.github.io/AutoSpectral/reference/create.heatmap.md)
is a more general purpose heatmap that offers more options, including
triangular output. There is also a function to plot raw data for cells
across all detectors, like you get during the unmixing process on most
spectral cytometers. This can be accessed from
[`spectral.ribbon.plot()`](https://drcytometer.github.io/AutoSpectral/reference/spectral.ribbon.plot.md).

See the details for these functions here:
[spectral.trace](https://drcytometer.github.io/AutoSpectral/reference/spectral.trace.html)
[spectral.heatmap](https://drcytometer.github.io/AutoSpectral/reference/spectral.heatmap.html)
[create.heatmap](https://drcytometer.github.io/AutoSpectral/reference/create.heatmap.html)
[spectral.ribbon](https://drcytometer.github.io/AutoSpectral/reference/spectral.ribbon.plot.html)

``` r
spectral.heatmap( chorus.spectra, title = "OMIP102", color.palette = "mako",
                  show.legend = FALSE )
```

OMIP-102 spectra: ![OMIP102
spectra](figures/OMIP102_spectral_heatmap.jpg) Note that you can set
`save = FALSE` to just have the plots appear in the Viewer, but the
saved versions are scaled automatically to try to adjust for the amount
of data present.

The other option is a set of traces for each fluorophore.

``` r
spectral.trace( chorus.spectra, title = "OMIP102", split.lasers = FALSE,
                show.legend = FALSE )
```

![Spectral Trace without legend](figures/OMIP102_trace_noLegend.jpg)

Spectral Trace without legend

You can add the legend back:

``` r
spectral.trace( chorus.spectra, title = "OMIP102", split.lasers = FALSE,
                show.legend = FALSE )
```

![Spectral Trace](figures/OMIP102.jpg)

Spectral Trace

It’s also better to split it up by laser when there are many
fluorophores:

``` r
spectral.trace( chorus.spectra, title = "OMIP102", split.lasers = TRUE )
```

![Spectral Trace](figures/OMIP102_by_laser.jpg) In this example, the
laser/channel data are out of order because OMIP-102 was produced on an
early version of the BD FACSDiscoverS8. This is now fixed, and the
channels should always be re-arranged from narrowest to longest emission
wavelength within a given laser.

We can create a cosine similarity matrix plot in one of two ways.

``` r
cosine.similarity.plot( chorus.spectra, title = "Cosine_OMIP102", 
                        output.dir = getwd() )
```

![Cosine Similarity](figures/Cosine_OMIP102.jpg)

Cosine Similarity

This way generates a square matrix (duplicated info across the diagonal)
showing the cosine similarity values between each pair. The mixing
matrix condition number is listed below for reference.

We can also call `create.heatmap` which is a more general purpose
function that’s more useful for things like plotting the “hotspot
matrix”. This gives you a bit more control, but there are aspects that
still need some work (sorry).

For this, let’s look at a small subset of the matrix (the first five
fluorophores). We’ll calculate the cosine similarity matrix by calling
`cosine.similarity` (hint, do the same with other functions to see other
aspects of your data).

``` r
matrix.subset <- chorus.spectra[ 1:5, ]
small.cosine.matrix <- cosine.similarity( matrix.subset )
create.heatmap( small.cosine.matrix, number.labels = TRUE,
                color.palette = "turbo", plot.dir = getwd() )
```

![Cosine heatmap1](figures/heatmap.jpg)

Cosine heatmap1

Let’s change a couple of things.

``` r
create.heatmap( small.cosine.matrix, 
                number.labels = TRUE,
                legend.label = "Cosine Similarity",
                triangular = TRUE,
                color.palette = "mako",
                plot.dir = getwd() )
```

![Cosine heatmap2](figures/heatmap1.jpg)

Cosine heatmap2

We have the option to pass min and max scale settings, which is useful
for comparing multiple sets of data on the same scale (not really
relevant in this case, though).

``` r
create.heatmap( small.cosine.matrix, 
                number.labels = FALSE,
                legend.label = "Cosine Similarity",
                title = "OMIP102_Cosine_Similarity",
                triangular = TRUE,
                fixed.scale = TRUE,
                scale.min = 0, scale.max = 2,
                color.palette = "plasma",
                plot.dir = getwd() )
```

![Cosine heatmap3](figures/OMIP102_Cosine_Similarity_heatmap.jpg)

Cosine heatmap3

Finally, we can make x-y biplots of FCS data, which can be useful for
quickly checking what’s going on. To be fair, this is where programs
with graphical user interfaces like FlowJo, FCS Express, SpectroFlo,
etc. really excel, so doing it here is in R not really the best option.

To start, we can load in an FCS file–let’s use the OMIP-102 file–and we
need to call up the AutoSpectral parameter list because I haven’t gotten
around to removing that dependency yet.

``` r
asp <- get.autospectral.param( cytometer = "s8" )
omip102.file <- "~/AutoSpectral_data/S8_data/BP0323502_1.fcs"
omip102.ff <- suppressWarnings( flowCore::read.FCS( omip102.file, 
                                          transformation = FALSE,
                                          truncate_max_range = FALSE,
                                          emptyValue = FALSE ) )
omip102.data <- flowCore::exprs( omip102.ff )
```

Now we can call the plotting function:

``` r
create.biplot( omip102.data,
               x.dim = "BUV395-A", y.dim = "BUV496-A",
               asp,
               title = "OMIP102_biplot1" )
```

![OMIP102 biplot](figures/OMIP102_biplot1.jpg)

OMIP102 biplot

Note that this is currently ungated. You can play around with the
automated gating functions in AutoSpectral if you want (these work okay
with OMIP-102), but you’ll be better off using a graphical interface.

We can change the plot limits, color palette, etc. I haven’t implemented
a rainbow palette like FlowJo has yet.

``` r
create.biplot( omip102.data,
               x.dim = "BUV395-A", y.dim = "BUV496-A",
               x.min = -15000, x.max = 1e6,
               color.palette = "inferno",
               asp, 
               title = "OMIP102_biplot2" )
```

![OMIP102 biplot](figures/OMIP102_biplot2.jpg)

OMIP102 biplot

Change the width basis parameter to compress the data more or less
around zero. As we’re using a replica of the FlowJo biexponential
transform in R, the width basis caps out at -1000. I’ve rigged it so
that putting in values less than -1000 will increase compression by
altering the positive log decades, the other parameter you can use in
FlowJo.

``` r
create.biplot( omip102.data,
               x.dim = "BUV395-A", y.dim = "BUV496-A",
               x.min = -15000, x.max = 1e6,
               y.width.basis = -250,
               color.palette = "inferno",
               asp,  
               title = "OMIP102_biplot3" )
```

![OMIP102 biplot](figures/OMIP102_biplot3.jpg)

OMIP102 biplot

The maximum number of points to be plotted is capped by default at 1e5
to speed up the plotting. Setting an arbitrarily big number will plot
everything.

``` r
create.biplot( omip102.data,
               x.dim = "BUV395-A", y.dim = "BUV496-A",
               x.min = -15000, x.max = 1e6,
               y.max = 2e6,
               max.points = 1e8,
               color.palette = "viridis",
               asp,   
               title = "OMIP102_biplot4" )
```

![OMIP102 biplot](figures/OMIP102_biplot4.jpg)

OMIP102 biplot
