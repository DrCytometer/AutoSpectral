# save_pipeline_panels.R
#
# Extracts individual panels from the results list returned by
# plot.spectra.automated.steps(), plot.spectra.legacy.steps(), or
# plot.spectra.standard.workflow(), and saves each as its own image file
# in output.dir, named "<fluorophore>_<workflow>_<panel>.<file.type>".
#
# Usage:
#   res.auto <- plot.spectra.automated.steps( ... )
#   save.pipeline.panels( res.auto, workflow = "automated", output.dir = "figures/panels" )
#
#   res.legacy <- plot.spectra.legacy.steps( ... )
#   save.pipeline.panels( res.legacy, workflow = "legacy", output.dir = "figures/panels" )
#
#   res.std <- plot.spectra.standard.workflow( ... )
#   save.pipeline.panels( res.std, workflow = "standard", output.dir = "figures/panels" )

## Field name -> output label for each workflow's per-fluorophore results
## list. Order here is the order panels are saved in (not meaningful beyond
## that). "scatter.match.file" is handled separately since it's a path, not
## a plot object.
.panel.map <- list(
  automated = c(
    singlet.plot     = "A_singlet_gate",
    trace.plot       = "B_peak_trace",
    brightest.panel  = "C_brightest_selection",
    cosine.panel     = "D_cosine_filter",
    subtraction.plot = "F_final_spectrum"
  ),
  legacy = c(
    gate.panel          = "A_gate",
    af.panel            = "B_af_exclusion",
    scatter.match.panel = "C_scatter_match",
    rlm.panel           = "D_rlm_diagnostic",
    subtraction.plot    = "E_final_spectrum"
  ),
  standard = c(
    gate.panel        = "A_gate",
    selection.panel   = "B_brightest_selection",
    cosine.panel      = "C_cosine_filter",
    subtraction.plot  = "D_final_spectrum"
  )
)

## Extracts and saves individual panels from a results list produced by one
## of the three plot.spectra.*() helpers.
##
## @param results          The (invisibly returned) list from
##                          plot.spectra.automated.steps(),
##                          plot.spectra.legacy.steps(), or
##                          plot.spectra.standard.workflow(), keyed by
##                          fluorophore.
## @param workflow          One of "automated", "legacy", "standard".
## @param output.dir        Directory to write the individual panel files
##                           into. Created if it doesn't exist.
## @param file.type         Image format for ggsave()'d panels (default
##                           "png"). Ignored for the automated workflow's
##                           scatter-match panel, which is copied verbatim
##                           in its original (JPEG) format.
## @param width,height,dpi  Passed to ggplot2::ggsave() for every plot
##                           panel.
## @param include.composite If TRUE, also writes out each fluorophore's
##                           full composite figure (field "composite") as
##                           "<fluor>_<workflow>_composite.<file.type>".
## @param verbose           Print each file written.
##
## @return Invisibly, a character vector of the file paths written.
save.pipeline.panels <- function(
    results, workflow = c( "automated", "legacy", "standard" ),
    output.dir, file.type = "png",
    width = 6, height = 5, dpi = 300,
    include.composite = FALSE, verbose = TRUE
) {
  workflow <- match.arg( workflow )
  
  if ( !dir.exists( output.dir ) )
    dir.create( output.dir, recursive = TRUE )
  
  panel.fields <- .panel.map[[ workflow ]]
  written <- character( 0 )
  
  for ( fluor in names( results ) ) {
    
    entry <- results[[ fluor ]]
    fluor.safe <- gsub( "[^A-Za-z0-9._-]+", "_", fluor )
    
    for ( field in names( panel.fields ) ) {
      if ( is.null( entry[[ field ]] ) ) next
      
      out.file <- file.path(
        output.dir,
        sprintf( "%s_%s_%s.%s", fluor.safe, workflow, panel.fields[[ field ]], file.type )
      )
      
      ggplot2::ggsave(
        out.file, plot = entry[[ field ]],
        width = width, height = height, dpi = dpi, limitsize = FALSE
      )
      written <- c( written, out.file )
      if ( verbose ) message( "  Saved: ", out.file )
    }
    
    # -- automated workflow only: scatter-match panel is a saved JPEG path,
    # not a plot object, so it's copied rather than ggsave()'d
    if ( workflow == "automated" && !is.null( entry$scatter.match.file ) &&
         file.exists( entry$scatter.match.file ) ) {
      
      ext <- tools::file_ext( entry$scatter.match.file )
      out.file <- file.path(
        output.dir, sprintf( "%s_automated_E_scatter_match.%s", fluor.safe, ext )
      )
      file.copy( entry$scatter.match.file, out.file, overwrite = TRUE )
      written <- c( written, out.file )
      if ( verbose ) message( "  Saved: ", out.file )
    }
    
    if ( include.composite && !is.null( entry$composite ) ) {
      out.file <- file.path(
        output.dir, sprintf( "%s_%s_composite.%s", fluor.safe, workflow, file.type )
      )
      ggplot2::ggsave(
        out.file, plot = entry$composite,
        width = width * 1.5, height = height * length( panel.fields ),
        dpi = dpi, limitsize = FALSE
      )
      written <- c( written, out.file )
      if ( verbose ) message( "  Saved: ", out.file )
    }
  }
  
  invisible( written )
}