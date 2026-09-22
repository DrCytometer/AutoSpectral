# compensaid_patch.R
#
# Runtime patch for CompensAID::.UpdateMatrixInfo().
#
# Confirmed root cause (from CompensAID's real source in UpdateMatrixInfo.R,
# cross-checked against summaryRprof() output on the 42-channel reproducing
# file): for every marker pair, .UpdateMatrixInfo() calls the internal
# helper .GetSegment() four separate times per segment -- once each to pull
# out "min", "max", "count", and "mfi.pos" -- even though a single
# .GetSegment() call already returns a one-row data.frame containing all
# four values. .GetSegment() re-filters the entire `primary.positive`
# population matrix from scratch and recomputes a median on every call, so
# this quadruples the cost of the single most expensive step in CompensAID
# for no benefit: 16 calls to .GetSegment() per marker pair instead of 4.
#
# This patch replaces .UpdateMatrixInfo() in the loaded CompensAID namespace
# with a version that calls .GetSegment() once per segment and reuses its
# four returned columns. The values assigned to si.input are identical to
# the original -- this deduplicates identical, deterministic calls, it does
# not change what gets computed.
#
# Usage: library(CompensAID); source("compensaid_patch.R")
# Applies for the current R session only (assignInNamespace() does not
# persist across sessions or touch the installed package).

if ( !requireNamespace( "CompensAID", quietly = TRUE ) ) {
  stop( "CompensAID must be installed and loadable before sourcing this patch." )
}

.patched.UpdateMatrixInfo <- function(si.input, range.input = NULL, primary, secondary, output, rv.input, population = NULL) {

  checkmate::assertDataFrame(si.input)
  checkmate::assertCharacter(primary)
  checkmate::assertCharacter(secondary)
  checkmate::assertCharacter(output)
  checkmate::assertNumeric(rv.input)

  aditVal <- !is.null(range.input)
  if (aditVal) {
    checkmate::assertNumeric(range.input)
    checkmate::assertList(population)
  }

  channel.primary <- si.input$primary.channel[si.input$primary.marker == primary & si.input$secondary.marker == secondary][1]
  channel.secondary <- si.input$secondary.channel[si.input$primary.marker == primary & si.input$secondary.marker == secondary][1]

  if (output == "No positive/negative population") {
    col.adjust <- c("segment.min", "segment.max", "event.count", "event.count.merge", "mfi.neg", "sd.neg", "mfi.pos", "ssi")
    si.input[si.input$primary.marker == primary & si.input$secondary.marker == secondary, col.adjust] <- NA
    si.input[si.input$primary.marker == primary & si.input$secondary.marker == secondary, "message"] <- "No positive/negative population"
  }

  if (output == "PASS") {

    si.input[si.input$primary.marker == primary & si.input$secondary.marker == secondary, "message"] <- "PASS"
    population.positive <- population$primary.positive

    for (s in seq_len(rv.input)) {

      segment.info <- .GetSegment(population = population.positive,
                                   primary.channel = channel.primary,
                                   secondary.channel = channel.secondary,
                                   segment = s,
                                   range.input = range.input)

      si.input$segment.min[si.input$primary.marker == primary &
                             si.input$secondary.marker == secondary &
                             si.input$segment == s] <- segment.info[,"min"]

      si.input$segment.max[si.input$primary.marker == primary &
                             si.input$secondary.marker == secondary &
                             si.input$segment == s] <- segment.info[,"max"]

      si.input$event.count[si.input$primary.marker == primary &
                             si.input$secondary.marker == secondary &
                             si.input$segment == s] <- segment.info[,"count"]

      si.input$mfi.pos[si.input$primary.marker == primary &
                         si.input$secondary.marker == secondary &
                         si.input$segment == s] <- segment.info[,"mfi.pos"]
    }

    population.negative <- population$primary.negative

    si.input$mfi.neg[si.input$primary.marker == primary &
                       si.input$secondary.marker == secondary] <- stats::median(population.negative[, channel.secondary])

    si.input$sd.neg[si.input$primary.marker == primary &
                      si.input$secondary.marker == secondary] <- stats::sd(population.negative[, channel.secondary])

    for (s in seq_len(rv.input)) {

      si.input$ssi[si.input$primary.marker == primary &
                     si.input$secondary.marker == secondary &
                     si.input$segment == s] <- .CalculateSSI(si.input = si.input,
                                                            primary.channel = channel.primary,
                                                            secondary.channel = channel.secondary,
                                                            segment = s)
    }
  }

  return(si.input)
}

# Rebind the patched function's enclosing environment to CompensAID's own
# namespace before installing it, so its unqualified calls to internal
# helpers (.GetSegment, .CalculateSSI) resolve the same way the original
# function's did -- without this, those calls would fail to find
# unexported objects.
environment( .patched.UpdateMatrixInfo ) <- asNamespace( "CompensAID" )
utils::assignInNamespace( ".UpdateMatrixInfo", .patched.UpdateMatrixInfo, ns = "CompensAID" )

message( "CompensAID::.UpdateMatrixInfo() patched: 16 redundant .GetSegment() calls per marker pair reduced to 4." )
