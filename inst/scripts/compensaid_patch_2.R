# compensaid_patch_2.R
#
# Second-round runtime patch for CompensAID. Source AFTER compensaid_patch.R
# (this patch touches different functions -- .GetPopulations, .WithinLimit,
# and the exported CompensAID() itself -- so the two stack without
# conflict).
#
# Root cause: .GetPopulations() re-derives the same two things from the
# flowFrame on every call --
#   cp <- names(flowCore::markernames(og))[flowCore::markernames(og) == primary]
#   cs <- names(flowCore::markernames(og))[flowCore::markernames(og) == secondary]
#   sub <- flowCore::exprs(og)[, c(cp, cs)]
# -- even though the flowFrame (`ff` inside CompensAID()) never changes
# after the one-time zero-row padding at the top of the function. Every
# call copies the FULL expression matrix out of the S4 object just to keep
# two columns, and re-runs the markernames() generic (a linear scan over
# all parameters) twice. .GetPopulations() is called 2-4 times per marker
# pair (twice in the main loop, plus up to two more inside .WithinLimit()
# when its percentage thresholds trigger), so for a 42-channel panel
# (1,722 ordered pairs) that's several thousand redundant full-matrix
# extractions of an object that is identical every single time.
#
# This patch precomputes the expression matrix and a marker->channel
# lookup ONCE in CompensAID(), and threads them down through
# .WithinLimit() into .GetPopulations(), which now indexes directly into
# the precomputed matrix instead of re-extracting it. Same subsetting
# logic, same four population matrices out -- this removes repeated work
# on an unchanging input, it does not change what gets computed.
#
# It also replaces the dplyr::filter()/dplyr::filter()/dplyr::pull() chain
# CompensAID() uses to pull the max-segment SSI value into the output
# matrix with the equivalent base R subsetting. Same result, without
# paying dplyr's per-call tidy-eval overhead 1,722 times.
#
# SAFETY: the functions below are hand-transcribed from a clone of the
# CompensAID source and modified. If the installed CompensAID differs from
# that source in these three functions, this patch could compute something
# different from the stock package. Don't trust this patch on its own --
# trust test_compensaid_patch_speedup.R's all.equal() check of $matrix and
# $matrixInfo between the stock and patched runs. If that check fails,
# don't use this patch; the mismatch means the installed version has
# diverged from what this was written against.
#
# Usage: library(CompensAID); source("compensaid_patch.R"); source("compensaid_patch_2.R")
# Applies for the current R session only.

if ( !requireNamespace( "CompensAID", quietly = TRUE ) ) {
  stop( "CompensAID must be installed and loadable before sourcing this patch." )
}

.patched.GetPopulations <- function(og, primary, secondary, co.input, sd.input, exprs.mat, channel.of.marker) {

  checkmate::assert(methods::is(og, "flowFrame"), "Object is not a flowFrame.")
  checkmate::assertCharacter(primary)
  checkmate::assertCharacter(secondary)
  checkmate::assertDataFrame(co.input)
  checkmate::assertNumeric(sd.input)

  populations <- list()

  cp <- channel.of.marker[[primary]]
  cs <- channel.of.marker[[secondary]]
  sub <- exprs.mat[, c(cp, cs)]

  co.s <- co.input$opt[co.input$channel == cs]
  co.p.n <- co.input$opt[co.input$channel == cp] - sd.input
  co.p.p <- co.input$opt[co.input$channel == cp] + sd.input

  populations[["secondary.negative"]] <- sub[sub[, cs] <= co.s, , drop = FALSE]
  sub.s <- populations[["secondary.negative"]]
  populations[["primary.negative"]] <- sub.s[sub.s[, cp] <= co.p.n, , drop = FALSE]
  populations[["primary.positive"]] <- sub.s[sub.s[, cp] >= co.p.p, , drop = FALSE]

  return(populations)
}

.patched.WithinLimit <- function(population, og, primary, secondary, min = 10, max = 90, si.input, sd.input, co.input, cp.value, exprs.mat, channel.of.marker) {

  checkmate::assertList(population)
  checkmate::assert(methods::is(og, "flowFrame"), "Object is not a flowFrame.")
  checkmate::assertCharacter(primary)
  checkmate::assertCharacter(secondary)
  checkmate::assertNumeric(min)
  checkmate::assertNumeric(max)
  checkmate::assertDataFrame(si.input)
  checkmate::assertNumeric(sd.input)
  checkmate::assertDataFrame(co.input)
  checkmate::assertNumeric(cp.value)

  percentage <- round((nrow(population[["secondary.negative"]])*100)/nrow(og))

  cp <- channel.of.marker[[primary]]
  cs <- channel.of.marker[[secondary]]

  if (percentage <= min | percentage >= max) {

    co.adjust <- flowDensity::deGate(og, channel = cs, all.cuts = TRUE, tinypeak.removal = 0.0001, verbose = FALSE, upper = TRUE)

    co.adjust <- .GetClosestLimit(old.limit = co.input$opt[co.input$channel == cs],
                                 new.limit = .GetClosestCenter(co.adjust, cp.value),
                                 center.plot = cp.value)

    si.input$secondary.cutoff[si.input$primary.marker == primary & si.input$secondary.marker == secondary] <- .GetClosestCenter(co.adjust, cp.value)
    co.input$opt[co.input$channel == cs] <- .GetClosestCenter(co.adjust, cp.value)

    pop <- .GetPopulations(og = og,
                          primary = primary,
                          secondary = secondary,
                          co.input = co.input,
                          sd.input = sd.input,
                          exprs.mat = exprs.mat,
                          channel.of.marker = channel.of.marker)
  }

  percentage.positive <- round((nrow(population[["primary.positive"]])*100)/nrow(population[["secondary.negative"]]))
  percentage.negative <- round((nrow(population[["primary.negative"]])*100)/nrow(population[["secondary.negative"]]))

  if (percentage.positive <= min | percentage.positive >= max | percentage.negative <= min | percentage.negative >= max) {

    co.adjust <- flowDensity::deGate(og, channel = cp, all.cuts = TRUE, tinypeak.removal = 0.0001, verbose = FALSE, upper = TRUE)

    co.adjust <- .GetClosestLimit(old.limit = co.input$opt[co.input$channel == cp],
                                 new.limit = .GetClosestCenter(co.adjust, cp.value),
                                 center.plot = cp.value)

    si.input$primary.cutoff.neg[si.input$primary.marker == primary & si.input$secondary.marker == secondary] <- .GetClosestCenter(co.adjust, cp.value) - sd.input
    si.input$primary.cutoff.pos[si.input$primary.marker == primary & si.input$secondary.marker == secondary] <- .GetClosestCenter(co.adjust, cp.value) + sd.input
    co.input$opt[co.input$channel == cp] <- .GetClosestCenter(co.adjust, cp.value)

    pop <- .GetPopulations(og = og,
                          primary = primary,
                          secondary = secondary,
                          co.input = co.input,
                          sd.input = sd.input,
                          exprs.mat = exprs.mat,
                          channel.of.marker = channel.of.marker)
  }

  output <- list(co.dat = co.input,
                 si.dat = si.input)

  return(output)
}

.patched.CompensAID <- function(ff, segment.value = 4, events.value = 50) {

  checkmate::assert(methods::is(ff, "flowFrame"), "Object is not a flowFrame.")
  checkmate::assertNumeric(segment.value)
  checkmate::assertNumeric(events.value)
  sprintf('Importing sample: %s', ff@description[["FILENAME"]]) |> ParallelLogger::logInfo()

  ff@exprs <- rbind(ff@exprs, rep(0, ncol(ff@exprs)))

  mc <- .GetMarkerCombinations(og = ff)

  center <- .DetermineCenter(og = ff)
  center.plot <- center$center.plot
  separation.distance <- center$separation.distance

  co <- .DensityGating(og = ff,
                      cp.value = center.plot)

  sm <- .EmptyMatrix(og = ff)
  si <- .EmptyMatrixInfo(og = ff,
                        rv.input = segment.value,
                        mc.input = mc,
                        co.input = co,
                        sd.input = separation.distance)

  # Precompute what .GetPopulations()/.WithinLimit() otherwise re-derive
  # from `ff` on every call, even though `ff` never changes below this
  # point.
  exprs.mat <- flowCore::exprs(ff)
  marker.vec <- flowCore::markernames(ff)
  channel.of.marker <- stats::setNames(names(marker.vec), unname(marker.vec))

  for (i in seq_len(nrow(mc))) {

    co.renew <- co

    pop <- .GetPopulations(og = ff,
                          primary = mc$primary.marker[i],
                          secondary = mc$secondary.marker[i],
                          co.input = co.renew,
                          sd.input = separation.distance,
                          exprs.mat = exprs.mat,
                          channel.of.marker = channel.of.marker)

    pop.limit <- .WithinLimit(population = pop,
                             og = ff,
                             primary = mc$primary.marker[i],
                             secondary = mc$secondary.marker[i],
                             min = 10, max = 90,
                             si.input = si,
                             sd.input = separation.distance,
                             co.input = co.renew,
                             cp.value = center.plot,
                             exprs.mat = exprs.mat,
                             channel.of.marker = channel.of.marker)

    co.renew <- pop.limit$co.dat
    si <- pop.limit$si.dat

    pop <- .GetPopulations(og = ff,
                          primary = mc$primary.marker[i],
                          secondary = mc$secondary.marker[i],
                          co.input = co.renew,
                          sd.input = separation.distance,
                          exprs.mat = exprs.mat,
                          channel.of.marker = channel.of.marker)

    if (.EventRequirement(pop$primary.negative, pop$primary.positive, events.value)) {

      si <- .UpdateMatrixInfo(si.input = si,
                             rv.input = segment.value,
                             primary = mc$primary.marker[i],
                             secondary = mc$secondary.marker[i],
                             output = "No positive/negative population")

    } else {

      range <- (max(pop$primary.positive) - si$primary.cutoff.pos[si$primary.marker == mc$primary.marker[i] & si$secondary.marker == mc$secondary.marker[i]][1])/segment.value
      si <- .UpdateMatrixInfo(si.input = si,
                             rv.input = segment.value,
                             range.input = range,
                             primary = mc$primary.marker[i],
                             secondary = mc$secondary.marker[i],
                             output = "PASS",
                             population = pop)
    }

    si <- .UpdateSegments(si.input = si,
                         primary = mc$primary.marker[i],
                         secondary = mc$secondary.marker[i],
                         ev.input = events.value,
                         rv.input = segment.value)

    si <- .UpdateSSI(si.input = si,
                    rv.input = segment.value,
                    primary = mc$primary.marker[i],
                    secondary = mc$secondary.marker[i],
                    population = pop)

    pc <- channel.of.marker[[mc$primary.marker[i]]]
    sc <- channel.of.marker[[mc$secondary.marker[i]]]

    if (all(is.na(si$ssi[si$primary.channel == pc & si$secondary.channel == sc]))) {
      sm[sc, pc] <- NA
    } else {
      si.sub <- si[si$primary.channel == pc & si$secondary.channel == sc & !is.na(si$ssi), , drop = FALSE]
      sm[sc, pc] <- si.sub$ssi[si.sub$segment == max(si.sub$segment)]
    }
  }

  result <- list("matrix" = sm,
                 "matrixInfo" = si)

  return(result)
}

environment( .patched.GetPopulations ) <- asNamespace( "CompensAID" )
environment( .patched.WithinLimit )    <- asNamespace( "CompensAID" )
environment( .patched.CompensAID )     <- asNamespace( "CompensAID" )

utils::assignInNamespace( ".GetPopulations", .patched.GetPopulations, ns = "CompensAID" )
utils::assignInNamespace( ".WithinLimit",    .patched.WithinLimit,    ns = "CompensAID" )
utils::assignInNamespace( "CompensAID",      .patched.CompensAID,     ns = "CompensAID" )

message( "CompensAID::CompensAID()/.GetPopulations()/.WithinLimit() patched: expression matrix and marker->channel lookup are now computed once per file instead of once per .GetPopulations() call, and the SSI-matrix fill step no longer goes through dplyr." )
