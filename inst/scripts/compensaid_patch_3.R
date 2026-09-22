# compensaid_patch_3.R
#
# Source AFTER compensaid_patch.R and compensaid_patch_2.R.
#
# Patch 2's caching didn't move .GetPopulations()'s self-time the way
# predicted (88.28s self / 37.25%, essentially flat vs patch 1's
# 86.88s/34.30%) -- the markernames()/exprs() re-extraction wasn't actually
# the dominant cost there, so removing it bought almost nothing. What DID
# work in patch 2 was replacing the dplyr filter/pull chain in
# CompensAID()'s SSI-matrix-fill step (filter_impl/filter_rows/filter_eval
# dropped out of the profile entirely). This patch applies that same
# dplyr-removal fix to a second, still-dplyr-driven hot spot, and adds one
# more genuinely redundant call this profiling run surfaced:
#
# 1) CompensAID()'s main loop calls .GetPopulations() twice per marker
#    pair -- once before .WithinLimit(), once after, always, regardless of
#    whether .WithinLimit() changed anything. .WithinLimit() only modifies
#    co.renew inside its two conditional blocks (population percentage
#    outside the 10-90% band). When neither fires, the co.renew coming out
#    of .WithinLimit() is identical to what went in, and the second
#    .GetPopulations() call recomputes the same four population matrices
#    the first call already produced, from the same inputs. This patch has
#    .WithinLimit() report whether it actually adjusted the cut-offs, and
#    CompensAID() skips the second call and reuses the first result when it
#    didn't.
#
# 2) .CalculateSSI() -- called at least once per segment per pair
#    (segment.value x n.pairs = 4 x 1,722 = 6,888 times on this file) --
#    does one line of arithmetic through dplyr::mutate()/dplyr::pull().
#    Because `si` keeps the dplyr grouping .EmptyMatrixInfo() puts on it
#    via group_by() for the whole run, every one of those calls pays
#    dplyr's per-call grouping overhead -- visible in the profile as
#    "grouped_df"/"compute_groups" (13.38s + 12.54s here). This patch
#    replaces it with the equivalent base R arithmetic.
#
# SAFETY: as with patch 2, this touches the exported CompensAID() function
# again (superseding patch 2's version, with the same two changes plus the
# ones below). Don't trust it on its own -- trust
# test_compensaid_patch_speedup.R's all.equal() check of $matrix and
# $matrixInfo between the stock and fully-patched runs. If that fails,
# don't use this patch.
#
# Usage: library(CompensAID); source("compensaid_patch.R"); source("compensaid_patch_2.R"); source("compensaid_patch_3.R")
# Applies for the current R session only.

if ( !requireNamespace( "CompensAID", quietly = TRUE ) ) {
  stop( "CompensAID must be installed and loadable before sourcing this patch." )
}

.patched.CalculateSSI <- function(si.input, primary.channel, secondary.channel, segment) {

  checkmate::assertDataFrame(si.input)
  checkmate::assertCharacter(primary.channel)
  checkmate::assertCharacter(secondary.channel)
  checkmate::assertNumeric(segment)

  row.idx <- which(si.input$primary.channel == primary.channel &
                      si.input$secondary.channel == secondary.channel &
                      si.input$segment == segment)

  ssi <- round((si.input$mfi.pos[row.idx] - si.input$mfi.neg[row.idx]) / (2 * si.input$sd.neg[row.idx]), digits = 2)

  return(ssi)
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

  adjusted <- FALSE

  percentage <- round((nrow(population[["secondary.negative"]])*100)/nrow(og))

  cp <- channel.of.marker[[primary]]
  cs <- channel.of.marker[[secondary]]

  if (percentage <= min | percentage >= max) {

    adjusted <- TRUE

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

    adjusted <- TRUE

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
                 si.dat = si.input,
                 adjusted = adjusted)

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

    # Only re-fetch populations when .WithinLimit() actually adjusted a
    # cut-off -- when it didn't, co.renew is unchanged and this would
    # recompute exactly what `pop` already holds.
    if (isTRUE(pop.limit$adjusted)) {
      pop <- .GetPopulations(og = ff,
                            primary = mc$primary.marker[i],
                            secondary = mc$secondary.marker[i],
                            co.input = co.renew,
                            sd.input = separation.distance,
                            exprs.mat = exprs.mat,
                            channel.of.marker = channel.of.marker)
    }

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

environment( .patched.CalculateSSI ) <- asNamespace( "CompensAID" )
environment( .patched.WithinLimit )  <- asNamespace( "CompensAID" )
environment( .patched.CompensAID )   <- asNamespace( "CompensAID" )

utils::assignInNamespace( ".CalculateSSI", .patched.CalculateSSI, ns = "CompensAID" )
utils::assignInNamespace( ".WithinLimit",  .patched.WithinLimit,  ns = "CompensAID" )
utils::assignInNamespace( "CompensAID",    .patched.CompensAID,   ns = "CompensAID" )

message( "CompensAID::.CalculateSSI() patched to base R (removes per-call dplyr grouping overhead); CompensAID()/.WithinLimit() patched to skip the redundant second .GetPopulations() call when cut-offs weren't adjusted." )
