# assign_gates.r

#' @title Assign Gates
#'
#' @description
#' Assigns gate labels to samples based on whether they share characteristics,
#' such as being the same for `viability.gate`, `large.gate` or `control.type`
#' in the control file and whether the universal negative used is the same. Updates
#' the supplied control table and returns it with gate assignments.
#'
#' @param control.table Dataframe or table of the control file, read in via
#' `define.flow.control()` and cleaned up by that function.
#' @param gating.system Character string selecting the automated gating system to
#' employ in defining the initial scatter gates for identifying cells in the FCS
#' files. Options are `landmarks` and `density`. The `density` option uses the
#' original gating system from AutoSpill, picking cell populations based on dense
#' regions on FSC/SSC. The default `landmarks` system picks out the cell regions
#' using the brightest events in the peak channel for the single-stained control(s).
#' This approach is generally more robust. A good way to use this is to utilize
#' `landmarks` in combination with specifying which single-stained control files
#' should be used for defining the gates. For example, abundant, bright cell
#' markers such as CD4 will reliably identify the lymphocyte region, as CD14 will
#' identify the monocyte region, etc. For more instructions, see the help pages
#' on GitHub.
#' @param gate Logical, default is `TRUE`, in which case, automated gating will
#' be performed. If `FALSE`, the FCS files will be imported without automatically
#' generated gates applied. That is, all data in the files will be used. This is
#' intended to allow the user to pre-gate the files in commercial software.
#' @param verbose Logical, default is `TRUE`. Set to `FALSE` to suppress messages.
#'
#' @seealso
#' * [tune.gate()]
#' * [define.flow.control()]
#'
#' @return A modified `control.table`, complete with gate.name and gate.define
#'
#' @export

assign.gates <- function(
    control.table,
    gating.system,
    gate,
    verbose = TRUE
  ) {
  # define gates needed
  if ( verbose ) message( "\033[34mDetermining the gates that will be needed \033[0m" )

  # fill in missing elements
  if ( !"gate.name" %in% colnames( control.table ) ) control.table$gate.name <- NA
  if ( !"gate.define" %in% colnames( control.table ) ) control.table$gate.define <- TRUE
  if ( all( is.na( control.table$gate.define ) ) ) {
    control.table$gate.define <- TRUE
  } else {
    control.table$gate.define[ is.na( control.table$gate.define ) ] <- FALSE
  }
  control.table$is.viability[ is.na( control.table$is.viability ) ] <- FALSE
  control.table$large.gate[ is.na( control.table$large.gate ) ] <- FALSE

  # universal.negative must be strictly filename or NA
  control.table$universal.negative[
    control.table$universal.negative == "" | is.na( control.table$universal.negative )
  ] <- NA

  # Identify indices for "Unstained/AF/Negative" samples
  neg.pattern <- "(^AF$)|(\\bAF\\b)|(\\bnegative\\b)|(\\bunstained\\b)"
  is.neg.sample <- grepl( neg.pattern, control.table$fluorophore, ignore.case = TRUE )

  get.gate.label <- function(vi, lg) {
    if (vi) return("viabilityGate")
    if (lg) return("largeGate")
    return("smallGate")
  }

  gate.types <- data.frame(
    negative    = control.table$universal.negative,
    type        = control.table$control.type,
    viability   = control.table$is.viability,
    large.gate  = control.table$large.gate,
    stringsAsFactors = FALSE
  )
  unique.gate.combos <- unique( gate.types )

  # Assign names to stained samples only
  for (i in which(!is.neg.sample)) {
    # If user provided a name, keep it.
    # Otherwise, generate one based on parameters.
    if (is.na(control.table$gate.name[i])) {
      match_idx <- which(
        ( ( gate.types$negative[i] %in% unique.gate.combos$negative ) |
            ( is.na( gate.types$negative[i] ) & is.na( unique.gate.combos$negative ) ) ) &
          ( gate.types$type[i]       == unique.gate.combos$type ) &
          ( gate.types$viability[i]  == unique.gate.combos$viability ) &
          ( gate.types$large.gate[i] == unique.gate.combos$large.gate )
      )[1]

      label <- get.gate.label(control.table$is.viability[i], control.table$large.gate[i])
      control.table$gate.name[i] <- paste0(label, "_", match_idx)
    }
  }

  # replicate negatives as needed per gate
  if ( verbose ) message( "\033[34mMatching and replicating negatives per gate \033[0m" )

  # Identify unique combinations of negative files and the gates they must belong to
  neg.map <- unique(
    control.table[
      !is.neg.sample & !is.na(control.table$universal.negative),
      c("universal.negative", "gate.name", "large.gate", "is.viability")
    ]
  )

  if (nrow(neg.map) > 0) {
    for (i in seq_len(nrow(neg.map))) {
      target.file <- neg.map$universal.negative[i]
      target.gate <- neg.map$gate.name[i]

      # Find original rows for this negative file
      neg.idx <- which(is.neg.sample & control.table$filename == target.file )

      if (length(neg.idx) > 0) {
        # Check if one of these negative rows already matches this gate
        if (!any(control.table$gate.name[neg.idx] == target.gate, na.rm = TRUE)) {

          unassigned <- neg.idx[is.na(control.table$gate.name[neg.idx])]

          if (length(unassigned) > 0) {
            # Assign the gate to the existing NA row
            idx <- unassigned[1]
            control.table$gate.name[idx] <- target.gate

            # The "AF" sample should never be modified as there are hard-coded uses of "AF"
            control.table$sample[idx] <- if ( control.table$fluorophore[ idx ]  == "AF" ) {
              "AF"
            } else {
              paste0( control.table$fluorophore[idx], "_", target.gate )
            }
          } else {
            # Replicate the row if all existing rows for this file are busy in other gates
            new.row <- control.table[neg.idx[1], ]
            new.row$gate.name <- target.gate
            new.row$large.gate <- neg.map$large.gate[i]
            new.row$is.viability <- neg.map$is.viability[i]
            new.row$fluorophore <- if ( grepl( "negative", new.row$fluorophore, ignore.case = TRUE ) ) {
              paste( new.row$fluorophore, target.gate )
            } else if ( new.row$fluorophore == "AF" ) {
              "AF"
            } else {
              paste( new.row$fluorophore, "Negative", target.gate )
            }
            new.row$sample <- new.row$fluorophore
            new.row$gate.define <- FALSE
            control.table <- rbind(control.table, new.row)
          }
        }
      }
      # update negative sample boolean--could be with c()
      is.neg.sample <- grepl(neg.pattern, control.table$fluorophore, ignore.case = TRUE)
      # update the map with the new rows
      neg.map <- unique(
        control.table[
          !is.neg.sample & !is.na(control.table$universal.negative),
          c("universal.negative", "gate.name", "large.gate", "is.viability")
        ]
      )

    }
  }

  # If a negative sample still has no gate.name (meaning no stained sample uses it as a universal negative)
  is.neg.sample <- grepl(neg.pattern, control.table$fluorophore, ignore.case = TRUE)
  orphans <- which(is.neg.sample & is.na(control.table$gate.name))

  if (length(orphans) > 0) {
    for (idx in orphans) {
      # Fallback to density gating nomenclature
      label <- get.gate.label(control.table$is.viability[idx], control.table$large.gate[idx])
      control.table$gate.name[idx] <- paste0(label, "_density_orphan")

      # The "AF" sample should never be modified as there are hard-coded uses of "AF"
      control.table$sample[idx] <- if ( control.table$fluorophore[ idx ]  == "AF" ) {
        "AF"
      } else {
        paste0( control.table$fluorophore[idx], "_", target.gate )
      }

      # Note: The gating loop later will call define.gate.density for this name
      # set gate.define = TRUE for these
      control.table$gate.define[idx] <- TRUE
    }
  }

  # Update universal.negative to point to the correct Sample Name within that gate
  control.table$universal.negative <- sapply(seq_len(nrow(control.table)), function(i) {
    neg.file <- control.table$universal.negative[i]
    if (is.na(neg.file)) return(NA)

    # Match the filename AND the gate name to ensure the correct "Negative" is picked
    match.row <- control.table[
      control.table$filename == neg.file & control.table$gate.name == control.table$gate.name[i],
    ]
    if (nrow(match.row) > 0) return(match.row$sample[1])

    # Fallback to first matching filename if no gate-specific match (safety)
    match.fallback <- control.table[control.table$filename == neg.file, ]
    if (nrow(match.fallback) > 0) return(match.fallback$sample[1])
    return(NA)
  })

  if ( gate ) {
    # check through unique combos and ensure that at least one stained sample is
    # set as gate.define = TRUE
    # if not, set all samples in the gate as TRUE and warn user
    unique.gates <- unique( control.table$gate.name )

    for ( g in unique.gates ) {
      # which samples use this gate?
      active.sample.idx <- which( control.table$gate.name == g )
      # how many stained samples do we have for this gate?
      stained.idx <- which( control.table$gate.name == g & !is.neg.sample )
      stained.n <- length( stained.idx )
      # how many are set to be used to define the gate?
      define.idx <- which( control.table$gate.define[ active.sample.idx ] )
      define.n <- length( define.idx )
      # which samples are both stained and being used to define the gate?
      stained.define.idx <- which(
        control.table$gate.name == g & !is.neg.sample & control.table$gate.define
      )

      # if none, set gate.define = TRUE for all samples--slightly redundant, but safer
      if ( stained.n == 0 & define.n == 0 ) {
        control.table$gate.define[ active.sample.idx ] <- TRUE
      } else if ( !length( stained.define.idx ) & gating.system == "landmarks" ) {
        # if some are stained, check that at least one must be TRUE for gate.define if using landmarks
        # if not, set TRUE for all and warn
        control.table$gate.define[ active.sample.idx ] <- TRUE

        warning(
          paste0(
            "For gate ", g, ", no stained samples have been marked as TRUE for `gate.define`.
            At least one stained sample must be used to define a gate using the landmarks system.
            All samples will be used to define the gate. Change the control file if this is not what you want."
          )
        )
      }

      # enforce consistency for large.gate and is.viability
      # sample.type consistency will be enforced via check.control.file()
      if ( any( control.table$large.gate[ active.sample.idx ] ) ) {
        control.table$large.gate[ active.sample.idx ] <- TRUE
      }
      if ( any( control.table$is.viability[ active.sample.idx ] ) ) {
        control.table$is.viability[ active.sample.idx ] <- TRUE
      }

    }
  }

  return( control.table )
}
