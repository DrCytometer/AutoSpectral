library(AutoSpectral)

# 1. Configuration paths
input_dir   <- "./Controls"  # Folder containing FlowJo exported FCS files ("export_...")
output_dir  <- "./Controls_renamed"   # Destination folder for updated FCS files
channel_csv <- "flowjo_channels.csv"  # CSV mapping old channel names to new names


if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 2. Load Channel / Parameter Mapping
channel_map   <- read.csv(channel_csv, stringsAsFactors = FALSE)
channel_map$FlowJo <- paste0("FJComp-", channel_map$FlowJo)
# Lookup: FlowJo channel -> Dye/fluorophore name
name_lookup <- setNames(channel_map$Dye, channel_map$FlowJo)

# 3. Locate FCS files to process
fcs_files <- list.files(input_dir, pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE)
message(sprintf("Found %d FCS file(s) to process.", length(fcs_files)))

# 4. Process files using AutoSpectral functions
for (fcs_file in fcs_files) {
  file_name <- basename(fcs_file)
  message(sprintf("Processing: %s ...", file_name))
  
  # Read matrix and keywords using AutoSpectral's fast reader
  fcs_data <- AutoSpectral::readFCS(fcs_file, return.keywords = TRUE)
  mat      <- fcs_data$data
  keys     <- fcs_data$keywords
  
  # Update column names on the expression matrix
  current_cols <- colnames(mat)
  renamed_cols <- current_cols
  
  # Determine total number of parameters from keywords ($PAR)
  n_par <- as.numeric(keys[["$PAR"]])
  
  for (i in seq_len(n_par)) {
    p_name_key <- paste0("$P", i, "N")
    p_desc_key <- paste0("$P", i, "S")
    
    orig_name <- keys[[p_name_key]]
    
    # Check if this channel exists in the mapping table
    if (!is.null(orig_name) && orig_name %in% names(name_lookup)) {
      new_name <- name_lookup[[orig_name]]
      
      # 1. Update matrix column header
      renamed_cols[current_cols == orig_name] <- new_name
      
      # 2. Update channel parameter name keyword ($PnN)
      keys[[p_name_key]] <- new_name
      
      # 3. Optionally align stain/description keyword ($PnS) if present
      if (!is.null(keys[[p_desc_key]])) {
        keys[[p_desc_key]] <- new_name
      }
    }
  }
  
  # Reassign updated column names back to matrix
  colnames(mat) <- renamed_cols
  
  # Write updated file back out with AutoSpectral's fast binary writer
  AutoSpectral::writeFCS(
    mat        = mat,
    keys       = keys,
    file.name  = file_name,
    output.dir = output_dir
  )
  
  message(sprintf("Successfully written to: %s", file.path(output_dir, file_name)))
}

message("Batch channel renaming complete!")
