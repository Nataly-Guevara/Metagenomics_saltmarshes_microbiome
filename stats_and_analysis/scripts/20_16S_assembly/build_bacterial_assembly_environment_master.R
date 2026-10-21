# Build pairwise bacterial assembly–environment master table using base R

assembly_file <- "assembly_process_pairwise_with_metadata.csv"
metadata_file <- "metadata_16S.csv"
output_file   <- "bacterial_assembly_environment_master_table.csv"

assembly <- read.csv(
  assembly_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

metadata <- read.csv(
  metadata_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

env_vars <- c(
  "pH_avg",
  "Temperature_avg",
  "Redox_avg",
  "Moisture_avg",
  "Salinity",
  "NO3",
  "NO2",
  "PO4",
  "NH4"
)

# Basic checks
required_assembly <- c("Sample1", "Sample2", "betaNTI", "RCbray", "Process")
required_metadata <- c(
  "SampleID", "Site", "Treatment", "Block", "SoilPortion", env_vars
)

missing_assembly <- setdiff(required_assembly, names(assembly))
missing_metadata <- setdiff(required_metadata, names(metadata))

if (length(missing_assembly) > 0) {
  stop("Missing assembly columns: ", paste(missing_assembly, collapse = ", "))
}

if (length(missing_metadata) > 0) {
  stop("Missing metadata columns: ", paste(missing_metadata, collapse = ", "))
}

if (anyDuplicated(metadata$SampleID) > 0) {
  stop("Duplicated SampleID values exist in metadata_16S.csv")
}

# Check matching before merging
samples_in_assembly <- unique(c(assembly$Sample1, assembly$Sample2))
missing_from_metadata <- setdiff(samples_in_assembly, metadata$SampleID)
unused_metadata <- setdiff(metadata$SampleID, samples_in_assembly)

cat("Assembly pairs:", nrow(assembly), "\n")
cat("Unique assembly samples:", length(samples_in_assembly), "\n")
cat("Metadata samples:", nrow(metadata), "\n")
cat("Assembly samples missing from metadata:",
    length(missing_from_metadata), "\n")

if (length(missing_from_metadata) > 0) {
  cat("Missing sample names:\n")
  print(missing_from_metadata)
  stop("Cannot proceed until all assembly samples match the metadata.")
}

cat("Metadata samples not used in assembly:",
    length(unused_metadata), "\n")

# Standardise environmental variables across samples
metadata_scaled <- metadata

for (v in env_vars) {
  metadata_scaled[[paste0(v, "_z")]] <-
    as.numeric(scale(metadata_scaled[[v]]))
}

# Prepare metadata copies for Sample1 and Sample2
keep_cols <- c(
  "SampleID",
  "Site",
  "EcologicalStatus",
  "Treatment",
  "Block",
  "SoilPortion",
  env_vars,
  paste0(env_vars, "_z")
)

meta1 <- metadata_scaled[, keep_cols]
meta2 <- metadata_scaled[, keep_cols]

names(meta1) <- c(
  "Sample1",
  paste0(names(meta1)[-1], "_env1")
)

names(meta2) <- c(
  "Sample2",
  paste0(names(meta2)[-1], "_env2")
)

# Preserve original pairwise order
assembly$Pair_order <- seq_len(nrow(assembly))

master <- merge(
  assembly,
  meta1,
  by = "Sample1",
  all.x = TRUE,
  sort = FALSE
)

master <- merge(
  master,
  meta2,
  by = "Sample2",
  all.x = TRUE,
  sort = FALSE
)

master <- master[order(master$Pair_order), ]

# Absolute environmental differences in original measurement units
for (v in env_vars) {
  master[[paste0("Delta_", v)]] <- abs(
    master[[paste0(v, "_env1")]] -
    master[[paste0(v, "_env2")]]
  )
}

# Differences in standard-deviation units
for (v in env_vars) {
  master[[paste0("Delta_z_", v)]] <- abs(
    master[[paste0(v, "_z_env1")]] -
    master[[paste0(v, "_z_env2")]]
  )
}

z_difference_columns <- paste0("Delta_z_", env_vars)

master$Environmental_variables_available <- rowSums(
  !is.na(master[, z_difference_columns, drop = FALSE])
)

# Euclidean environmental distance using all nine standardised variables.
# Rows missing any environmental variable receive NA.
complete_environment <- complete.cases(
  master[, z_difference_columns, drop = FALSE]
)

master$Environmental_distance_z <- NA_real_

master$Environmental_distance_z[complete_environment] <- sqrt(
  rowSums(
    master[
      complete_environment,
      z_difference_columns,
      drop = FALSE
    ]^2
  )
)

# Remove temporary ordering column
master$Pair_order <- NULL

write.csv(master, output_file, row.names = FALSE)

cat("\nCreated:", output_file, "\n")
cat("Rows:", nrow(master), "\n")
cat("Columns:", ncol(master), "\n")
cat("Pairs with complete environmental distance:",
    sum(!is.na(master$Environmental_distance_z)), "\n")
cat("Pairs missing complete environmental data:",
    sum(is.na(master$Environmental_distance_z)), "\n")

cat("\nAssembly-process counts:\n")
print(table(master$Process, useNA = "ifany"))

cat("\nEnvironmental-distance summary:\n")
print(summary(master$Environmental_distance_z))
