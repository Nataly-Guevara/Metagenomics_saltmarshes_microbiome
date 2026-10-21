x <- read.csv("assembly_process_pairwise.csv", check.names = FALSE)

parse_sample <- function(samples) {
  parts <- strsplit(as.character(samples), "_")

  out <- data.frame(
    Sample      = samples,
    SampleID    = sapply(parts, `[`, 1),
    Run         = sapply(parts, `[`, 2),
    Site        = sapply(parts, `[`, 3),
    Treatment   = sapply(parts, `[`, 4),
    Compartment = sapply(parts, `[`, 5),
    stringsAsFactors = FALSE
  )

  # Use clearer treatment labels
  treatment_names <- c(
    NAT  = "NATURAL",
    BESE = "BESE",
    BARE = "BARE",
    DEGR = "DEGRADED",
    DEG  = "DEGRADED"
  )

  out$Treatment_full <- unname(
    ifelse(
      out$Treatment %in% names(treatment_names),
      treatment_names[out$Treatment],
      out$Treatment
    )
  )

  out
}

samples <- sort(unique(c(x$Sample1, x$Sample2)))
metadata <- parse_sample(samples)

m1 <- metadata
names(m1) <- paste0(names(m1), "1")

m2 <- metadata
names(m2) <- paste0(names(m2), "2")

out <- merge(x, m1, by.x = "Sample1", by.y = "Sample1", all.x = TRUE)
out <- merge(out, m2, by.x = "Sample2", by.y = "Sample2", all.x = TRUE)

out$Same_site        <- out$Site1 == out$Site2
out$Same_treatment   <- out$Treatment_full1 == out$Treatment_full2
out$Same_compartment <- out$Compartment1 == out$Compartment2

out$Site_comparison <- ifelse(
  out$Same_site,
  out$Site1,
  paste(sort(c(out$Site1, out$Site2)), collapse = "-")
)

# Correct pair-specific site comparisons
out$Site_comparison <- mapply(
  function(a, b) paste(sort(c(a, b)), collapse = "-"),
  out$Site1, out$Site2
)

out$Treatment_comparison <- mapply(
  function(a, b) paste(sort(c(a, b)), collapse = "-"),
  out$Treatment_full1, out$Treatment_full2
)

out$Compartment_comparison <- mapply(
  function(a, b) paste(sort(c(a, b)), collapse = "-"),
  out$Compartment1, out$Compartment2
)

write.csv(
  metadata,
  "sample_metadata_from_names.csv",
  row.names = FALSE
)

write.csv(
  out,
  "assembly_process_pairwise_with_metadata.csv",
  row.names = FALSE
)

cat("\nSample metadata dimensions:\n")
print(dim(metadata))

cat("\nSites:\n")
print(table(metadata$Site, useNA = "ifany"))

cat("\nTreatments:\n")
print(table(metadata$Treatment_full, useNA = "ifany"))

cat("\nCompartments:\n")
print(table(metadata$Compartment, useNA = "ifany"))

cat("\nAssembly processes by same versus different site:\n")
print(addmargins(table(out$Process, out$Same_site)))

cat("\nAssembly processes by same versus different treatment:\n")
print(addmargins(table(out$Process, out$Same_treatment)))

cat("\nAssembly processes by same versus different compartment:\n")
print(addmargins(table(out$Process, out$Same_compartment)))
