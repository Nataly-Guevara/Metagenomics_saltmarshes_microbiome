x <- read.csv(
  "bacterial_assembly_environment_master_table.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

x$Env_distance_bin <- cut(
  x$Environmental_distance_z,
  breaks = quantile(
    x$Environmental_distance_z,
    probs = seq(0, 1, 0.1),
    na.rm = TRUE
  ),
  include.lowest = TRUE,
  labels = paste0("D", 1:10)
)

bins <- levels(x$Env_distance_bin)

beta_summary <- do.call(
  rbind,
  lapply(bins, function(b) {
    z <- x$betaNTI[x$Env_distance_bin == b]
    data.frame(
      Env_distance_bin = b,
      n = length(z),
      mean = mean(z, na.rm = TRUE),
      median = median(z, na.rm = TRUE),
      sd = sd(z, na.rm = TRUE)
    )
  })
)

process_counts <- as.data.frame(
  table(x$Env_distance_bin, x$Process)
)

names(process_counts) <- c(
  "Env_distance_bin",
  "Process",
  "Count"
)

process_counts$Proportion <- ave(
  process_counts$Count,
  process_counts$Env_distance_bin,
  FUN = function(z) z / sum(z)
)

write.csv(
  beta_summary,
  "betaNTI_by_environmental_distance_decile.csv",
  row.names = FALSE
)

write.csv(
  process_counts,
  "assembly_process_by_environmental_distance_decile.csv",
  row.names = FALSE
)

cat("\nβNTI by environmental-distance decile:\n")
print(beta_summary)

cat("\nAssembly-process proportions by decile:\n")
print(process_counts[process_counts$Count > 0, ])
