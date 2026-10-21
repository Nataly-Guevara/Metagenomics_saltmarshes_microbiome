x <- read.csv(
  "bacterial_assembly_environment_master_table.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

delta_vars <- c(
  "Delta_pH_avg",
  "Delta_Temperature_avg",
  "Delta_Redox_avg",
  "Delta_Moisture_avg",
  "Delta_Salinity",
  "Delta_NO3",
  "Delta_NO2",
  "Delta_PO4",
  "Delta_NH4"
)

analyse_subset <- function(d, comparison_name) {
  out <- do.call(rbind, lapply(delta_vars, function(v) {
    test <- suppressWarnings(
      cor.test(d[[v]], d$betaNTI, method = "spearman", exact = FALSE)
    )

    data.frame(
      Comparison = comparison_name,
      Variable = v,
      N = sum(complete.cases(d[[v]], d$betaNTI)),
      Spearman_rho = unname(test$estimate),
      P_value = test$p.value
    )
  }))

  out$FDR <- p.adjust(out$P_value, method = "BH")
  out
}

results <- analyse_subset(x, "All pairs")

results <- rbind(
  results,
  analyse_subset(x[x$Same_site, ], "Within-site pairs"),
  analyse_subset(x[x$Site_comparison == "ACP-ACP", ], "ACP"),
  analyse_subset(x[x$Site_comparison == "BASS-BASS", ], "BASS"),
  analyse_subset(x[x$Site_comparison == "TRUG-TRUG", ], "TRUG")
)

results <- results[
  order(results$Comparison, -abs(results$Spearman_rho)),
]

write.csv(
  results,
  "environmental_drivers_betaNTI_spearman.csv",
  row.names = FALSE
)

cat("\nEnvironmental differences associated with betaNTI:\n")
print(results, row.names = FALSE)

cat("\nMedian environmental differences by assembly process:\n")

process_summary <- do.call(
  rbind,
  lapply(split(x, x$Process), function(d) {
    values <- sapply(delta_vars, function(v) median(d[[v]], na.rm = TRUE))
    data.frame(
      Process = unique(d$Process),
      t(values),
      check.names = FALSE
    )
  })
)

print(process_summary, row.names = FALSE)

write.csv(
  process_summary,
  "median_environmental_difference_by_process.csv",
  row.names = FALSE
)
