pairs <- read.csv(
  "bacterial_assembly_environment_master_table.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

scores <- read.csv(
  "environmental_PCA_sample_scores.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

pc_vars <- c("PC1", "PC2", "PC3")

# Match PCA scores to both samples in every pair
i1 <- match(pairs$Sample1, scores$SampleID)
i2 <- match(pairs$Sample2, scores$SampleID)

if (anyNA(i1) || anyNA(i2)) {
  stop("Some sample IDs in the pairwise table do not match the PCA score table.")
}

for (pc in pc_vars) {
  pairs[[paste0("Delta_", pc)]] <- abs(
    scores[[pc]][i1] - scores[[pc]][i2]
  )
}

# Combined Euclidean distance across PC1-PC3
pairs$PCA_distance_PC1_PC3 <- sqrt(
  pairs$Delta_PC1^2 +
  pairs$Delta_PC2^2 +
  pairs$Delta_PC3^2
)

write.csv(
  pairs,
  "bacterial_assembly_environment_PCA_master_table.csv",
  row.names = FALSE
)

# Spearman associations with betaNTI
distance_vars <- c(
  "Delta_PC1",
  "Delta_PC2",
  "Delta_PC3",
  "PCA_distance_PC1_PC3"
)

analyse_correlations <- function(d, comparison) {
  do.call(rbind, lapply(distance_vars, function(v) {
    test <- suppressWarnings(
      cor.test(
        d[[v]],
        d$betaNTI,
        method = "spearman",
        exact = FALSE
      )
    )

    data.frame(
      Comparison = comparison,
      Variable = v,
      N = sum(complete.cases(d[[v]], d$betaNTI)),
      Spearman_rho = unname(test$estimate),
      P_value = test$p.value
    )
  }))
}

results <- analyse_correlations(pairs, "All pairs")

results <- rbind(
  results,
  analyse_correlations(
    pairs[pairs$Same_site, ],
    "Within-site pairs"
  ),
  analyse_correlations(
    pairs[pairs$Site_comparison == "ACP-ACP", ],
    "ACP"
  ),
  analyse_correlations(
    pairs[pairs$Site_comparison == "BASS-BASS", ],
    "BASS"
  ),
  analyse_correlations(
    pairs[pairs$Site_comparison == "TRUG-TRUG", ],
    "TRUG"
  )
)

results$FDR <- ave(
  results$P_value,
  results$Comparison,
  FUN = function(z) p.adjust(z, method = "BH")
)

write.csv(
  results,
  "PCA_distances_betaNTI_spearman.csv",
  row.names = FALSE
)

cat("\nPCA-axis distances associated with betaNTI:\n")
print(results, row.names = FALSE)

# Median PCA-axis distance by assembly process
process_summary <- do.call(
  rbind,
  lapply(split(pairs, pairs$Process), function(d) {
    data.frame(
      Process = unique(d$Process),
      N = nrow(d),
      Median_Delta_PC1 = median(d$Delta_PC1, na.rm = TRUE),
      Median_Delta_PC2 = median(d$Delta_PC2, na.rm = TRUE),
      Median_Delta_PC3 = median(d$Delta_PC3, na.rm = TRUE),
      Median_PCA_distance = median(
        d$PCA_distance_PC1_PC3,
        na.rm = TRUE
      )
    )
  })
)

write.csv(
  process_summary,
  "PCA_distance_by_assembly_process.csv",
  row.names = FALSE
)

cat("\nMedian PCA distances by assembly process:\n")
print(process_summary, row.names = FALSE)

# Process proportions across combined PCA-distance deciles
pairs$PCA_distance_decile <- cut(
  pairs$PCA_distance_PC1_PC3,
  breaks = quantile(
    pairs$PCA_distance_PC1_PC3,
    probs = seq(0, 1, 0.1),
    na.rm = TRUE
  ),
  include.lowest = TRUE,
  labels = paste0("D", 1:10)
)

process_deciles <- as.data.frame(
  table(pairs$PCA_distance_decile, pairs$Process)
)

names(process_deciles) <- c(
  "PCA_distance_decile",
  "Process",
  "Count"
)

process_deciles$Proportion <- ave(
  process_deciles$Count,
  process_deciles$PCA_distance_decile,
  FUN = function(z) z / sum(z)
)

write.csv(
  process_deciles,
  "assembly_process_by_PCA_distance_decile.csv",
  row.names = FALSE
)

cat("\nAssembly-process proportions by PCA-distance decile:\n")
print(process_deciles[process_deciles$Count > 0, ], row.names = FALSE)
