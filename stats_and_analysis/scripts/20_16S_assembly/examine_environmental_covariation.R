m <- read.csv(
  "metadata_16S.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

env_vars <- c(
  "pH_avg", "Temperature_avg", "Redox_avg", "Moisture_avg",
  "Salinity", "NO3", "NO2", "PO4", "NH4"
)

env <- m[, env_vars]

# Spearman correlations at the independent sample level
cor_all <- cor(env, method = "spearman", use = "pairwise.complete.obs")
write.csv(cor_all, "environmental_spearman_correlation_all_samples.csv")

cat("\nOverall environmental Spearman correlations:\n")
print(round(cor_all, 2))

# Correlations within each site
for (s in unique(m$Site)) {
  cs <- cor(
    env[m$Site == s, , drop = FALSE],
    method = "spearman",
    use = "pairwise.complete.obs"
  )

  write.csv(
    cs,
    paste0("environmental_spearman_correlation_", s, ".csv")
  )

  cat("\nEnvironmental correlations —", s, ":\n")
  print(round(cs, 2))
}

# PCA after standardising all environmental variables
pca <- prcomp(env, center = TRUE, scale. = TRUE)

variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)

pca_variance <- data.frame(
  PC = paste0("PC", seq_along(variance)),
  Variance_percent = variance,
  Cumulative_percent = cumsum(variance)
)

pca_loadings <- data.frame(
  Variable = rownames(pca$rotation),
  pca$rotation,
  row.names = NULL
)

pca_scores <- data.frame(
  SampleID = m$SampleID,
  Site = m$Site,
  Treatment = m$Treatment,
  SoilPortion = m$SoilPortion,
  pca$x,
  check.names = FALSE
)

write.csv(pca_variance, "environmental_PCA_variance.csv", row.names = FALSE)
write.csv(pca_loadings, "environmental_PCA_loadings.csv", row.names = FALSE)
write.csv(pca_scores, "environmental_PCA_sample_scores.csv", row.names = FALSE)

cat("\nPCA explained variance:\n")
print(round(pca_variance, 2))

cat("\nPCA loadings for PC1–PC3:\n")
print(round(pca$rotation[, 1:3], 3))
