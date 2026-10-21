# 16S zeta diversity by site, treatment and habitat
# Uses original 103 individual community samples.
# Presence/absence of ASVs.
# Exact enumeration of sample combinations, orders 1-5.
# Base R only.

comm <- readRDS("betaNTI/comm_16S_betaNTI.rds")

meta <- read.csv(
  "betaNTI/metadata_16S.csv",
  fileEncoding="UTF-8-BOM",
  check.names=FALSE,
  stringsAsFactors=FALSE
)

names(meta)[1] <- "SampleID"

stopifnot(
  nrow(comm) == 103,
  ncol(comm) == 14060,
  setequal(rownames(comm), meta$SampleID)
)

meta <- meta[match(rownames(comm), meta$SampleID), ]

presence <- comm > 0

results <- list()

for (site in c("ACP", "BASS", "TRUG")) {
  for (habitat in c("BS", "RZ")) {
    for (treatment in c("NATURAL","BESE","BARE","DEGRADED")) {

      ix <- which(
        meta$Site == site &
        meta$SoilPortion == habitat &
        meta$Treatment == treatment
      )

      n <- length(ix)
      if (n == 0) next

      X <- presence[ix,,drop=FALSE]

      zeta1 <- mean(rowSums(X))

      for (order in seq_len(min(5,n))) {

        combinations <- combn(n, order)
        shared <- numeric(ncol(combinations))

        for (j in seq_len(ncol(combinations))) {
          selected <- combinations[,j]
          shared[j] <- sum(
            colSums(X[selected,,drop=FALSE]) == order
          )
        }

        results[[length(results)+1]] <- data.frame(
          Site=site,
          Habitat=habitat,
          Treatment=treatment,
          Zeta_order=order,
          N_samples=n,
          N_combinations=length(shared),
          Zeta_mean=mean(shared),
          Zeta_SD=if(length(shared)>1) sd(shared) else NA_real_,
          Zeta_normalised=mean(shared)/zeta1
        )
      }
    }
  }
}

results <- do.call(rbind,results)

output <- "betaNTI/zeta_16S_site_treatment_habitat.tsv"

write.table(
  results,output,
  sep="\t",row.names=FALSE,quote=FALSE
)

cat("\nSUCCESS: 16S zeta diversity calculated\n")
cat("Site-treatment-habitat groups:",
    length(unique(paste(
      results$Site,results$Habitat,results$Treatment
    ))),"\n")
cat("Result rows:",nrow(results),"\n")
cat("Output:",output,"\n")

cat("\n=== ZETA ORDER 2 ===\n")
print(
  subset(
    results,
    Zeta_order==2,
    select=c("Site","Habitat","Treatment",
             "N_samples","Zeta_mean","Zeta_normalised")
  ),
  row.names=FALSE
)
