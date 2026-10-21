.libPaths(c("/path/to/local/Rlibs/4.4", .libPaths()))

library(iCAMP)

comm <- readRDS("comm_16S_betaNTI.rds")
pd_names <- read.csv(
  "pdist_16S/pd.taxon.name.csv",
  row.names = 1,
  check.names = FALSE
)[, 1]

stopifnot(
  nrow(comm) == 103,
  ncol(comm) == 14060,
  identical(colnames(comm), pd_names)
)

set.seed(20260720)

bnti_999 <- iCAMP::bNTI.big(
  comm = comm,
  pd.desc = "pd.desc",
  pd.spname = pd_names,
  pd.wd = "pdist_16S",
  spname.check = TRUE,
  nworker = 4,
  memo.size.GB = 48,
  weighted = TRUE,
  rand = 999,
  output.dtail = FALSE,
  trace = TRUE
)

saveRDS(bnti_999, "betaNTI_weighted_rand999.rds")
write.csv(bnti_999, "betaNTI_weighted_rand999.csv")

cat("Full weighted βNTI run completed successfully\n")
print(dim(bnti_999))
print(summary(as.vector(bnti_999[upper.tri(bnti_999)])))
