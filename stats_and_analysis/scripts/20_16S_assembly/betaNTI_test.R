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

bnti_test <- iCAMP::bNTI.big(
  comm = comm,
  pd.desc = "pd.desc",
  pd.spname = pd_names,
  pd.wd = "pdist_16S",
  spname.check = TRUE,
  nworker = 4,
  memo.size.GB = 48,
  weighted = TRUE,
  rand = 20,
  output.dtail = FALSE,
  trace = TRUE
)

saveRDS(bnti_test, "betaNTI_test_rand20.rds")
write.csv(bnti_test, "betaNTI_test_rand20.csv")

cat("βNTI test completed successfully\n")
print(dim(bnti_test))
print(summary(as.vector(bnti_test[upper.tri(bnti_test)])))
