.libPaths(c("/path/to/local/Rlibs/4.4", .libPaths()))

library(iCAMP)

comm <- readRDS("comm_16S_betaNTI.rds")

stopifnot(
  nrow(comm) == 103,
  ncol(comm) == 14060,
  all(rowSums(comm) > 0)
)

set.seed(20260721)

rc_out <- iCAMP::RC.pc(
  comm = comm,
  rand = 999,
  nworker = 4,
  memory.G = 24,
  weighted = TRUE,
  sig.index = "RC",
  detail.null = FALSE,
  output.bray = TRUE,
  silent = FALSE
)

saveRDS(rc_out, "RCbray_weighted_rand999_full.rds")
saveRDS(rc_out$index, "RCbray_weighted_rand999.rds")
write.csv(rc_out$index, "RCbray_weighted_rand999.csv")

cat("Full weighted RCbray run completed successfully\n")
print(dim(rc_out$index))
print(summary(rc_out$index[upper.tri(rc_out$index)]))
