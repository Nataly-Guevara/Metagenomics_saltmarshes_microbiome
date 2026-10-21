bnti <- readRDS("betaNTI_weighted_rand999.rds")
rc <- readRDS("RCbray_weighted_rand999.rds")

stopifnot(
  identical(dim(bnti), dim(rc)),
  identical(rownames(bnti), rownames(rc)),
  identical(colnames(bnti), colnames(rc))
)

idx <- which(upper.tri(bnti), arr.ind = TRUE)

out <- data.frame(
  Sample1 = rownames(bnti)[idx[,1]],
  Sample2 = colnames(bnti)[idx[,2]],
  betaNTI = bnti[idx],
  RCbray = rc[idx]
)

out$Process <- ifelse(
  out$betaNTI > 2, "Heterogeneous selection",
  ifelse(
    out$betaNTI < -2, "Homogeneous selection",
    ifelse(
      out$RCbray > 0.95, "Dispersal limitation",
      ifelse(
        out$RCbray < -0.95, "Homogenizing dispersal",
        "Undominated"
      )
    )
  )
)

write.csv(out, "assembly_process_pairwise.csv", row.names = FALSE)

print(table(out$Process))
print(round(100 * prop.table(table(out$Process)), 2))
