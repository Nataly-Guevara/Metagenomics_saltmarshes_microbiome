.libPaths(c("/path/to/local/Rlibs/4.4", .libPaths()))

library(iCAMP)
library(ape)

# Existing 16S inputs
comm <- readRDS("betaNTI/comm_16S_betaNTI.rds")
tree <- read.tree("bact_fasttree.nwk")

pd_names <- read.csv(
  "betaNTI/pdist_16S/pd.taxon.name.csv",
  row.names = 1,
  check.names = FALSE
)[, 1]

stopifnot(
  nrow(comm) == 103,
  ncol(comm) == 14060,
  identical(colnames(comm), pd_names),
  setequal(tree$tip.label, pd_names)
)

# Technical pilot: 12 samples only
comm <- comm[seq_len(12), , drop = FALSE]

out_dir <- normalizePath("betaNTI/icamp_pilot")

set.seed(20261010)

result <- iCAMP::icamp.big(
  comm = comm,
  tree = tree,
  pd.desc = "pd.desc",
  pd.spname = pd_names,
  pd.wd = normalizePath("betaNTI/pdist_16S"),
  rand = 19,
  prefix = "icamp_pilot",
  ds = 0.2,
  bin.size.limit = 24,
  phylo.metric = "bMPD",
  sig.index = "Confidence",
  nworker = 4,
  memory.G = 40,
  rtree.save = TRUE,
  detail.save = TRUE,
  qp.save = TRUE,
  output.wd = out_dir
)

saveRDS(result, file.path(out_dir, "icamp_pilot_result.rds"))

cat("\nSUCCESS: iCAMP pilot completed\n")
cat("Samples:", nrow(comm), "\n")
cat("Output:", out_dir, "\n")
