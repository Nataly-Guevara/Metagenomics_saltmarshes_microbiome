.libPaths(c("/path/to/local/Rlibs/4.4", .libPaths()))

library(iCAMP)
library(ape)

# Existing validated inputs
comm <- readRDS("betaNTI/comm_16S_betaNTI.rds")

tree <- read.tree(
  "betaNTI/icamp_pilot/icamp_pilot.tree.rooted.nwk"
)

pd_names <- read.csv(
  "betaNTI/pdist_16S/pd.taxon.name.csv",
  row.names = 1,
  check.names = FALSE
)[, 1]

stopifnot(
  nrow(comm) == 103,
  ncol(comm) == 14060,
  identical(colnames(comm), pd_names),
  setequal(tree$tip.label, pd_names),
  is.rooted(tree)
)

# Separate directory: preserve all pilot results
out_dir <- "betaNTI/icamp_full_103"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir)

set.seed(20261010)

cat("Starting full iCAMP analysis\n")
cat("Samples:", nrow(comm), "\n")
cat("Taxa:", ncol(comm), "\n")

result <- iCAMP::icamp.big(
  comm = comm,
  tree = tree,
  pd.desc = "pd.desc",
  pd.spname = pd_names,
  pd.wd = normalizePath("betaNTI/pdist_16S"),
  rand = 999,
  prefix = "icamp_full_103",
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

saveRDS(
  result,
  file.path(out_dir, "icamp_full_103_result.rds")
)

cat("\nSUCCESS: Full iCAMP completed\n")
cat("Output:", out_dir, "\n")
