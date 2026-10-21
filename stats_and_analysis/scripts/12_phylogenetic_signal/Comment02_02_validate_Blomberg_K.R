# Independent validation of Blomberg's K
# Uses a phylogenetic covariance matrix and GLS.
# Does not use the original Python contrasts algorithm.
# Base R only; no external packages.

base <- "results/ecological_strategies/phylogenetic_signal/"

treefile <- paste0(
  "results/phylogeny_bacteria_352/",
  "bacteria_352_bac120.tree"
)

scores <- read.delim(
  paste0(base, "GRESP_six_strategies_352_MAGs.tsv"),
  check.names=FALSE
)

reported <- read.delim(
  paste0(base, "GRESP_phylogenetic_signal_K.tsv")
)

# Independently parse the Newick tree.
s <- paste(readLines(treefile, warn=FALSE), collapse="")

tokens <- regmatches(
  s,
  gregexpr(
    "'(?:[^']|'')*'|[(),:;]|[^(),:;[:space:]]+",
    s, perl=TRUE
  )
)[[1]]

at <- 1L
parents <- integer()
edges <- numeric()
labels <- character()
is_tip <- logical()
child_count <- integer()

new_node <- function() {
  id <- length(parents) + 1L
  parents <<- c(parents, NA_integer_)
  edges <<- c(edges, 0)
  labels <<- c(labels, "")
  is_tip <<- c(is_tip, TRUE)
  child_count <<- c(child_count, 0L)
  id
}

parse_node <- function() {
  id <- new_node()

  if (tokens[at] == "(") {
    is_tip[id] <<- FALSE
    at <<- at + 1L

    repeat {
      child <- parse_node()
      parents[child] <<- id
      child_count[id] <<- child_count[id] + 1L

      if (tokens[at] == ",") {
        at <<- at + 1L
      } else if (tokens[at] == ")") {
        at <<- at + 1L
        break
      } else {
        stop("Invalid Newick structure")
      }
    }

    if (!(tokens[at] %in% c(":",",",")",";"))) {
      labels[id] <<- tokens[at]
      at <<- at + 1L
    }

  } else {
    labels[id] <<- tokens[at]
    at <<- at + 1L
  }

  if (tokens[at] == ":") {
    at <<- at + 1L
    edges[id] <<- as.numeric(tokens[at])
    at <<- at + 1L
  }

  id
}

root <- parse_node()
stopifnot(tokens[at] == ";")

tip_ids <- which(is_tip)
tips <- labels[tip_ids]
tips <- gsub("''", "'", sub("'$", "", sub("^'", "", tips)))

stopifnot(
  length(tips) == 352,
  !anyDuplicated(tips),
  setequal(tips, scores$Final_MAG)
)

# Match original treatment of zero-length branches.
nonroot <- setdiff(seq_along(edges), root)

stopifnot(
  all(is.finite(edges[nonroot])),
  all(edges[nonroot] >= 0),
  any(edges[nonroot] > 0)
)

epsilon <- max(
  min(edges[nonroot][edges[nonroot] > 0]) * 1e-6,
  1e-12
)

zero <- nonroot[edges[nonroot] == 0]
edges[zero] <- epsilon

# Build shared-ancestry Brownian covariance matrix.
# Rows: tips; columns: branches.
B <- matrix(0, length(tip_ids), length(nonroot))

for (i in seq_along(tip_ids)) {
  v <- tip_ids[i]

  while (v != root) {
    B[i, match(v, nonroot)] <- sqrt(edges[v])
    v <- parents[v]
  }
}

C <- tcrossprod(B)
n <- nrow(C)

# GLS phylogenetic mean and expected Brownian ratio.
one <- rep(1, n)
weights <- solve(C, one)
denominator <- sum(weights)

expected_ratio <- (
  sum(diag(C)) - n / denominator
) / (n - 1)

stopifnot(expected_ratio > 0)

calculate_K <- function(y) {
  mu <- sum(weights * y) / denominator
  residual <- y - mu

  ordinary_ss <- sum(residual^2)
  phylo_ss <- as.numeric(
    crossprod(residual, solve(C, residual))
  )

  (ordinary_ss / phylo_ss) / expected_ratio
}

results <- data.frame(
  Strategy=reported$Strategy,
  Original_K=reported$Blombergs_K,
  Independent_K=NA_real_,
  Absolute_difference=NA_real_
)

for (i in seq_len(nrow(results))) {
  strategy <- results$Strategy[i]

  y <- scores[[strategy]][
    match(tips, scores$Final_MAG)
  ]

  k <- calculate_K(y)

  results$Independent_K[i] <- k
  results$Absolute_difference[i] <-
    abs(k - results$Original_K[i])
}

cat("\n=== INDEPENDENT BLOMBERG K VALIDATION ===\n")
cat("Matched bacterial MAGs:", n, "\n")
cat("Root descendants:", child_count[root], "\n")
cat("Zero-length branches regularised:", length(zero), "\n\n")

print(results, row.names=FALSE, digits=9)

cat(
  "\nMaximum K difference:",
  max(results$Absolute_difference), "\n"
)

# Independently verify FDR adjustment.
q <- p.adjust(reported$Permutation_p, method="BH")

cat(
  "Maximum FDR difference:",
  max(abs(q - reported$FDR_q)), "\n"
)

output <- paste0(
  base,
  "Comment02_02_independent_K_validation.tsv"
)

write.table(
  results, output,
  sep="\t", row.names=FALSE, quote=FALSE
)

if (max(results$Absolute_difference) < 1e-5 &&
    max(abs(q - reported$FDR_q)) < 1e-10) {
  cat("\nSUCCESS: Numerical K and FDR validation passed\n")
} else {
  cat("\nCHECK NEEDED: Numerical discrepancy detected\n")
}

cat("Output:", output, "\n")
