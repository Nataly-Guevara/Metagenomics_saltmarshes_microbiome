# H2: Site-adjusted multivariate dispersion diagnostic
# Bray-Curtis distance to treatment centroids
# Site-restricted permutation of reduced-model residuals
# Base R only; no packages required.

input <- paste0(
  "results/functional_annotation/",
  "figure2_functional_reorganisation/",
  "KO_BrayCurtis_distances.tsv"
)

output <- paste0(
  "results/functional_annotation/",
  "figure2_functional_reorganisation/",
  "H2_KO_PERMDISP_site_blocked.tsv"
)

d <- read.delim(input, stringsAsFactors=FALSE)

permutations <- function(x) {
  if (length(x) == 1L)
    return(matrix(x, nrow=1L))

  do.call(rbind, lapply(seq_along(x), function(i) {
    cbind(x[i], permutations(x[-i]))
  }))
}

analyse <- function(habitat) {

  x <- subset(d, Habitat == habitat)

  meta <- unique(rbind(
    data.frame(
      Sample=x$Sample1,
      Site=x$Site1,
      Treatment=x$Treatment1
    ),
    data.frame(
      Sample=x$Sample2,
      Site=x$Site2,
      Treatment=x$Treatment2
    )
  ))

  stopifnot(!anyDuplicated(meta$Sample))

  meta <- meta[order(meta$Sample), ]
  n <- nrow(meta)

  stopifnot(n == if (habitat == "BS") 11 else 9)

  D <- matrix(
    NA_real_, n, n,
    dimnames=list(meta$Sample, meta$Sample)
  )

  diag(D) <- 0

  for (i in seq_len(nrow(x))) {
    D[x$Sample1[i], x$Sample2[i]] <- x$BrayCurtis[i]
    D[x$Sample2[i], x$Sample1[i]] <- x$BrayCurtis[i]
  }

  stopifnot(!anyNA(D))

  # Principal coordinates with Cailliez correction
  pc <- cmdscale(
    as.dist(D),
    k=n-1,
    eig=TRUE,
    add=TRUE
  )

  Y <- as.matrix(pc$points)

  # Distances to treatment centroids
  # Bias correction for small and unequal group sizes
  distances <- numeric(n)

  for (tr in unique(meta$Treatment)) {

    ix <- which(meta$Treatment == tr)
    nt <- length(ix)

    stopifnot(nt >= 2)

    group <- Y[ix, , drop=FALSE]
    centroid <- colMeans(group)

    deviations <- sweep(group, 2, centroid, "-")

    distances[ix] <- sqrt(rowSums(deviations^2)) *
      sqrt(nt / (nt - 1))
  }

  meta$Site <- factor(meta$Site)
  meta$Treatment <- factor(meta$Treatment)

  X0 <- model.matrix(~ Site, data=meta)
  X1 <- model.matrix(~ Site + Treatment, data=meta)

  q0 <- qr(X0)
  q1 <- qr(X1)

  df1 <- q1$rank - q0$rank
  df2 <- n - q1$rank

  Fstat <- function(y) {
    rss0 <- sum(qr.resid(q0, y)^2)
    rss1 <- sum(qr.resid(q1, y)^2)
    ((rss0-rss1)/df1)/(rss1/df2)
  }

  observed <- Fstat(distances)

  # Freedman-Lane residual permutations within sites
  fitted0 <- qr.fitted(q0, distances)
  residual0 <- distances - fitted0

  blocks <- split(seq_len(n), meta$Site)
  options <- lapply(blocks, permutations)

  grid <- expand.grid(
    lapply(options, function(z) seq_len(nrow(z))),
    KEEP.OUT.ATTRS=FALSE
  )

  null_F <- numeric(nrow(grid))

  for (i in seq_len(nrow(grid))) {

    yperm <- fitted0

    for (j in seq_along(blocks)) {
      ix <- blocks[[j]]
      permutation <- options[[j]][grid[i,j], ]
      yperm[ix] <- fitted0[ix] + residual0[permutation]
    }

    null_F[i] <- Fstat(yperm)
  }

  p <- mean(null_F >= observed - 1e-10)

  means <- tapply(
    distances,
    meta$Treatment,
    mean
  )

  cat("\n", habitat, "\n", sep="")
  cat("Site-adjusted dispersion F:", round(observed,4), "\n")
  cat("Restricted-permutation p:", round(p,6), "\n")
  cat("Permutations:", length(null_F), "\n")
  cat("Mean distance to centroid:\n")
  print(round(means,5))

  data.frame(
    Habitat=habitat,
    N=n,
    F=observed,
    Permutation_p=p,
    N_permutations=length(null_F),
    Treatment_mean_distances=paste(
      paste(names(means), round(means,6), sep="="),
      collapse=";"
    )
  )
}

results <- rbind(
  analyse("BS"),
  analyse("RZ")
)

write.table(
  results, output,
  sep="\t",
  row.names=FALSE,
  quote=FALSE
)

cat("\nSUCCESS\n")
cat("Output:", output, "\n")
