# H2: KO functional composition
# Exact site-restricted distance-based PERMANOVA
# Uses the existing Figure 2A Bray-Curtis distance matrix.
# Base R only. No external packages.
#
# IMPORTANT:
# Permutation p-values assume condition-label exchangeability
# within sites under the null. Natural and Degraded are
# observational conditions, not randomised treatments.

input <- paste0(
  "results/functional_annotation/",
  "figure2_functional_reorganisation/",
  "KO_BrayCurtis_distances.tsv"
)

output <- paste0(
  "results/functional_annotation/",
  "figure2_functional_reorganisation/",
  "H2_KO_exact_site_blocked_PERMANOVA.tsv"
)

d <- read.delim(input, stringsAsFactors=FALSE)

all_perms <- function(v) {
  if (length(v) == 1)
    return(matrix(v, nrow=1))

  do.call(rbind, lapply(seq_along(v), function(i) {
    cbind(v[i], all_perms(v[-i]))
  }))
}

load_habitat <- function(h) {
  x <- subset(d, Habitat == h)

  m <- unique(rbind(
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

  stopifnot(!anyDuplicated(m$Sample))

  m <- m[order(m$Sample), ]
  n <- nrow(m)

  stopifnot(n == if (h == "BS") 11 else 9)

  D <- matrix(
    NA_real_, n, n,
    dimnames=list(m$Sample, m$Sample)
  )
  diag(D) <- 0

  for (i in seq_len(nrow(x))) {
    D[x$Sample1[i], x$Sample2[i]] <- x$BrayCurtis[i]
    D[x$Sample2[i], x$Sample1[i]] <- x$BrayCurtis[i]
  }

  stopifnot(!anyNA(D))
  stopifnot(isTRUE(all.equal(D, t(D))))

  list(meta=m, D=D)
}

run_test <- function(h, name, keep=NULL, complete=FALSE) {
  obj <- load_habitat(h)
  m <- obj$meta

  if (!is.null(keep)) {
    m <- m[m$Treatment %in% keep, ]

    if (complete) {
      sites <- unique(m$Site)
      valid <- sites[sapply(sites, function(s) {
        all(keep %in% m$Treatment[m$Site == s])
      })]
      m <- m[m$Site %in% valid, ]
    }
  }

  m <- m[order(m$Sample), ]
  m$Site <- factor(m$Site)
  m$Treatment <- factor(m$Treatment)

  D <- obj$D[m$Sample, m$Sample, drop=FALSE]
  n <- nrow(D)

  J <- diag(n) - matrix(1/n, n, n)
  G <- -0.5 * J %*% (D^2) %*% J
  totalSS <- sum(diag(G))

  modelSS <- function(X) {
    Q <- qr.Q(qr(X))
    sum(Q * (G %*% Q))
  }

  X0 <- model.matrix(~Site, data=m)
  ss0 <- modelSS(X0)
  rank0 <- qr(X0)$rank

  original <- as.character(m$Treatment)
  lev <- levels(m$Treatment)

  statistic <- function(labels) {
    design <- data.frame(
      Site=m$Site,
      Treatment=factor(labels, levels=lev)
    )

    X1 <- model.matrix(
      ~Site + Treatment, data=design
    )

    rank1 <- qr(X1)$rank
    ss1 <- modelSS(X1)

    df1 <- rank1 - rank0
    df2 <- n - rank1

    if (df1 <= 0 || df2 <= 0)
      stop("Insufficient residual degrees of freedom")

    effectSS <- ss1 - ss0
    residualSS <- totalSS - ss1

    if (residualSS <= 0)
      stop("Non-positive residual sum of squares")

    c(
      F=(effectSS/df1)/(residualSS/df2),
      R2_total=effectSS/totalSS,
      Partial_R2=effectSS/(totalSS-ss0),
      df1=df1,
      df2=df2
    )
  }

  observed <- statistic(original)

  # Enumerate every permitted within-site label assignment.
  idx <- split(seq_len(n), m$Site)

  possibilities <- lapply(idx, function(z) {
    all_perms(original[z])
  })

  choices <- lapply(possibilities, function(p) {
    seq_len(nrow(p))
  })

  grid <- expand.grid(
    choices,
    KEEP.OUT.ATTRS=FALSE
  )

  null_F <- numeric(nrow(grid))

  for (i in seq_len(nrow(grid))) {
    labels <- original

    for (j in seq_along(idx)) {
      labels[idx[[j]]] <-
        possibilities[[j]][grid[i,j], ]
    }

    null_F[i] <- statistic(labels)["F"]
  }

  p <- mean(null_F >= observed["F"] - 1e-9)

  cat(
    h, name,
    "| n =", n,
    "| sites =", nlevels(m$Site),
    "| permutations =", length(null_F),
    "| F =", round(observed["F"], 4),
    "| R2 =", round(observed["R2_total"], 4),
    "| p =", round(p, 6),
    "\n"
  )

  data.frame(
    Habitat=h,
    Analysis=name,
    N_samples=n,
    N_sites=nlevels(m$Site),
    N_conditions=nlevels(m$Treatment),
    Df_treatment=unname(observed["df1"]),
    Df_residual=unname(observed["df2"]),
    Pseudo_F=unname(observed["F"]),
    R2_after_site=unname(observed["R2_total"]),
    Partial_R2=unname(observed["Partial_R2"]),
    Exact_permutation_p=unname(p),
    N_permutations=length(null_F)
  )
}

cat("\nH2: Exact site-blocked KO PERMANOVA\n\n")

results <- rbind(
  run_test("BS", "Primary_all_conditions"),
  run_test("RZ", "Primary_all_conditions"),

  run_test(
    "BS", "Sensitivity_complete_sites",
    c("NAT","BESE","BARE","DEGR"),
    complete=TRUE
  ),

  run_test(
    "BS", "Sensitivity_three_conditions",
    c("NAT","BESE","DEGR"),
    complete=TRUE
  )
)

# Confirm reproduction of existing Figure 2A R2 values.
stopifnot(
  abs(results$R2_after_site[1] - 0.2261273264) < 1e-5,
  abs(results$R2_after_site[2] - 0.182) < 0.005
)

results$Primary_FDR_q <- NA_real_

primary <- which(
  results$Analysis == "Primary_all_conditions"
)

results$Primary_FDR_q[primary] <- p.adjust(
  results$Exact_permutation_p[primary],
  method="BH"
)

write.table(
  results, output,
  sep="\t",
  row.names=FALSE,
  quote=FALSE
)

cat("\nSUCCESS\n")
cat("Output:", output, "\n")
