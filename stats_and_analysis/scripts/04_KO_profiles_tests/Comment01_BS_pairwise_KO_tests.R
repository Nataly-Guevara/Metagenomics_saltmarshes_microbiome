# Comment 1: Site-matched bulk-soil KO comparisons
# Existing Bray-Curtis distances; base R
# Exact within-site treatment-label swaps
# Observational NAT/DEGR comparisons are exploratory.

input <- paste0(
  "results/functional_annotation/",
  "figure2_functional_reorganisation/",
  "KO_BrayCurtis_distances.tsv"
)

outdir <- paste0(
  "results/functional_annotation/",
  "figure2_functional_reorganisation/"
)

d <- read.delim(input, stringsAsFactors=FALSE)
d <- subset(d, Habitat=="BS")

meta <- unique(rbind(
  data.frame(Sample=d$Sample1, Site=d$Site1,
             Treatment=d$Treatment1),
  data.frame(Sample=d$Sample2, Site=d$Site2,
             Treatment=d$Treatment2)
))

stopifnot(nrow(meta)==11)
stopifnot(!anyDuplicated(meta$Sample))

D <- matrix(
  NA_real_, 11, 11,
  dimnames=list(meta$Sample,meta$Sample)
)
diag(D) <- 0

for (i in seq_len(nrow(d))) {
  D[d$Sample1[i],d$Sample2[i]] <- d$BrayCurtis[i]
  D[d$Sample2[i],d$Sample1[i]] <- d$BrayCurtis[i]
}
stopifnot(!anyNA(D))

comparisons <- list(
  c("NAT","BESE"),
  c("NAT","BARE"),
  c("BESE","BARE"),
  c("NAT","DEGR"),
  c("BESE","DEGR"),
  c("BARE","DEGR")
)

results <- list()
site_rows <- list()

for (pair in comparisons) {

  valid <- sort(unique(meta$Site[
    meta$Treatment==pair[1]
  ]))

  valid <- valid[vapply(valid,function(s)
    all(pair %in% meta$Treatment[meta$Site==s]),
    logical(1))]

  m <- subset(
    meta, Site %in% valid & Treatment %in% pair
  )
  m <- m[order(m$Site,m$Treatment),]

  stopifnot(nrow(m)==2*length(valid))

  A <- D[m$Sample,m$Sample,drop=FALSE]
  n <- nrow(A)

  site <- factor(m$Site)
  original <- as.character(m$Treatment)

  J <- diag(n)-matrix(1/n,n,n)
  G <- -0.5*J%*%(A^2)%*%J
  total <- sum(diag(G))

  modelSS <- function(X) {
    Q <- qr.Q(qr(X))
    sum(Q*(G%*%Q))
  }

  X0 <- model.matrix(~site)
  ss0 <- modelSS(X0)
  rank0 <- qr(X0)$rank

  statistic <- function(labels) {
    tr <- factor(labels,levels=pair)
    X1 <- model.matrix(~site+tr)
    ss1 <- modelSS(X1)

    df1 <- qr(X1)$rank-rank0
    df2 <- n-qr(X1)$rank

    c(
      F=((ss1-ss0)/df1)/((total-ss1)/df2),
      R2=(ss1-ss0)/total,
      Partial_R2=(ss1-ss0)/(total-ss0)
    )
  }

  observed <- statistic(original)

  blocks <- split(seq_len(n),site)

  # All exact swap combinations: 2^N_sites
  grid <- expand.grid(
    rep(list(c(FALSE,TRUE)),length(blocks)),
    KEEP.OUT.ATTRS=FALSE
  )

  null <- numeric(nrow(grid))

  for (i in seq_len(nrow(grid))) {
    labels <- original
    for (j in seq_along(blocks)) {
      if (grid[i,j]) {
        ix <- blocks[[j]]
        labels[ix] <- rev(original[ix])
      }
    }
    null[i] <- statistic(labels)["F"]
  }

  p <- mean(null >= observed["F"]-1e-9)

  name <- paste(pair,collapse="_vs_")

  site_dist <- sapply(valid,function(s) {
    samples <- m$Sample[m$Site==s]
    A[samples[1],samples[2]]
  })

  for (s in valid) {
    site_rows[[length(site_rows)+1]] <- data.frame(
      Comparison=name,
      Site=s,
      BrayCurtis=unname(site_dist[s])
    )
  }

  results[[length(results)+1]] <- data.frame(
    Comparison=name,
    N_sites=length(valid),
    Sites=paste(valid,collapse=","),
    Pseudo_F=unname(observed["F"]),
    R2_after_site=unname(observed["R2"]),
    Partial_R2=unname(observed["Partial_R2"]),
    Exact_p=p,
    N_permutations=length(null),
    Mean_pair_distance=mean(site_dist)
  )

  cat(
    name,
    "| sites:",paste(valid,collapse=","),
    "| F:",round(observed["F"],4),
    "| R2:",round(observed["R2"],4),
    "| p:",round(p,6),
    "| mean distance:",round(mean(site_dist),5),
    "\n"
  )
}

results <- do.call(rbind,results)
site_rows <- do.call(rbind,site_rows)

results$Exploratory_FDR_q <- p.adjust(
  results$Exact_p,method="BH"
)

write.table(
  results,
  paste0(outdir,"Comment01_BS_pairwise_tests.tsv"),
  sep="\t",quote=FALSE,row.names=FALSE
)

write.table(
  site_rows,
  paste0(outdir,"Comment01_BS_site_matched_distances.tsv"),
  sep="\t",quote=FALSE,row.names=FALSE
)

cat("\nSUCCESS: Comment 1 bulk-soil contrasts saved\n")
