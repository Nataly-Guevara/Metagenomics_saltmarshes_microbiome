# H2: Restoration contrasts and leave-one-site-out sensitivity
# Exact enumeration of within-site KO PERMANOVA permutations
# Base R only; no additional packages.

input <- "results/functional_annotation/figure2_functional_reorganisation/KO_BrayCurtis_distances.tsv"
output <- "results/functional_annotation/figure2_functional_reorganisation/H2_KO_treatment_site_sensitivity.tsv"

d <- read.delim(input, stringsAsFactors=FALSE)

perms <- function(v) {
  if (length(v) == 1)
    return(matrix(v, nrow=1))
  do.call(rbind, lapply(seq_along(v), function(i)
    cbind(v[i], perms(v[-i]))))
}

analyse <- function(habitat, treatments, omit=NULL, name="") {
  x <- subset(d, Habitat == habitat)

  m <- unique(rbind(
    data.frame(Sample=x$Sample1, Site=x$Site1,
               Treatment=x$Treatment1),
    data.frame(Sample=x$Sample2, Site=x$Site2,
               Treatment=x$Treatment2)
  ))

  stopifnot(!anyDuplicated(m$Sample))

  m <- m[m$Treatment %in% treatments, ]

  if (!is.null(omit))
    m <- m[m$Site != omit, ]

  sites <- unique(m$Site)
  valid <- sites[sapply(sites, function(s)
    all(treatments %in% m$Treatment[m$Site == s]))]

  m <- m[m$Site %in% valid, ]
  m <- m[order(m$Sample), ]

  stopifnot(nrow(m) == length(valid) * length(treatments))

  n <- nrow(m)
  D <- matrix(NA_real_, n, n,
              dimnames=list(m$Sample, m$Sample))
  diag(D) <- 0

  for (i in seq_len(nrow(x))) {
    a <- x$Sample1[i]
    b <- x$Sample2[i]
    if (a %in% m$Sample && b %in% m$Sample) {
      D[a,b] <- x$BrayCurtis[i]
      D[b,a] <- x$BrayCurtis[i]
    }
  }
  stopifnot(!anyNA(D))

  J <- diag(n) - matrix(1/n,n,n)
  G <- -0.5 * J %*% (D^2) %*% J
  total <- sum(diag(G))

  m$Site <- factor(m$Site)
  X0 <- model.matrix(~Site, m)

  modelSS <- function(X) {
    Q <- qr.Q(qr(X))
    sum(Q * (G %*% Q))
  }

  ss0 <- modelSS(X0)
  rank0 <- qr(X0)$rank

  stat <- function(labels) {
    mm <- data.frame(
      Site=m$Site,
      Treatment=factor(labels, levels=treatments)
    )
    X1 <- model.matrix(~Site+Treatment, mm)
    ss1 <- modelSS(X1)
    df1 <- qr(X1)$rank - rank0
    df2 <- n - qr(X1)$rank

    c(
      F=((ss1-ss0)/df1)/((total-ss1)/df2),
      R2=(ss1-ss0)/total,
      Partial_R2=(ss1-ss0)/(total-ss0)
    )
  }

  original <- as.character(m$Treatment)
  observed <- stat(original)

  blocks <- split(seq_len(n),m$Site)
  options <- lapply(blocks, function(ix)
    perms(original[ix]))

  grid <- expand.grid(
    lapply(options, function(z) seq_len(nrow(z))),
    KEEP.OUT.ATTRS=FALSE
  )

  null <- numeric(nrow(grid))

  for (i in seq_len(nrow(grid))) {
    labels <- original
    for (j in seq_along(blocks))
      labels[blocks[[j]]] <- options[[j]][grid[i,j],]
    null[i] <- stat(labels)["F"]
  }

  p <- mean(null >= observed["F"] - 1e-9)

  paired_distances <- NA_character_

  if (length(treatments) == 2) {
    vals <- sapply(names(blocks), function(s) {
      ix <- blocks[[s]]
      D[m$Sample[ix[1]],m$Sample[ix[2]]]
    })
    paired_distances <- paste(
      paste(names(vals),round(vals,5),sep="="),
      collapse=";"
    )
  }

  result <- data.frame(
    Habitat=habitat,
    Analysis=name,
    Sites=paste(valid,collapse=","),
    N_samples=n,
    F=unname(observed["F"]),
    R2_after_site=unname(observed["R2"]),
    Partial_R2=unname(observed["Partial_R2"]),
    Exact_p=p,
    N_permutations=length(null),
    Within_site_pair_distances=paired_distances
  )

  cat(name,
      "| sites:",paste(valid,collapse=","),
      "| F:",round(result$F,4),
      "| R2:",round(result$R2_after_site,4),
      "| p:",round(p,6),
      "| permutations:",length(null),"\n")

  if (!is.na(paired_distances))
    cat("  Within-site Bray-Curtis:",paired_distances,"\n")

  result
}

cat("\n=== Direct BESE versus BARE ===\n")

results <- list(
  analyse("RZ",c("BARE","BESE"),name="RZ_BESE_vs_BARE"),
  analyse("BS",c("BARE","BESE"),name="BS_BESE_vs_BARE")
)

cat("\n=== RZ omnibus leave-one-site-out ===\n")

for (site in c("ACP","BASS","TRUG")) {
  results[[length(results)+1]] <- analyse(
    "RZ",c("NAT","BESE","BARE"),
    omit=site,
    name=paste0("RZ_without_",site)
  )
}

result <- do.call(rbind,results)

write.table(result,output,sep="\t",
            row.names=FALSE,quote=FALSE)

cat("\nSUCCESS\nOutput:",output,"\n")
