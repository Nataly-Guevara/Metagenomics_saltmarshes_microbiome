# H2: Site-adjusted capability dispersion diagnostic
# Recalculates treatment centroids on every permutation
# Base R only; exact within-site permutations
# Observational-condition exchangeability remains an assumption

input <- paste0(
  "results/ecological_strategies/functional_redundancy/",
  "functional_redundancy_20samples_37capabilities.tsv"
)

output <- paste0(
  "results/ecological_strategies/functional_redundancy/",
  "H2_capability_dispersion.tsv"
)

d <- read.delim(input, stringsAsFactors=FALSE)
d$Feature <- paste(d$GRESP, d$Capability, sep="::")

perms <- function(x) {
  if (length(x) == 1L)
    return(matrix(x, nrow=1L))
  do.call(rbind, lapply(seq_along(x), function(i) {
    cbind(x[i], perms(x[-i]))
  }))
}

analyse <- function(filter, habitat) {

  a <- subset(d, Filter == filter & Habitat == habitat)
  meta <- unique(a[,c("Sample","Site","Treatment")])
  meta <- meta[order(meta$Sample),]

  features <- sort(unique(a$Feature))
  stopifnot(length(features) == 37)

  Y <- matrix(
    NA_real_, nrow(meta), length(features),
    dimnames=list(meta$Sample, features)
  )

  Y[cbind(
    match(a$Sample, meta$Sample),
    match(a$Feature, features)
  )] <- a$Weighted_prevalence_pct

  stopifnot(!anyNA(Y))

  Y <- asin(sqrt(pmin(pmax(Y,0),100)/100))
  keep <- apply(Y,2,sd) > 1e-10
  Y <- scale(Y[,keep,drop=FALSE])

  # Remove additive site effects before dispersion analysis.
  Z <- Y
  blocks <- split(seq_len(nrow(meta)),meta$Site)

  for (ix in blocks) {
    Z[ix,] <- sweep(
      Y[ix,,drop=FALSE],2,
      colMeans(Y[ix,,drop=FALSE]),"-"
    )
  }

  site <- factor(meta$Site)
  original <- as.character(meta$Treatment)
  levels_tr <- sort(unique(original))

  X0 <- model.matrix(~site)
  q0 <- qr(X0)

  statistic <- function(labels) {
    tr <- factor(labels,levels=levels_tr)
    distances <- numeric(nrow(Z))

    for (g in levels_tr) {
      ix <- which(labels == g)
      ng <- length(ix)
      stopifnot(ng >= 2)

      centroid <- colMeans(Z[ix,,drop=FALSE])
      residual <- sweep(
        Z[ix,,drop=FALSE],2,centroid,"-"
      )

      distances[ix] <-
        sqrt(rowSums(residual^2)) *
        sqrt(ng/(ng-1))
    }

    X1 <- model.matrix(~site+tr)
    q1 <- qr(X1)

    ss0 <- sum(qr.resid(q0,distances)^2)
    ss1 <- sum(qr.resid(q1,distances)^2)

    df1 <- q1$rank-q0$rank
    df2 <- length(distances)-q1$rank

    F <- ((ss0-ss1)/df1)/(ss1/df2)

    list(
      F=F,
      means=tapply(distances,tr,mean)
    )
  }

  observed <- statistic(original)

  options <- lapply(blocks,function(ix)
    perms(original[ix]))

  grid <- expand.grid(
    lapply(options,function(z)seq_len(nrow(z))),
    KEEP.OUT.ATTRS=FALSE
  )

  null <- numeric(nrow(grid))

  for (i in seq_len(nrow(grid))) {
    labels <- original

    for (j in seq_along(blocks)) {
      labels[blocks[[j]]] <-
        options[[j]][grid[i,j],]
    }

    null[i] <- statistic(labels)$F
  }

  p <- mean(null >= observed$F-1e-9)

  cat("\n",filter,habitat,"\n")
  cat("Dispersion F:",round(observed$F,4),"\n")
  cat("Exact permutation p:",round(p,6),"\n")
  cat("Permutations:",length(null),"\n")
  print(round(observed$means,4))

  data.frame(
    Filter=filter,
    Habitat=habitat,
    N_samples=nrow(Y),
    N_capabilities=sum(keep),
    Dispersion_F=observed$F,
    Exact_p=p,
    N_permutations=length(null),
    Mean_distances=paste(
      paste(names(observed$means),
            round(observed$means,5),sep="="),
      collapse=";"
    )
  )
}

results <- do.call(rbind,lapply(
  c("Breadth10","Breadth50"),function(f)
    do.call(rbind,lapply(c("BS","RZ"),function(h)
      analyse(f,h)))
))

write.table(
  results,output,sep="\t",
  row.names=FALSE,quote=FALSE
)

cat("\nSUCCESS\nOutput:",output,"\n")
