# H2: Site-adjusted multivariate capability analysis
# Exact within-site permutation test
# 37 curated capabilities; 10% primary, 50% sensitivity
# Base R only.

input <- paste0(
  "results/ecological_strategies/functional_redundancy/",
  "functional_redundancy_20samples_37capabilities.tsv"
)

output <- paste0(
  "results/ecological_strategies/functional_redundancy/",
  "H2_capability_global_site_blocked.tsv"
)

d <- read.delim(input, stringsAsFactors=FALSE)
d$Feature <- paste(d$GRESP, d$Capability, sep="::")

perms <- function(v) {
  if (length(v) == 1)
    return(matrix(v, nrow=1))
  do.call(rbind, lapply(seq_along(v), function(i)
    cbind(v[i], perms(v[-i]))))
}

analyse <- function(filter, habitat) {

  a <- subset(d, Filter == filter & Habitat == habitat)

  meta <- unique(a[, c("Sample","Site","Treatment")])
  meta <- meta[order(meta$Sample), ]
  features <- sort(unique(a$Feature))

  stopifnot(length(features) == 37)
  stopifnot(!anyDuplicated(
    paste(a$Sample, a$Feature, sep="::")
  ))

  Y <- matrix(
    NA_real_, nrow(meta), length(features),
    dimnames=list(meta$Sample,features)
  )

  Y[cbind(
    match(a$Sample,meta$Sample),
    match(a$Feature,features)
  )] <- a$Weighted_prevalence_pct

  stopifnot(!anyNA(Y))
  stopifnot(nrow(Y) == if (habitat=="BS") 11 else 9)

  # Transform percentages, then standardise capabilities
  # so widespread traits do not automatically dominate.
  Y <- asin(sqrt(pmin(pmax(Y, 0), 100) / 100))
  keep <- apply(Y,2,sd) > 1e-10
  Y <- scale(Y[,keep,drop=FALSE])

  meta$Site <- factor(meta$Site)
  meta$Treatment <- factor(meta$Treatment)

  X0 <- model.matrix(~Site,data=meta)
  q0 <- qr(X0)
  rss0 <- sum(qr.resid(q0,Y)^2)

  original <- as.character(meta$Treatment)
  lev <- levels(meta$Treatment)

  stat <- function(labels) {
    mm <- data.frame(
      Site=meta$Site,
      Treatment=factor(labels,levels=lev)
    )
    X1 <- model.matrix(~Site+Treatment,data=mm)
    q1 <- qr(X1)
    rss1 <- sum(qr.resid(q1,Y)^2)

    df1 <- q1$rank-q0$rank
    df2 <- nrow(Y)-q1$rank

    c(
      F=((rss0-rss1)/df1)/(rss1/df2),
      Partial_R2=(rss0-rss1)/rss0,
      Df1=df1,
      Df2=df2
    )
  }

  observed <- stat(original)

  blocks <- split(seq_len(nrow(meta)),meta$Site)
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
    null[i] <- stat(labels)["F"]
  }

  p <- mean(null >= observed["F"]-1e-9)

  result <- data.frame(
    Filter=filter,
    Habitat=habitat,
    N_samples=nrow(Y),
    N_capabilities=sum(keep),
    F=unname(observed["F"]),
    Partial_R2=unname(observed["Partial_R2"]),
    Df1=unname(observed["Df1"]),
    Df2=unname(observed["Df2"]),
    Exact_p=p,
    N_permutations=length(null)
  )

  cat(
    filter, habitat,
    "| F =",round(result$F,4),
    "| partial R2 =",round(result$Partial_R2,4),
    "| p =",round(p,6),
    "| permutations =",length(null),"\n"
  )

  result
}

results <- do.call(rbind,lapply(
  c("Breadth10","Breadth50"),function(f)
    do.call(rbind,lapply(c("BS","RZ"),function(h)
      analyse(f,h)))
))

results$Primary_FDR_q <- NA_real_
idx <- which(results$Filter=="Breadth10")
results$Primary_FDR_q[idx] <- p.adjust(
  results$Exact_p[idx],method="BH"
)

write.table(
  results,output,sep="\t",
  row.names=FALSE,quote=FALSE
)

cat("\nSUCCESS\nOutput:",output,"\n")
