# Comment 1: Site-matched treatment x habitat analysis
# Functional profiles: 6,679 MAG-derived KOs
# Exact site-restricted permutations; base R only.
#
# Response = Hellinger(RZ) - Hellinger(BS)
# One paired response per Site x Treatment.
# Missing TRUG BARE BS is never imputed.
#
# IMPORTANT: Permutation inference involving NAT assumes
# exchangeability of observational and restored conditions.
# These p-values must therefore be treated as exploratory.

base <- "results/functional_annotation/"
ko_file <- paste0(
  base, "sample_KO_profiles/",
  "sample_KO_relative_MAG_abundance_weighted.tsv"
)
dist_file <- paste0(
  base, "figure2_functional_reorganisation/",
  "KO_BrayCurtis_distances.tsv"
)
outdir <- paste0(base, "figure2_functional_reorganisation/")

k <- read.delim(ko_file, check.names=FALSE)
d <- read.delim(dist_file, stringsAsFactors=FALSE)

X <- as.matrix(k[,-1,drop=FALSE])
rownames(X) <- k$Sample
storage.mode(X) <- "double"

stopifnot(nrow(X)==20, ncol(X)==6679)
stopifnot(all(is.finite(X)), all(X>=0))
stopifnot(all(rowSums(X)>0))

# Relative KO abundance followed by Hellinger transform
H <- sqrt(sweep(X, 1, rowSums(X), "/"))

meta <- unique(rbind(
  data.frame(
    Sample=d$Sample1, Site=d$Site1,
    Treatment=d$Treatment1, Habitat=d$Habitat
  ),
  data.frame(
    Sample=d$Sample2, Site=d$Site2,
    Treatment=d$Treatment2, Habitat=d$Habitat
  )
))

stopifnot(nrow(meta)==20)
stopifnot(!anyDuplicated(meta$Sample))

# Match BS and RZ at the same site and treatment
pairs <- list()

for (site in c("ACP","BASS","TRUG")) {
  for (tr in c("NAT","BESE","BARE")) {
    bs <- meta$Sample[
      meta$Site==site &
      meta$Treatment==tr &
      meta$Habitat=="BS"
    ]
    rz <- meta$Sample[
      meta$Site==site &
      meta$Treatment==tr &
      meta$Habitat=="RZ"
    ]

    if (length(bs)==1 && length(rz)==1) {
      pairs[[length(pairs)+1]] <- data.frame(
        Site=site,
        Treatment=tr,
        BS_sample=bs,
        RZ_sample=rz
      )
    }
  }
}

pairs <- do.call(rbind,pairs)
stopifnot(nrow(pairs)==8)

# Signed multivariate change from BS to RZ
Delta <- H[pairs$RZ_sample,,drop=FALSE] -
         H[pairs$BS_sample,,drop=FALSE]

pairs$Hellinger_BS_RZ <- sqrt(rowSums(Delta^2))

permute <- function(v) {
  if (length(v)==1)
    return(matrix(v,nrow=1))
  do.call(rbind,lapply(seq_along(v),function(i)
    cbind(v[i],permute(v[-i]))))
}

test <- function(treatments, name) {

  ix <- which(pairs$Treatment %in% treatments)
  m <- pairs[ix,,drop=FALSE]
  Y <- Delta[ix,,drop=FALSE]

  # Retain only sites containing both treatments
  # for pairwise BS-to-RZ interaction comparisons.
  if (length(treatments) == 2L) {
    valid_sites <- unique(m$Site)[
      vapply(unique(m$Site), function(z) {
        all(treatments %in% m$Treatment[m$Site == z])
      }, logical(1))
    ]

    retain <- m$Site %in% valid_sites
    m <- m[retain,,drop=FALSE]
    Y <- Y[retain,,drop=FALSE]

    stopifnot(
      nrow(m) == 2L * length(valid_sites)
    )
  }

  site <- factor(m$Site)
  original <- as.character(m$Treatment)

  X0 <- model.matrix(~site)
  q0 <- qr(X0)
  rss0 <- sum(qr.resid(q0,Y)^2)

  total <- sum(scale(Y,center=TRUE,scale=FALSE)^2)

  stat <- function(labels) {
    tr <- factor(labels,levels=treatments)
    X1 <- model.matrix(~site+tr)
    q1 <- qr(X1)

    rss1 <- sum(qr.resid(q1,Y)^2)
    df1 <- q1$rank-q0$rank
    df2 <- nrow(Y)-q1$rank

    effect <- rss0-rss1

    c(
      F=(effect/df1)/(rss1/df2),
      R2_total=effect/total,
      Partial_R2=effect/rss0,
      Df1=df1,
      Df2=df2
    )
  }

  observed <- stat(original)

  blocks <- split(seq_len(nrow(m)),m$Site)
  options <- lapply(blocks,function(z)
    permute(original[z]))

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
    Analysis=name,
    N_pairs=nrow(m),
    N_sites=nlevels(site),
    Df_treatment=unname(observed["Df1"]),
    Df_residual=unname(observed["Df2"]),
    Pseudo_F=unname(observed["F"]),
    R2_total=unname(observed["R2_total"]),
    Partial_R2=unname(observed["Partial_R2"]),
    Exact_p=p,
    N_permutations=length(null)
  )

  cat(
    name,
    "| pairs:",nrow(m),
    "| F:",round(result$Pseudo_F,4),
    "| R2:",round(result$R2_total,4),
    "| p:",round(p,6),
    "| permutations:",length(null),"\n"
  )

  result
}

cat("\n=== COMMENT 1: HABITAT INTERACTION ===\n")

results <- rbind(
  test(c("NAT","BESE","BARE"),
       "All_three_conditions_exploratory"),
  test(c("BESE","BARE"),
       "BESE_vs_BARE"),
  test(c("NAT","BESE"),
       "NAT_vs_BESE_exploratory"),
  test(c("NAT","BARE"),
       "NAT_vs_BARE_exploratory")
)

write.table(
  pairs,
  paste0(outdir,"Comment01_BS_RZ_paired_samples.tsv"),
  sep="\t",quote=FALSE,row.names=FALSE
)

write.table(
  results,
  paste0(outdir,"Comment01_KO_habitat_interaction.tsv"),
  sep="\t",quote=FALSE,row.names=FALSE
)

cat("\nSUCCESS: Comment 1 interaction analysis saved\n")
