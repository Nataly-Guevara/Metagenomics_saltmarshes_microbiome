# Paired habitat analysis of reconstructed KO profiles.
# Eight complete site x condition pairs, 16 metagenomes.
# Exact site-cluster habitat-label permutations.
# Base R only.

input <- paste0(
  "results/functional_annotation/sample_KO_profiles/",
  "sample_KO_relative_MAG_abundance_weighted.tsv"
)
out <- paste0(
  "results/functional_annotation/",
  "figure2_functional_reorganisation/",
  "Habitat_paired_KO_analysis.tsv"
)

tab <- read.delim(input, check.names=FALSE)

parts <- strsplit(tab$Sample, "_")

meta <- data.frame(
  Sample=tab$Sample,
  Site=vapply(parts, function(x) x[2], character(1)),
  Treatment=vapply(parts, function(x) x[3], character(1)),
  Habitat=vapply(parts, function(x) x[4], character(1)),
  stringsAsFactors=FALSE
)

meta$Pair <- paste(meta$Site, meta$Treatment, sep="_")

counts <- table(meta$Pair, meta$Habitat)
valid <- rownames(counts)[
  counts[,"BS"]==1 & counts[,"RZ"]==1
]

m <- meta[meta$Pair %in% valid, ]
m <- m[order(m$Pair,m$Habitat), ]

stopifnot(
  nrow(m)==16,
  length(valid)==8,
  identical(
    as.integer(table(m$Site)[c("ACP","BASS","TRUG")]),
    c(6L,6L,4L)
  )
)

Y <- as.matrix(
  tab[match(m$Sample,tab$Sample),-1,drop=FALSE]
)
storage.mode(Y) <- "double"

stopifnot(!anyNA(Y),all(rowSums(Y)>0))

# Row-normalised KO Bray-Curtis dissimilarities
Y <- Y / rowSums(Y)
n <- nrow(Y)
D <- matrix(0,n,n)

for (i in seq_len(n-1)) {
  for (j in seq.int(i+1,n)) {
    v <- sum(abs(Y[i,]-Y[j,]))/2
    D[i,j] <- D[j,i] <- v
  }
}

J <- diag(n)-matrix(1/n,n,n)
G <- -0.5*J%*%(D^2)%*%J
total <- sum(diag(G))

modelSS <- function(X) {
  Q <- qr.Q(qr(X))
  sum(Q*(G%*%Q))
}

m$Pair <- factor(m$Pair)
ss0 <- modelSS(model.matrix(~Pair,data=m))

habitat_stat <- function(labels) {
  mm <- m
  mm$Habitat <- factor(labels,levels=c("BS","RZ"))

  ss1 <- modelSS(model.matrix(~Pair+Habitat,data=mm))
  df1 <- 1
  df2 <- n-9

  c(
    F=((ss1-ss0)/df1)/((total-ss1)/df2),
    R2=(ss1-ss0)/total,
    Partial_R2=(ss1-ss0)/(total-ss0)
  )
}

# Dispersion diagnostic in a Cailliez-corrected
# principal-coordinate representation.
pc <- cmdscale(
  as.dist(D),k=n-1,eig=TRUE,add=TRUE
)
coordinates <- as.matrix(pc$points)

dispersion <- function(labels) {
  distances <- function(h) {
    x <- coordinates[labels==h,,drop=FALSE]
    sqrt(rowSums(
      sweep(x,2,colMeans(x),"-")^2
    ))
  }

  c(
    BS=mean(distances("BS")),
    RZ=mean(distances("RZ"))
  )
}

original <- m$Habitat
observed <- habitat_stat(original)
observed_disp <- dispersion(original)
observed_d <- observed_disp["RZ"]-observed_disp["BS"]

sites <- c("ACP","BASS","TRUG")
flips <- expand.grid(
  rep(list(c(FALSE,TRUE)),length(sites)),
  KEEP.OUT.ATTRS=FALSE
)

null_F <- numeric(nrow(flips))
null_disp <- numeric(nrow(flips))

for (k in seq_len(nrow(flips))) {
  labels <- original

  for (j in seq_along(sites)) {
    if (flips[k,j]) {
      ix <- m$Site==sites[j]
      labels[ix] <- ifelse(
        original[ix]=="BS","RZ","BS"
      )
    }
  }

  null_F[k] <- habitat_stat(labels)["F"]

  d <- dispersion(labels)
  null_disp[k] <- d["RZ"]-d["BS"]
}

p_F <- mean(null_F >= observed["F"]-1e-10)
p_disp <- mean(
  abs(null_disp)>=abs(observed_d)-1e-10
)

result <- data.frame(
  N_metagenomes=n,
  N_complete_pairs=length(valid),
  N_sites=length(sites),
  Habitat_pseudoF=unname(observed["F"]),
  Habitat_R2_after_pair=unname(observed["R2"]),
  Habitat_partial_R2=unname(observed["Partial_R2"]),
  Habitat_exact_site_p=p_F,
  Mean_dispersion_BS=unname(observed_disp["BS"]),
  Mean_dispersion_RZ=unname(observed_disp["RZ"]),
  Dispersion_exact_site_p=p_disp,
  N_site_permutations=nrow(flips)
)

write.table(
  result,out,sep="\t",
  row.names=FALSE,quote=FALSE
)

cat("\n=== PAIRED HABITAT KO ANALYSIS ===\n")
cat("Complete pairs:",length(valid),"\n")
cat("Metagenomes:",n,"\n")
cat("Sites:",paste(sites,collapse=", "),"\n")
cat("Habitat pseudo-F:",round(observed["F"],4),"\n")
cat("Habitat R2 after pair:",
    round(observed["R2"],4),"\n")
cat("Habitat partial R2:",
    round(observed["Partial_R2"],4),"\n")
cat("Exact site-restricted p:",p_F,"\n")

cat("\n=== DISPERSION DIAGNOSTIC ===\n")
cat("BS mean distance:",
    round(observed_disp["BS"],4),"\n")
cat("RZ mean distance:",
    round(observed_disp["RZ"],4),"\n")
cat("Exact site-restricted p:",p_disp,"\n")

cat("\nSUCCESS\nOutput:",out,"\n")
cat("NOTE: Only 8 site-level swap configurations;\n")
cat("permutation inference assumes habitat-label\n")
cat("exchangeability within the matched design.\n")
