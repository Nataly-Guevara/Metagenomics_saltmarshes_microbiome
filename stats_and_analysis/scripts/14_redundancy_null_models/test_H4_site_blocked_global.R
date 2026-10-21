# H4: Exact site-blocked permutation tests
# Multivariate functional prevalence (37 curated capabilities)
# No external R packages required.

input <- "results/ecological_strategies/functional_redundancy/functional_redundancy_20samples_37capabilities.tsv"
output <- "results/ecological_strategies/functional_redundancy/H4_site_blocked_permutation_tests.tsv"

d <- read.delim(input, check.names=FALSE)
d$Feature <- paste(d$GRESP, d$Capability, sep="::")

permute <- function(x) {
  if (length(x) == 1) return(list(x))
  unlist(lapply(seq_along(x), function(i)
    lapply(permute(x[-i]), function(z) c(x[i], z))
  ), recursive=FALSE)
}

test_group <- function(filter, habitat, comparison, treatments=NULL) {
  a <- d[d$Filter == filter & d$Habitat == habitat, ]

  if (!is.null(treatments)) {
    a <- a[a$Treatment %in% treatments, ]
    tab <- table(a$Site, a$Treatment)
    valid_sites <- rownames(tab)[rowSums(tab > 0) == length(treatments)]
    a <- a[a$Site %in% valid_sites, ]
  }

  if (is.null(treatments)) treatments <- sort(unique(a$Treatment))

  meta <- unique(a[,c("Sample","Site","Treatment")])
  features <- sort(unique(a$Feature))
  Y <- matrix(NA_real_, nrow(meta), length(features))

  Y[cbind(match(a$Sample, meta$Sample),
          match(a$Feature, features))] <- a$Weighted_prevalence_pct

  stopifnot(!anyNA(Y), length(features) == 37)

  # Arcsine-square-root transform of capability prevalence.
  Y <- asin(sqrt(pmin(1, pmax(0, Y / 100))))

  # Equalise feature contributions; exclude invariant features.
  keep <- apply(Y, 2, sd) > 1e-12
  Y <- scale(Y[,keep,drop=FALSE])

  site <- factor(meta$Site)
  X0 <- model.matrix(~site)
  rss0 <- sum(qr.resid(qr(X0), Y)^2)

  statistic <- function(labels) {
    treatment <- factor(labels, levels=treatments)
    X1 <- model.matrix(~site + treatment)
    rss1 <- sum(qr.resid(qr(X1), Y)^2)
    df1 <- qr(X1)$rank - qr(X0)$rank
    df2 <- nrow(Y) - qr(X1)$rank
    F <- ((rss0-rss1)/df1)/(rss1/df2)
    c(F=F, partial_R2=(rss0-rss1)/rss0)
  }

  original <- as.character(meta$Treatment)
  observed <- statistic(original)

  indices <- split(seq_len(nrow(meta)), meta$Site)
  options <- lapply(indices, function(idx) permute(original[idx]))
  grid <- expand.grid(lapply(options, seq_along),
                      KEEP.OUT.ATTRS=FALSE)

  null_F <- numeric(nrow(grid))
  for (i in seq_len(nrow(grid))) {
    labels <- original
    for (j in seq_along(indices)) {
      labels[indices[[j]]] <- options[[j]][[grid[i,j]]]
    }
    null_F[i] <- statistic(labels)["F"]
  }

  data.frame(
    Filter=filter,
    Habitat=habitat,
    Comparison=comparison,
    N_samples=nrow(meta),
    N_sites=length(indices),
    N_capabilities=sum(keep),
    Observed_F=unname(observed["F"]),
    Partial_R2=unname(observed["partial_R2"]),
    Exact_p=mean(null_F >= observed["F"] - 1e-9),
    N_permutations=length(null_F)
  )
}

results <- list()

for (f in c("Breadth10","Breadth50")) {
  for (h in c("BS","RZ")) {
    results[[length(results)+1]] <-
      test_group(f,h,"All_treatments")

    results[[length(results)+1]] <-
      test_group(f,h,"BESE_vs_BARE",c("BESE","BARE"))
  }
}

result <- do.call(rbind, results)

# FDR across the two primary omnibus tests.
result$Primary_FDR_q <- NA_real_
idx <- which(result$Filter=="Breadth10" &
             result$Comparison=="All_treatments")
result$Primary_FDR_q[idx] <- p.adjust(result$Exact_p[idx], "BH")

write.table(result, output, sep="\t",
            row.names=FALSE, quote=FALSE)

print(result, row.names=FALSE, digits=4)
cat("\nSUCCESS:", output, "\n")
