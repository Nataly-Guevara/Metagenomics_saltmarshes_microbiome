# Comment 4: Environmental differences and microbial assembly
# Sample-level, site- and field-block-aware permutations
# Existing 999-randomisation betaNTI results
# Base R only. Exploratory observational associations.

set.seed(20261012)
nperm <- 999L

base <- "betaNTI/"
metadata <- read.csv(
  paste0(base, "metadata_16S.csv"),
  fileEncoding="UTF-8-BOM",
  check.names=FALSE,
  stringsAsFactors=FALSE
)

names(metadata)[1] <- "SampleID"

pairs <- read.csv(
  paste0(base, "assembly_process_pairwise_with_metadata.csv"),
  stringsAsFactors=FALSE
)

predictors <- c("Redox_avg", "Moisture_avg", "NH4")
metrics <- c("Signed_betaNTI", "Absolute_betaNTI")
sites <- c("ACP", "BASS", "TRUG")
habitats <- c("BS", "RZ")

stopifnot(nrow(metadata)==103, nrow(pairs)==5253)
stopifnot(!anyDuplicated(metadata$SampleID))

for (v in predictors)
  metadata[[v]] <- as.numeric(metadata[[v]])

stopifnot(all(is.finite(as.matrix(metadata[,predictors]))))

prepare <- function(site, habitat) {
  m <- metadata[
    metadata$Site==site &
    metadata$SoilPortion==habitat, ,
    drop=FALSE
  ]

  p <- pairs[
    pairs$Sample1 %in% m$SampleID &
    pairs$Sample2 %in% m$SampleID, ,
    drop=FALSE
  ]

  stopifnot(nrow(p)==choose(nrow(m),2))

  i <- match(p$Sample1,m$SampleID)
  j <- match(p$Sample2,m$SampleID)

  groups <- paste(
    pmin(m$Treatment[i],m$Treatment[j]),
    pmax(m$Treatment[i],m$Treatment[j]),
    sep="::"
  )

  y <- as.numeric(p$betaNTI)
  stopifnot(all(is.finite(y)))

  outcomes <- list(
    Signed_betaNTI=y,
    Absolute_betaNTI=abs(y)
  )

  residuals <- lapply(outcomes,function(z) {
    z - ave(z,groups,FUN=mean)
  })

  list(
    site=site,
    habitat=habitat,
    meta=m,
    i=i,
    j=j,
    groups=groups,
    residuals=residuals,
    n_pairs=nrow(p)
  )
}

# Permute environmental profiles, not pairwise observations.
# Natural plots are exchanged only with Natural plots.
# Restored block labels are exchanged synchronously
# across treatments, preserving field-block matching.
# Incomplete blocks remain fixed.

permuted_environment <- function(obj) {
  m <- obj$meta
  index <- seq_len(nrow(m))

  natural <- which(m$Treatment=="NATURAL")
  if (length(natural)>1)
    index[natural] <- sample(natural)

  restored <- if (obj$habitat=="BS")
    c("BESE","BARE","DEGRADED") else
    c("BESE","BARE")

  blocks <- sort(unique(
    m$Block[m$Treatment %in% restored]
  ))

  complete <- blocks[vapply(blocks,function(b) {
    all(restored %in% m$Treatment[m$Block==b])
  },logical(1))]

  if (length(complete)>1) {
    shuffled <- sample(complete)

    for (k in seq_along(complete)) {
      for (tr in restored) {
        target <- which(
          m$Block==complete[k] & m$Treatment==tr
        )
        source <- which(
          m$Block==shuffled[k] & m$Treatment==tr
        )

        stopifnot(length(target)==1,length(source)==1)
        index[target] <- source
      }
    }
  }

  m[index,predictors,drop=FALSE]
}

calculate <- function(obj, environment) {
  answer <- numeric(length(metrics)*length(predictors))
  names(answer) <- as.vector(outer(
    metrics,predictors,paste,sep="::"
  ))

  for (v in predictors) {
    z <- as.numeric(environment[[v]])
    difference <- abs(z[obj$i]-z[obj$j]) / sd(z)

    # Remove treatment-pair means before association testing.
    x <- difference -
      ave(difference,obj$groups,FUN=mean)

    for (metric in metrics) {
      y <- obj$residuals[[metric]]
      key <- paste(metric,v,sep="::")

      answer[key] <- if (sd(x)>0 && sd(y)>0)
        cor(x,y) else NA_real_
    }
  }

  answer
}

datasets <- list()
for (h in habitats) {
  for (s in sites) {
    key <- paste(h,s,sep="::")
    datasets[[key]] <- prepare(s,h)
  }
}

keys <- as.vector(outer(metrics,predictors,paste,sep="::"))

observed_sites <- list()
for (name in names(datasets)) {
  obj <- datasets[[name]]
  observed_sites[[name]] <- calculate(
    obj,obj$meta[,predictors,drop=FALSE]
  )
}

observed <- c()
for (h in habitats) {
  values <- do.call(rbind,lapply(sites,function(s)
    observed_sites[[paste(h,s,sep="::")]]))
  observed <- c(observed,colMeans(values))
}

names(observed) <- unlist(lapply(
  habitats,function(h)paste(h,keys,sep="::")
))

stopifnot(all(is.finite(observed)))

null <- matrix(
  NA_real_,nperm,length(observed),
  dimnames=list(NULL,names(observed))
)

for (iteration in seq_len(nperm)) {
  values <- c()

  for (h in habitats) {
    site_values <- do.call(rbind,lapply(sites,function(s) {
      obj <- datasets[[paste(h,s,sep="::")]]
      calculate(obj,permuted_environment(obj))
    }))

    values <- c(values,colMeans(site_values))
  }

  null[iteration,] <- values
}

stopifnot(all(is.finite(null)))

results <- list()
site_results <- list()

for (h in habitats) {
  for (metric in metrics) {
    for (v in predictors) {
      key <- paste(h,metric,v,sep="::")
      obs <- observed[key]
      distribution <- null[,key]
      centre <- mean(distribution)

      p <- (1+sum(
        abs(distribution-centre) >=
          abs(obs-centre)-1e-12
      ))/(nperm+1)

      site_r <- sapply(sites,function(s) {
        observed_sites[[paste(h,s,sep="::")]][
          paste(metric,v,sep="::")
        ]
      })

      # Descriptive site bootstrap; only three sites.
      boot <- replicate(1999,mean(
        sample(site_r,length(site_r),replace=TRUE)
      ))

      interval <- quantile(boot,c(0.025,0.975))

      results[[length(results)+1]] <- data.frame(
        Habitat=h,
        Response=metric,
        Environmental_predictor=v,
        Mean_site_partial_r=unname(obs),
        Site_bootstrap_95_low=unname(interval[1]),
        Site_bootstrap_95_high=unname(interval[2]),
        Permutation_p=p,
        N_permutations=nperm,
        N_sites=length(sites)
      )

      for (s in sites) {
        obj <- datasets[[paste(h,s,sep="::")]]
        site_results[[length(site_results)+1]] <-
          data.frame(
            Habitat=h,
            Site=s,
            Response=metric,
            Environmental_predictor=v,
            Partial_r=unname(site_r[s]),
            N_samples=nrow(obj$meta),
            N_pairs=obj$n_pairs
          )
      }
    }
  }
}

results <- do.call(rbind,results)
site_results <- do.call(rbind,site_results)

# Exploratory FDR family: 3 predictors x 2 habitats
# x 2 betaNTI response definitions.
results$FDR_q <- p.adjust(
  results$Permutation_p,method="BH"
)

write.table(
  results,
  paste0(base,"Comment04_environment_assembly_permutation.tsv"),
  sep="\t",row.names=FALSE,quote=FALSE
)

write.table(
  site_results,
  paste0(base,"Comment04_environment_assembly_site_effects.tsv"),
  sep="\t",row.names=FALSE,quote=FALSE
)

cat("\n=== COMMENT 4: ENVIRONMENT AND ASSEMBLY ===\n")
print(results,row.names=FALSE,digits=4)

cat("\nSUCCESS: Analysis and supplementary tables saved\n")
cat("NOTE: Site bootstrap intervals are exploratory (3 sites).\n")
cat("NOTE: Tests assess within-condition environmental associations.\n")
