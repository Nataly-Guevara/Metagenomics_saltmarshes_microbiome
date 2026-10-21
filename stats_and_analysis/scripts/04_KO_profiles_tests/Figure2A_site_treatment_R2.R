# Descriptive distance-based variance partitioning
# Based on KO Bray-Curtis distances; base R only.

input <- "results/functional_annotation/figure2_functional_reorganisation/KO_BrayCurtis_distances.tsv"
output <- "results/functional_annotation/figure2_functional_reorganisation/Figure2A_site_treatment_R2.tsv"

d <- read.delim(input, check.names = FALSE)

analyse <- function(habitat) {

  x <- subset(d, Habitat == habitat)

  meta <- unique(rbind(
    data.frame(Sample=x$Sample1, Site=x$Site1,
               Treatment=x$Treatment1),
    data.frame(Sample=x$Sample2, Site=x$Site2,
               Treatment=x$Treatment2)
  ))

  if (anyDuplicated(meta$Sample))
    stop("Conflicting metadata")

  samples <- sort(meta$Sample)
  n <- length(samples)

  stopifnot(n == if (habitat == "BS") 11 else 9)

  meta <- meta[match(samples, meta$Sample), ]
  meta$Site <- factor(meta$Site)
  meta$Treatment <- factor(meta$Treatment)

  # Reconstruct symmetric Bray-Curtis matrix
  D <- matrix(NA_real_, n, n,
              dimnames=list(samples, samples))
  diag(D) <- 0

  for (i in seq_len(nrow(x))) {
    a <- x$Sample1[i]
    b <- x$Sample2[i]
    D[a,b] <- x$BrayCurtis[i]
    D[b,a] <- x$BrayCurtis[i]
  }

  stopifnot(!anyNA(D))
  stopifnot(isTRUE(all.equal(D, t(D))))

  # Gower-centred squared-distance matrix
  J <- diag(n) - matrix(1/n, n, n)
  G <- -0.5 * J %*% (D^2) %*% J

  totalSS <- sum(diag(G))

  # Model-explained sum of squares
  modelSS <- function(formula) {
    X <- model.matrix(formula, data=meta)
    Q <- qr.Q(qr(X))
    sum(Q * (G %*% Q))
  }

  ss_site <- modelSS(~ Site)
  ss_treat <- modelSS(~ Treatment)
  ss_full <- modelSS(~ Site + Treatment)

  site_first <- ss_site / totalSS
  treat_after_site <- (ss_full - ss_site) / totalSS

  treat_first <- ss_treat / totalSS
  site_after_treat <- (ss_full - ss_treat) / totalSS

  residual <- (totalSS - ss_full) / totalSS

  result <- data.frame(
    Habitat = habitat,
    N = n,
    Site_R2_first = site_first,
    Treatment_R2_after_Site = treat_after_site,
    Residual_R2 = residual,
    Treatment_R2_first = treat_first,
    Site_R2_after_Treatment = site_after_treat
  )

  cat("\n", habitat, " | ", n, " metagenomes\n", sep="")
  cat(sprintf("Site R2:                 %.1f%%\n", 100 * site_first))
  cat(sprintf("Treatment R2 after Site: %.1f%%\n", 100 * treat_after_site))
  cat(sprintf("Residual R2:             %.1f%%\n", 100 * residual))
  cat(sprintf("Site R2 after Treatment: %.1f%%\n", 100 * site_after_treat))

  result
}

results <- rbind(analyse("BS"), analyse("RZ"))

write.table(
  results, output,
  sep="\t", row.names=FALSE, quote=FALSE
)

cat("\nSUCCESS\n")
cat("Output:", output, "\n")
