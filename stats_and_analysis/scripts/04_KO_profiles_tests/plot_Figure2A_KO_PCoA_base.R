# Figure 2A: KO-based Bray-Curtis PCoA
# Requires only base R.

input <- "results/functional_annotation/figure2_functional_reorganisation/KO_BrayCurtis_distances.tsv"
output <- "results/functional_annotation/figure2_functional_reorganisation/Figure2A_KO_PCoA.svg"

d <- read.delim(input, check.names = FALSE)

treatment_colors <- c(
  NATURAL = "#9BC1BC",
  BESE = "#FF5A5F",
  BARE = "#FFA630",
  DEGRADED = "#5D576B"
)

site_shapes <- c(
  ACP = 21,   # Circle
  BASS = 22,  # Square
  TRUG = 24   # Triangle
)

make_pcoa <- function(habitat) {

  x <- d[d$Habitat == habitat, ]

  meta <- unique(rbind(
    data.frame(Sample=x$Sample1, Site=x$Site1,
               Treatment=x$Treatment1),
    data.frame(Sample=x$Sample2, Site=x$Site2,
               Treatment=x$Treatment2)
  ))

  stopifnot(!anyDuplicated(meta$Sample))

  samples <- sort(unique(c(x$Sample1, x$Sample2)))
  n <- length(samples)

  stopifnot(n == if (habitat == "BS") 11 else 9)

  mat <- matrix(NA_real_, n, n,
                dimnames=list(samples, samples))
  diag(mat) <- 0

  for (i in seq_len(nrow(x))) {
    a <- x$Sample1[i]
    b <- x$Sample2[i]
    mat[a,b] <- x$BrayCurtis[i]
    mat[b,a] <- x$BrayCurtis[i]
  }

  stopifnot(!anyNA(mat))

  # Cailliez correction for non-Euclidean distances
  pc <- cmdscale(as.dist(mat), k=2, eig=TRUE, add=TRUE)

  explained <- 100 * pc$eig[1:2] /
    sum(pc$eig[pc$eig > 0])

  coords <- pc$points
  meta <- meta[match(rownames(coords), meta$Sample), ]

  treatment <- as.character(meta$Treatment)
  treatment[treatment == "NAT"] <- "NATURAL"
  treatment[treatment == "DEGR"] <- "DEGRADED"

  stopifnot(
    all(treatment %in% names(treatment_colors)),
    all(meta$Site %in% names(site_shapes))
  )

  title <- if (habitat == "BS") {
    "i. Bulk soil"
  } else {
    "ii. Rhizosphere"
  }

  plot(
    coords[,1], coords[,2],
    type="n", asp=1,
    main=title,
    xlab=sprintf("PCoA1 (%.1f%%)", explained[1]),
    ylab=sprintf("PCoA2 (%.1f%%)", explained[2]),
    cex.main=1.35,
    cex.lab=1.1,
    las=1
  )

  abline(h=0, v=0, col="#E5E5E5", lty=3)

  points(
    coords[,1], coords[,2],
    pch=unname(site_shapes[meta$Site]),
    bg=unname(treatment_colors[treatment]),
    col="#333333",
    cex=1.9,
    lwd=0.8
  )

  cat(
    habitat,
    "| Samples:", n,
    "| PCoA1:", round(explained[1],1), "%",
    "| PCoA2:", round(explained[2],1), "%\n"
  )
}

svg(output, width=12.5, height=6.3, pointsize=11)

layout(
  matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
  heights=c(5,0.85)
)

par(mar=c(5,5,3.5,1.5))
make_pcoa("BS")
make_pcoa("RZ")

# Shared legend
par(mar=c(0,0,0,0))
plot.new()
plot.window(xlim=c(0,12.5), ylim=c(0,1))

text(0.15,0.55,"Treatment:",font=2,adj=0)

tx <- c(1.5,3.25,4.75,6.05)
tn <- names(treatment_colors)

points(tx,rep(0.55,4),pch=15,
       col=treatment_colors,cex=1.4)
text(tx+0.18,rep(0.55,4),tn,adj=0,cex=0.85)

text(8.0,0.55,"Site:",font=2,adj=0)

sx <- c(8.75,10.05,11.35)
sn <- names(site_shapes)

points(sx,rep(0.55,3),pch=site_shapes,
       bg="white",col="#333333",cex=1.2)
text(sx+0.16,rep(0.55,3),sn,adj=0,cex=0.85)

dev.off()

cat("SUCCESS\n")
cat("Output:", output, "\n")
