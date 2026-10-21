# Figure 4A: Community-weighted MAG niche breadth by restoration treatment
# Base R only

base <- "results/ecological_strategies"
abdir <- "mapping/competitive_MAGs/MAG_abundance"
outdir <- file.path(base, "figure4_niche")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

niche <- read.delim(file.path(base, "MAG_niche_breadth_GRESP_373.tsv"))
depth <- read.delim(
  file.path(abdir, "final_373_MAGs_mean_depth_matrix.tsv"),
  check.names=FALSE
)
breadth <- read.delim(
  file.path(abdir, "final_373_MAGs_breadth_matrix.tsv"),
  check.names=FALSE
)

stopifnot(
  nrow(niche)==373,
  nrow(depth)==373,
  identical(names(depth), names(breadth)),
  identical(depth$MAG_ID, breadth$MAG_ID),
  !anyDuplicated(niche$MAG_ID)
)

idx <- match(depth$MAG_ID, niche$MAG_ID)
stopifnot(!anyNA(idx))

D <- as.matrix(depth[-1])
B <- as.matrix(breadth[-1])
storage.mode(D) <- "numeric"
storage.mode(B) <- "numeric"

samples <- colnames(D)
stopifnot(length(samples)==20)

records <- list()

for (threshold in c(10, 50)) {

  levins <- niche[[paste0("Levins_", threshold)]][idx]

  for (j in seq_along(samples)) {
    sample <- samples[j]
    parts <- strsplit(sample, "_", fixed=TRUE)[[1]]

    weights <- D[, j]
    weights[B[, j] < threshold/100] <- 0

    stopifnot(sum(weights) > 0)

    records[[length(records)+1]] <- data.frame(
      Sample=sample,
      Site=parts[2],
      Treatment=parts[3],
      Habitat=parts[4],
      Breadth_threshold=threshold,
      Weighted_Levins=weighted.mean(levins, weights),
      Detected_MAGs=sum(weights > 0)
    )
  }
}

df <- do.call(rbind, records)
stopifnot(nrow(df)==40)

write.table(
  df,
  file.path(outdir, "Figure4A_treatment_niche_data.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)

treatment_colors <- c(
  NAT="#9BC1BC",
  BESE="#FF5A5F",
  BARE="#FFA630",
  DEGR="#5D576B"
)

site_shapes <- c(ACP=21, BASS=22, TRUG=24)
site_offset <- c(ACP=-0.065, BASS=0, TRUG=0.065)

for (threshold in c(10,50)) {

  dat <- subset(df, Breadth_threshold==threshold)

  output <- file.path(
    outdir,
    paste0("Figure4A_restoration_niche_breadth", threshold, ".svg")
  )

  svg(output, width=11.5, height=5.8, pointsize=12)

  par(
    mfrow=c(1,2),
    mar=c(4.5,5,3,1),
    oma=c(2.5,0,0,0)
  )

  ylims <- range(dat$Weighted_Levins)
  ylims <- ylims + c(-0.7,0.7)

  for (hab in c("BS","RZ")) {

    x <- subset(dat, Habitat==hab)

    treatments <- if (hab=="BS")
      c("NAT","BESE","BARE","DEGR")
    else
      c("NAT","BESE","BARE")

    title <- if (hab=="BS")
      "A | Bulk soil"
    else
      "B | Rhizosphere"

    plot(
      NA,
      xlim=c(0.5,length(treatments)+0.5),
      ylim=ylims,
      xaxt="n",
      xlab="",
      ylab="Abundance-weighted MAG

cat > scripts/plot_Figure4A_restoration_niche.R <<'RS'
# Figure 4A: Community-weighted MAG niche breadth by restoration treatment
# Base R only

base <- "results/ecological_strategies"
abdir <- "mapping/competitive_MAGs/MAG_abundance"
outdir <- file.path(base, "figure4_niche")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

niche <- read.delim(file.path(base, "MAG_niche_breadth_GRESP_373.tsv"))
depth <- read.delim(
  file.path(abdir, "final_373_MAGs_mean_depth_matrix.tsv"),
  check.names=FALSE
)
breadth <- read.delim(
  file.path(abdir, "final_373_MAGs_breadth_matrix.tsv"),
  check.names=FALSE
)

stopifnot(
  nrow(niche)==373,
  nrow(depth)==373,
  identical(names(depth), names(breadth)),
  identical(depth$MAG_ID, breadth$MAG_ID),
  !anyDuplicated(niche$MAG_ID)
)

idx <- match(depth$MAG_ID, niche$MAG_ID)
stopifnot(!anyNA(idx))

D <- as.matrix(depth[-1])
B <- as.matrix(breadth[-1])
storage.mode(D) <- "numeric"
storage.mode(B) <- "numeric"

samples <- colnames(D)
stopifnot(length(samples)==20)

records <- list()

for (threshold in c(10, 50)) {

  levins <- niche[[paste0("Levins_", threshold)]][idx]

  for (j in seq_along(samples)) {
    sample <- samples[j]
    parts <- strsplit(sample, "_", fixed=TRUE)[[1]]

    weights <- D[, j]
    weights[B[, j] < threshold/100] <- 0

    stopifnot(sum(weights) > 0)

    records[[length(records)+1]] <- data.frame(
      Sample=sample,
      Site=parts[2],
      Treatment=parts[3],
      Habitat=parts[4],
      Breadth_threshold=threshold,
      Weighted_Levins=weighted.mean(levins, weights),
      Detected_MAGs=sum(weights > 0)
    )
  }
}

df <- do.call(rbind, records)
stopifnot(nrow(df)==40)

write.table(
  df,
  file.path(outdir, "Figure4A_treatment_niche_data.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)

treatment_colors <- c(
  NAT="#9BC1BC",
  BESE="#FF5A5F",
  BARE="#FFA630",
  DEGR="#5D576B"
)

site_shapes <- c(ACP=21, BASS=22, TRUG=24)
site_offset <- c(ACP=-0.065, BASS=0, TRUG=0.065)

for (threshold in c(10,50)) {

  dat <- subset(df, Breadth_threshold==threshold)

  output <- file.path(
    outdir,
    paste0("Figure4A_restoration_niche_breadth", threshold, ".svg")
  )

  svg(output, width=11.5, height=5.8, pointsize=12)

  par(
    mfrow=c(1,2),
    mar=c(4.5,5,3,1),
    oma=c(2.5,0,0,0)
  )

  ylims <- range(dat$Weighted_Levins)
  ylims <- ylims + c(-0.7,0.7)

  for (hab in c("BS","RZ")) {

    x <- subset(dat, Habitat==hab)

    treatments <- if (hab=="BS")
      c("NAT","BESE","BARE","DEGR")
    else
      c("NAT","BESE","BARE")

    title <- if (hab=="BS")
      "A | Bulk soil"
    else
      "B | Rhizosphere"

    plot(
      NA,
      xlim=c(0.5,length(treatments)+0.5),
      ylim=ylims,
      xaxt="n",
      xlab="",
      ylab="Abundance-weighted MAG niche breadth",
      main=title,
      bty="l",
      las=1
    )

    abline(h=pretty(ylims),col="#E8E8E8",lty=3)

    axis(
      1,
      at=seq_along(treatments),
      labels=c(
        NAT="NATURAL", BESE="BESE",
        BARE="BARE", DEGR="DEGRADED"
      )[treatments],
      tick=FALSE
    )

    # Connect BESE and BARE within each site
    for (site in names(site_shapes)) {
      pair <- x[
        x$Site==site &
        x$Treatment %in% c("BESE","BARE"),
      ]

      if (nrow(pair)==2) {
        xx <- match(pair$Treatment,treatments) +
              site_offset[[site]]

        lines(
          xx,
          pair$Weighted_Levins,
          col="#9BA1A8",
          lwd=1.1,
          lty=2
        )
      }
    }

    # One point per pooled metagenome
    for (i in seq_len(nrow(x))) {

      points(
        match(x$Treatment[i],treatments) +
          site_offset[[x$Site[i]],
        x$Weighted_Levins[i],
        pch=site_shapes[[x$Site[i]]],
        bg=treatment_colors[[x$Treatment[i]]],
        col="#30343B",
        cex=1.7,
        lwd=0.8
      )
    }
  }

  mtext(
    "Circle: ACP     Square: BASS     Triangle: TRUG     Dashed lines: site-matched BESE–BARE",
    side=1,outer=TRUE,line=0.6,cex=0.83,col="#526071"
  )

  dev.off()
  cat("Saved:",output,"\n")
}

cat("\nSUCCESS\n")
cat("Metagenomes:",length(samples),"\n")

cat("\nPrimary analysis: 10% genome breadth\n")
print(
  df[df$Breadth_threshold==10,
     c("Site","Treatment","Habitat","Weighted_Levins")],
  row.names=FALSE
)
