# Figure 4A: MAG niche breadth versus mean relative abundance
# Uses existing MAG depth, genome breadth, niche breadth, and taxonomy

base <- "results/ecological_strategies"
abdir <- "mapping/competitive_MAGs/MAG_abundance"

niche <- read.delim(
  file.path(base, "MAG_niche_breadth_GRESP_373.tsv"),
  check.names = FALSE
)

depth <- read.delim(
  file.path(abdir, "final_373_MAGs_mean_depth_matrix.tsv"),
  check.names = FALSE
)

breadth <- read.delim(
  file.path(abdir, "final_373_MAGs_breadth_matrix.tsv"),
  check.names = FALSE
)

tax <- read.delim(
  "results/final_summaries/combined_binners/final_373_MAGs_GTDB_taxonomy.tsv",
  check.names = FALSE
)

# Check MAG and sample matching
stopifnot(
  nrow(niche) == 373,
  nrow(depth) == 373,
  nrow(breadth) == 373,
  identical(names(depth), names(breadth))
)

idx <- match(niche$MAG_ID, depth$MAG_ID)

stopifnot(
  !anyNA(idx),
  identical(depth$MAG_ID, breadth$MAG_ID)
)

D <- as.matrix(depth[idx, -1, drop = FALSE])
B <- as.matrix(breadth[idx, -1, drop = FALSE])

storage.mode(D) <- "numeric"
storage.mode(B) <- "numeric"

stopifnot(
  ncol(D) == 20,
  all(D >= 0),
  all(B >= 0 & B <= 1),
  all(colSums(D) > 0)
)

# Relative abundance within the recovered MAG catalogue,
# calculated independently for each metagenome.
relative <- sweep(D, 2, colSums(D), "/")

# Apply consistent detection thresholds
relative10 <- relative * (B >= 0.10)
relative50 <- relative * (B >= 0.50)

# Mean abundance (%) across all 20 sampled environments
mean10 <- rowMeans(relative10) * 100
mean50 <- rowMeans(relative50) * 100

stopifnot(all(mean10 > 0), all(mean50 > 0))

# Extract GTDB phylum
tax_idx <- match(niche$MAG_ID, tax$Final_MAG)
stopifnot(!anyNA(tax_idx))

classification <- tax$GTDB_classification[tax_idx]

get_phylum <- function(x) {
  hit <- regmatches(x, regexpr("p__[^;]+", x))
  if (length(hit) == 0 || is.na(hit) || hit == "")
    return("Unclassified")
  sub("^p__", "", hit)
}

phyla <- vapply(classification, get_phylum, character(1))

# Seven most represented phyla + Other
top_phyla <- names(sort(table(phyla), decreasing = TRUE))
top_phyla <- head(top_phyla, 7)

groups <- ifelse(phyla %in% top_phyla, phyla, "Other")

palette <- c(
  "#3977AA",
  "#E29A4A",
  "#5E9E60",
  "#9672B5",
  "#B65D80",
  "#46A8A4",
  "#A47C50",
  "#A8ADB8"
)

levels_group <- c(top_phyla, "Other")
colours <- setNames(
  palette[seq_along(levels_group)],
  levels_group
)

# Build plotting data
df <- data.frame(
  MAG_ID = niche$MAG_ID,
  Phylum = phyla,
  Phylum_group = groups,
  Mean_abundance_10_pct = mean10,
  Mean_abundance_50_pct = mean50,
  Levins_10 = niche$Levins_10,
  Levins_50 = niche$Levins_50
)

outdir <- file.path(base, "figure4_niche")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

write.table(
  df,
  file.path(outdir, "Figure4A_MAG_scatter_data.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Common axes for the two panels
all_abundances <- c(mean10, mean50)

xmin <- 10^floor(log10(min(all_abundances)))
xmax <- 10^ceiling(log10(max(all_abundances)))

xticks <- 10^seq(
  floor(log10(xmin)),
  ceiling(log10(xmax))
)

ymax <- ceiling(max(niche$Levins_10, niche$Levins_50))
yticks <- seq(0, ymax, by = 2)

draw_panel <- function(x, y, title) {

  plot(
    x, y,
    type = "n",
    log = "x",
    xlim = c(xmin, xmax),
    ylim = c(0, ymax),
    axes = FALSE,
    xlab = "Mean relative MAG abundance (%)",
    ylab = "Levins' distribution breadth",
    main = title,
    cex.main = 1.2,
    cex.lab = 1.1
  )

  abline(v = xticks, col = "#E8EBEF")
  abline(h = yticks, col = "#E8EBEF")

  axis(
    1, at = xticks,
    labels = format(xticks, scientific = FALSE, trim = TRUE),
    cex.axis = 0.85
  )

  axis(2, at = yticks, las = 1, cex.axis = 0.9)
  box(col = "#777777")

  # One point per MAG, coloured by taxonomic group
  for (group in levels_group) {

    i <- which(df$Phylum_group == group)

    points(
      x[i], y[i],
      pch = 21,
      bg = adjustcolor(colours[group], alpha.f = 0.70),
      col = adjustcolor(colours[group], alpha.f = 0.85),
      cex = 0.95,
      lwd = 0.35
    )
  }
}

# Export publication-quality vector graphic
svg(
  file.path(outdir, "Figure4A_MAG_niche_vs_abundance.svg"),
  width = 13,
  height = 7,
  pointsize = 12
)

layout(
  matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE),
  heights = c(5.2, 1.1)
)

par(mar = c(5, 5, 3, 1))

draw_panel(
  mean10, niche$Levins_10,
  "A | Genome breadth >= 10%"
)

draw_panel(
  mean50, niche$Levins_50,
  "B | Genome breadth >= 50%"
)

# Shared phylum legend
par(mar = c(0, 1, 0, 1))
plot.new()

legend(
  "center",
  legend = levels_group,
  pch = 21,
  pt.bg = unname(colours),
  col = unname(colours),
  ncol = 4,
  bty = "n",
  cex = 0.85,
  pt.cex = 1.2
)

dev.off()

cat("\nSUCCESS\n")
cat("MAGs plotted:", nrow(df), "\n")
cat("Phylum groups:", paste(levels_group, collapse = ", "), "\n")
cat("Output:", file.path(outdir, "Figure4A_MAG_niche_vs_abundance.svg"), "\n")
