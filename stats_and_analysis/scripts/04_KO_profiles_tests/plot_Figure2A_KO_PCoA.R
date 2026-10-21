library(ggplot2)
library(grid)

input <- "results/functional_annotation/figure2_functional_reorganisation/KO_BrayCurtis_distances.tsv"
output <- "results/functional_annotation/figure2_functional_reorganisation/Figure2A_KO_BrayCurtis_PCoA_BS_RZ.svg"

d <- read.delim(input, check.names = FALSE)

stopifnot(all(c(
  "Sample1", "Sample2", "Site1", "Site2",
  "Treatment1", "Treatment2", "Habitat", "BrayCurtis"
) %in% names(d)))

treatment_colors <- c(
  "NATURAL" = "#9BC1BC",
  "BESE" = "#FF5A5F",
  "BARE" = "#FFA630",
  "DEGRADED" = "#5D576B"
)

site_shapes <- c(
  "ACP" = 16,
  "BASS" = 17,
  "TRUG" = 15
)

make_pcoa <- function(habitat) {

  x <- subset(d, Habitat == habitat)

  metadata <- unique(rbind(
    data.frame(
      Sample = x$Sample1,
      Site = x$Site1,
      Treatment = x$Treatment1
    ),
    data.frame(
      Sample = x$Sample2,
      Site = x$Site2,
      Treatment = x$Treatment2
    )
  ))

  if (anyDuplicated(metadata$Sample)) {
    stop("Conflicting metadata for samples")
  }

  samples <- sort(metadata$Sample)
  n <- length(samples)

  expected <- if (habitat == "BS") 11 else 9
  stopifnot(n == expected)

  mat <- matrix(
    NA_real_, n, n,
    dimnames = list(samples, samples)
  )
  diag(mat) <- 0

  for (i in seq_len(nrow(x))) {
    a <- x$Sample1[i]
    b <- x$Sample2[i]
    mat[a, b] <- x$BrayCurtis[i]
    mat[b, a] <- x$BrayCurtis[i]
  }

  if (anyNA(mat)) {
    stop(paste("Missing distances in", habitat))
  }

  # Cailliez correction handles non-Euclidean Bray-Curtis distances.
  pc <- cmdscale(as.dist(mat), k = 2, eig = TRUE, add = TRUE)

  explained <- 100 * pc$eig[1:2] /
    sum(pc$eig[pc$eig > 0])

  coords <- data.frame(
    Sample = rownames(pc$points),
    PCoA1 = pc$points[, 1],
    PCoA2 = pc$points[, 2]
  )

  coords <- merge(coords, metadata, by = "Sample", sort = FALSE)

  coords$Treatment <- recode(
    coords$Treatment,
    NAT = "NATURAL",
    DEGR = "DEGRADED"
  )

  coords$Treatment <- factor(
    coords$Treatment,
    levels = names(treatment_colors)
  )

  coords$Site <- factor(
    coords$Site,
    levels = names(site_shapes)
  )

  title <- if (habitat == "BS") {
    "i. Bulk soil"
  } else {
    "ii. Rhizosphere"
  }

  p <- ggplot(
    coords,
    aes(PCoA1, PCoA2, color = Treatment, shape = Site)
  ) +
    geom_point(size = 3.8, stroke = 0.7) +
    scale_color_manual(
      values = treatment_colors,
      drop = TRUE
    ) +
    scale_shape_manual(values = site_shapes) +
    coord_equal() +
    labs(
      title = paste0(title, " (n = ", n, ")"),
      x = sprintf("PCoA1 (%.1f%%)", explained[1]),
      y = sprintf("PCoA2 (%.1f%%)", explained[2]),
      color = "Treatment",
      shape = "Site"
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 15),
      axis.title = element_text(size = 12),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.text = element_text(size = 10),
      plot.margin = margin(12, 12, 12, 12)
    )

  cat(
    habitat, ": n =", n,
    "| PCoA1 =", round(explained[1], 1), "%",
    "| PCoA2 =", round(explained[2], 1), "%\n"
  )

  p
}

p_bs <- make_pcoa("BS")
p_rz <- make_pcoa("RZ")

svg(output, width = 13, height = 6.5)

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))

print(
  p_bs,
  vp = viewport(layout.pos.row = 1, layout.pos.col = 1),
  newpage = FALSE
)

print(
  p_rz,
  vp = viewport(layout.pos.row = 1, layout.pos.col = 2),
  newpage = FALSE
)

dev.off()

cat("SUCCESS\n")
cat("Output:", output, "\n")
