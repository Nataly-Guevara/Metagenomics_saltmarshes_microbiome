input <- "results/functional_annotation/figure2_functional_reorganisation/KO_same_site_distance_to_NAT.tsv"
output <- "results/functional_annotation/figure2_functional_reorganisation/Figure2C_distance_to_Natural.svg"

d <- read.delim(input, check.names = FALSE)

stopifnot(all(c("Site", "Habitat", "Treatment",
                "BrayCurtis_to_NAT") %in% names(d)))

stopifnot(!anyDuplicated(
  paste(d$Site, d$Habitat, d$Treatment)
))

colours <- c(
  BESE = "#FF5A5F",
  BARE = "#FFA630",
  DEGR = "#5D576B"
)

shapes <- c(
  ACP = 21,
  BASS = 22,
  TRUG = 24
)

draw_panel <- function(habitat, groups, ylim, title) {

  x <- d[d$Habitat == habitat, ]

  stopifnot(all(x$Treatment %in% groups))

  plot(
    NA,
    xlim = c(0.65, length(groups) + 0.35),
    ylim = ylim,
    xaxt = "n",
    xlab = "",
    ylab = "Bray-Curtis distance to Natural",
    main = title,
    las = 1,
    bty = "l"
  )

  abline(h = pretty(ylim), col = "#E8E8E8", lty = 3)

  axis(
    1,
    at = seq_along(groups),
    labels = ifelse(groups == "DEGR", "DEGRADED", groups),
    tick = FALSE,
    line = 0.5
  )

  for (site in names(shapes)) {

    sub <- x[x$Site == site, ]

    idx <- match(sub$Treatment, groups)
    ord <- order(idx)

    idx <- idx[ord]
    vals <- sub$BrayCurtis_to_NAT[ord]
    trt <- sub$Treatment[ord]

    if (length(idx) >= 2) {
      lines(
        idx, vals,
        col = "#A0A4AA",
        lty = 2,
        lwd = 1.2
      )
    }

    points(
      idx, vals,
      pch = shapes[[site]],
      bg = colours[trt],
      col = "#30343B",
      cex = 2.0,
      lwd = 0.9
    )
  }
}

svg(output, width = 11.5, height = 5.7, pointsize = 12)

layout(
  matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE),
  heights = c(4.5, 0.9)
)

par(mar = c(4.5, 4.8, 3, 1))
draw_panel(
  "BS",
  c("BESE", "BARE", "DEGR"),
  c(0.16, 0.33),
  "i. Bulk soil"
)

par(mar = c(4.5, 4.8, 3, 1))
draw_panel(
  "RZ",
  c("BESE", "BARE"),
  c(0.06, 0.13),
  "ii. Rhizosphere"
)

par(mar = c(0, 0, 0, 0))
plot.new()

legend(
  "top",
  legend = names(shapes),
  pch = shapes,
  pt.bg = "white",
  col = "#30343B",
  horiz = TRUE,
  bty = "n",
  cex = 0.95
)

text(
  0.5, 0.15,
  "Dashed lines connect site-matched observations; lower values indicate greater similarity to Natural. Axes differ.",
  cex = 0.75,
  col = "#526071"
)

dev.off()

cat("SUCCESS\n")
cat("Output:", output, "\n")
