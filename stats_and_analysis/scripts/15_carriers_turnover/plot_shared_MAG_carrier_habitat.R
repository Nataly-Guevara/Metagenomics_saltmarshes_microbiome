suppressPackageStartupMessages(library(ggplot2))

folder <- paste0(
  "results/ecological_strategies/",
  "functional_carrier_taxonomy/carrier_turnover/"
)

d <- read.delim(
  paste0(folder, "Whole_MAG_vs_capability_paired_habitats.tsv"),
  check.names = FALSE
)
d <- subset(d, Filter == "Breadth10")

stopifnot(nrow(d) == 7)

long <- rbind(
  data.frame(d[, c("Site", "Contrast")],
             Metric = "Whole MAGs",
             Habitat = "BS",
             Pct = d$BS_whole_MAG_pct),
  data.frame(d[, c("Site", "Contrast")],
             Metric = "Whole MAGs",
             Habitat = "RZ",
             Pct = d$RZ_whole_MAG_pct),
  data.frame(d[, c("Site", "Contrast")],
             Metric = "Capability-weighted",
             Habitat = "BS",
             Pct = d$BS_capability_pct),
  data.frame(d[, c("Site", "Contrast")],
             Metric = "Capability-weighted",
             Habitat = "RZ",
             Pct = d$RZ_capability_pct)
)

long$Site <- factor(long$Site,
                    levels = c("ACP", "BASS", "TRUG"))
long$Contrast <- factor(
  long$Contrast,
  levels = c("NAT_vs_BESE", "NAT_vs_BARE", "BESE_vs_BARE"),
  labels = c("Natural–BESE", "Natural–BARE", "BESE–BARE")
)
long$Habitat <- factor(
  long$Habitat,
  levels = c("BS", "RZ"),
  labels = c("Bulk soil", "Rhizosphere")
)
long$Metric <- factor(
  long$Metric,
  levels = c("Whole MAGs", "Capability-weighted")
)

p <- ggplot(
  long,
  aes(x = Habitat, y = Pct, group = Site, colour = Site)
) +
  geom_line(linewidth = 0.9, alpha = 0.85) +
  geom_point(size = 2.9) +
  facet_grid(Contrast ~ Metric) +
  scale_colour_manual(values = c(
    ACP = "#346B80",
    BASS = "#BB7D36",
    TRUG = "#8770A7"
  )) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    title = "Shared MAG populations and functional carriers",
    subtitle = "Site-matched habitat comparisons | 10% genome breadth",
    x = NULL,
    y = "Abundance contributed by shared MAGs (%)",
    colour = "Site",
    caption = paste(
      "One metagenome per site × condition × habitat.",
      "TRUG BARE bulk soil was unavailable.",
      "Capability-weighted values count MAG contributions across traits."
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "#F1F3F4"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom",
    plot.caption = element_text(size = 8, hjust = 0),
    axis.text.x = element_text(size = 9),
    plot.margin = margin(12, 12, 12, 12)
  )

pdf_file <- paste0(folder, "Shared_MAG_carrier_habitats.pdf")
svg_file <- paste0(folder, "Shared_MAG_carrier_habitats.svg")

pdf(pdf_file, width = 9.5, height = 8.5,
    useDingbats = FALSE)
print(p)
dev.off()

svg(svg_file, width = 9.5, height = 8.5)
print(p)
dev.off()

cat("\nSUCCESS: Site-matched carrier-sharing figure\n")
cat("PDF:", pdf_file, "\n")
cat("SVG:", svg_file, "\n")
