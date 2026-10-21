#!/usr/bin/env Rscript
# Integrated 37-capability heatmap for saltmarsh restoration metagenomes.
# One figure: four BS conditions and three RZ conditions.
# Source: sample-level abundance-weighted MAG capability prevalence.
# Authoritative detection rule: 10% MAG genome breadth (change below for sensitivity).
# Only ggplot2 is required; no re-annotation or re-mapping.

suppressPackageStartupMessages(library(ggplot2))

FILTER <- "Breadth10"       # Sensitivity: "Breadth50"
SHOW_VALUES <- TRUE         # FALSE gives an uncluttered, number-free version
SHOW_MAG_COUNTS <- TRUE     # Mean DETECTED MAGs/sample, not assembly-origin counts

source_file <- paste0(
  "results/ecological_strategies/community_weighted_GRESP/",
  "community_capabilities_20_samples_long.tsv"
)
depth_file <- paste0(
  "mapping/competitive_MAGs/MAG_abundance/",
  "final_373_MAGs_mean_depth_matrix.tsv"
)
breadth_file <- paste0(
  "mapping/competitive_MAGs/MAG_abundance/",
  "final_373_MAGs_breadth_matrix.tsv"
)
outdir <- "results/functional_annotation/figure2_functional_reorganisation"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

stopifnot(file.exists(source_file), file.exists(depth_file), file.exists(breadth_file))
raw <- read.delim(source_file, check.names=FALSE, stringsAsFactors=FALSE)
required <- c("Sample", "Site", "Treatment", "Habitat", "Filter",
              "GRESP", "Capability", "Weighted_prevalence_pct")
stopifnot(all(required %in% colnames(raw)))
raw <- raw[raw$Filter == FILTER, , drop=FALSE]
stopifnot(nrow(raw) == 20L * 37L, length(unique(raw$Sample)) == 20L)
stopifnot(!anyDuplicated(raw[, c("Sample", "GRESP", "Capability")]))
stopifnot(all(is.finite(raw$Weighted_prevalence_pct)),
          all(raw$Weighted_prevalence_pct >= 0 & raw$Weighted_prevalence_pct <= 100))

# Each name is paired with the original GRESP group; e.g., DNRA occurs
# twice, but the redox and nitrogen entries have DIFFERENT definitions.
groups <- list(
  Osmotic = c(
    Capability_Trehalose = "Trehalose-associated potential",
    Capability_Compatible_solute_uptake = "Compatible-solute uptake",
    Capability_Ectoine_biosynthesis = "Ectoine biosynthesis",
    Capability_Mechanosensitive_protection = "Mechanosensitive protection",
    Capability_Ion_homeostasis = "Broad ion homeostasis"
  ),
  Redox = c(
    Capability_Aerobic_respiration = "Aerobic respiration",
    Capability_LowO2_respiration = "Low-O2 respiration",
    Capability_Nitrate_respiration = "Nitrate respiration",
    Capability_Nitrite_to_NO = "Nitrite to NO",
    Capability_DNRA = "Nitrite-to-ammonium potential (NirBD/NrfAH)",
    Capability_NO_reduction = "NO reduction",
    Capability_N2O_reduction = "N2O reduction",
    Capability_Dissimilatory_sulfate_reduction_core = "Sat-AprAB-DsrAB marker core*",
    Capability_Reduced_sulfur_oxidation = "Reduced-sulfur oxidation",
    Capability_Sulfite_oxidation = "Sulfite oxidation"
  ),
  Nitrogen = c(
    Capability_Nitrogen_fixation = "Nitrogen fixation",
    Capability_Ammonium_uptake = "Ammonium uptake",
    Capability_Assimilatory_nitrate_reduction = "Assimilatory nitrate reduction",
    Capability_Urea_utilisation = "Urea utilisation",
    Capability_DNRA = "Complete nitrate-to-ammonium potential",
    Capability_Complete_denitrification = "Complete denitrification",
    Capability_Nitrification = "Nitrification",
    Capability_Anammox = "Anammox"
  ),
  Carbon_P = c(
    Capability_Carbon_fixation = "Carbon fixation",
    Capability_C1_utilisation = "C1 utilisation",
    Capability_High_affinity_phosphate_uptake = "High-affinity phosphate uptake",
    Capability_Organic_P_mobilisation = "Organic-P mobilisation",
    Capability_Phosphonate_acquisition = "Phosphonate acquisition",
    Capability_PolyP_turnover = "Polyphosphate turnover"
  ),
  Iron = c(
    Ferrous_iron_uptake = "Ferrous-iron uptake",
    Direct_ferric_iron_uptake = "Direct ferric-iron uptake",
    Siderophore_biosynthesis = "Siderophore biosynthesis (candidate)",
    Iron_complex_uptake = "Iron-complex uptake"
  ),
  Colonisation = c(
    Chemotaxis = "Chemotaxis",
    Motility = "Motility",
    Attachment_biofilm_matrix = "Attachment / biofilm-matrix potential",
    Direct_plant_interaction = "Curated plant-interaction markers"
  )
)

spec <- do.call(rbind, lapply(names(groups), function(group) {
  data.frame(GRESP=group, Capability=names(groups[[group]]),
             Display=unname(groups[[group]]), stringsAsFactors=FALSE)
}))
rownames(spec) <- NULL
stopifnot(nrow(spec) == 37, !anyDuplicated(spec[,c("GRESP","Capability")]),
          !anyDuplicated(spec$Display))
observed_keys <- unique(paste(raw$GRESP, raw$Capability, sep="::"))
expected_keys <- paste(spec$GRESP, spec$Capability, sep="::")
stopifnot(setequal(observed_keys, expected_keys))

columns <- data.frame(
  Habitat=c("BS", "BS", "BS", "BS", "RZ", "RZ", "RZ"),
  Treatment=c("NAT", "BESE", "BARE", "DEGR", "NAT", "BESE", "BARE"),
  Col=1:7,
  Short=c("NAT", "BESE", "BARE", "DEG", "NAT", "BESE", "BARE"),
  stringsAsFactors=FALSE
)

# One metagenome per site x treatment x habitat: calculate mean across sites.
averages <- aggregate(
  Weighted_prevalence_pct ~ Habitat + Treatment + GRESP + Capability,
  data=raw, FUN=mean
)
names(averages)[names(averages)=="Weighted_prevalence_pct"] <- "Mean_pct"

# Complete 37 x 7 plotting grid.  No DEGRADED RZ is fabricated.
plotdat <- merge(spec, columns, by=NULL)
plotdat <- merge(plotdat, averages,
                 by=c("GRESP", "Capability", "Habitat", "Treatment"),
                 all.x=TRUE, sort=FALSE)
stopifnot(nrow(plotdat)==37L*7L, !anyNA(plotdat$Mean_pct))
plotdat$GRESP <- factor(plotdat$GRESP, levels=names(groups))
plotdat$Display <- factor(plotdat$Display, levels=rev(spec$Display))
plotdat$Cell_text <- sprintf("%.1f", plotdat$Mean_pct)
plotdat$Text_tone <- ifelse(plotdat$Mean_pct > 58, "light", "dark")

# Detected MAG counts: unlike MAGs first assembled from a particular source,
# these are the MAGs with reads recruited in the actual metagenome.
depth <- read.delim(depth_file, check.names=FALSE, stringsAsFactors=FALSE)
breadth <- read.delim(breadth_file, check.names=FALSE, stringsAsFactors=FALSE)
stopifnot(nrow(depth)==373L, nrow(breadth)==373L,
          setequal(depth$MAG_ID, breadth$MAG_ID))
breadth <- breadth[match(depth$MAG_ID, breadth$MAG_ID),,drop=FALSE]
sample_cols <- setdiff(names(depth), "MAG_ID")
stopifnot(setequal(sample_cols, unique(raw$Sample)))
min_breadth <- switch(FILTER, Breadth10=0.10, Breadth50=0.50,
                      stop("Unrecognised breadth filter"))
count_by_sample <- setNames(vapply(sample_cols, function(s) {
  sum(depth[[s]] > 0 & breadth[[s]] >= min_breadth)
}, integer(1)), sample_cols)
meta <- unique(raw[,c("Sample","Site","Treatment","Habitat")])
stopifnot(nrow(meta)==20L)
meta$Detected_MAGs <- unname(count_by_sample[meta$Sample])
stopifnot(!anyNA(meta$Detected_MAGs))

column_labels <- vapply(seq_len(nrow(columns)), function(i) {
  x <- columns[i,]
  sub <- meta[meta$Habitat==x$Habitat & meta$Treatment==x$Treatment,,drop=FALSE]
  stopifnot(nrow(sub) %in% c(2L,3L))
  if (SHOW_MAG_COUNTS) {
    sprintf("%s\nn=%d | MAGs~%d", x$Short, nrow(sub),
            round(mean(sub$Detected_MAGs)))
  } else {
    sprintf("%s\nn=%d", x$Short, nrow(sub))
  }
}, character(1))

labels <- c(
  Osmotic="Osmotic",
  Redox="Redox / sulfur",
  Nitrogen="Nitrogen",
  Carbon_P="Carbon / P",
  Iron="Iron",
  Colonisation="Colonisation"
)

p <- ggplot(plotdat, aes(x=Col, y=Display)) +
  geom_tile(aes(fill=Mean_pct), width=0.98, height=0.98,
            color="white", linewidth=0.32) +
  geom_vline(xintercept=4.5, linewidth=2.1, colour="white") +
  scale_fill_gradientn(
    colours=c("#FFF8CE", "#CBE6CE", "#62B7C7", "#68339B"),
    values=c(0, 0.18, 0.49, 1),
    limits=c(0,100), breaks=seq(0,100,20),
    name="Mapped MAG abundance\nencoding capability (%)"
  ) +
  scale_x_continuous(
    breaks=1:7, labels=column_labels,
    limits=c(0.5,7.5), expand=expansion(mult=0, add=0.01),
    position="top"
  ) +
  scale_y_discrete(drop=TRUE, expand=expansion(add=0.04)) +
  facet_grid(GRESP ~ ., scales="free_y", space="free_y", switch="y",
             labeller=labeller(GRESP=as_labeller(labels))) +
  labs(
    title="Genomic survival and resource-acquisition potential",
    subtitle=paste0("Mean abundance-weighted capability prevalence across sites | ", FILTER,
                    " genome-breadth filter"),
    x=NULL, y=NULL,
    caption=paste(
      "Each cell is the unweighted mean of site-level metagenomic prevalence values (one pooled metagenome per site x condition x habitat).",
      "n = site-level pooled metagenomes; MAGs = mean number detected per metagenome (not number assembled from the condition).",
      "* DsrAB directionality not resolved. Curated genomic potential does not establish expression or activity.",
      sep="\n"
    )
  ) +
  theme_minimal(base_size=11) +
  theme(
    plot.title=element_text(face="bold", size=15, margin=margin(b=7)),
    plot.subtitle=element_text(size=10.5, colour="#4D5966", margin=margin(b=8)),
    plot.caption=element_text(hjust=0, size=8.3, colour="#5D6571", lineheight=1.1,
                              margin=margin(t=13)),
    axis.text.x.top=element_text(size=9.5, face="bold", colour="#1F2937",
                                 lineheight=1.12, margin=margin(b=10)),
    axis.text.y=element_text(size=10.0, colour="#1F2937", margin=margin(r=5)),
    axis.ticks=element_blank(),
    panel.grid=element_blank(),
    panel.spacing.y=grid::unit(0.08,"cm"),
    strip.placement="outside",
    strip.background=element_rect(fill="#EFF2F5", colour=NA),
    strip.text.y.left=element_text(angle=0, size=9.5, face="bold", colour="#324359",
                                   margin=margin(l=7,r=7)),
    legend.title=element_text(size=10),
    legend.text=element_text(size=9),
    plot.margin=margin(t=20,r=20,b=15,l=17)
  )

if (SHOW_VALUES) {
  p <- p + geom_text(
    aes(label=Cell_text, colour=Text_tone),
    size=2.9, fontface="plain", show.legend=FALSE
  ) + scale_colour_manual(values=c(dark="#18232F",light="#FFFFFF"),guide="none")
}

# Add genuine group headings aligned with the 7 equal-width heatmap columns.
# This changes only presentation. ggplot2 already depends on gtable.
make_figure <- function(plot_object) {
  gt <- ggplot2::ggplotGrob(plot_object)
  ii <- grep("^panel", gt$layout$name)
  stopifnot(length(ii) == 6L)
  left <- min(gt$layout$l[ii]); right <- max(gt$layout$r[ii])
  top <- min(gt$layout$t[ii])
  # One row directly above the plot panels, below the column labels.
  gt <- gtable::gtable_add_rows(gt, heights=grid::unit(0.35, "in"), pos=top-1L)
  grp <- grid::grobTree(
    grid::rectGrob(x=2/7, width=4/7, y=0.5, height=0.87,
                   gp=grid::gpar(fill="#E6F0F0", col="white", lwd=2)),
    grid::rectGrob(x=5.5/7, width=3/7, y=0.5, height=0.87,
                   gp=grid::gpar(fill="#E9ECF5", col="white", lwd=2)),
    grid::textGrob("BULK SOIL  (BS)", x=2/7, y=0.5,
                   gp=grid::gpar(fontsize=12, fontface="bold", col="#22494B")),
    grid::textGrob("RHIZOSPHERE  (RZ)", x=5.5/7, y=0.5,
                   gp=grid::gpar(fontsize=12, fontface="bold", col="#3F486B"))
  )
  gt <- gtable::gtable_add_grob(gt, grp,
                                t=top, b=top, l=left, r=right,
                                name="soil-habitat-superheaders", clip="off")
  gt
}

combined <- make_figure(p)
prefix <- file.path(outdir, paste0("Integrated_37_capability_heatmap_",FILTER,"_v2"))
# Publication-ready vector output; show exactly the same values as v1.
grDevices::svg(paste0(prefix,".svg"), width=14.0,height=13.0,bg="white")
grid::grid.newpage(); grid::grid.draw(combined)
grDevices::dev.off()
grDevices::pdf(paste0(prefix,".pdf"),width=14.0,height=13.0,
               onefile=FALSE,useDingbats=FALSE)
grid::grid.newpage(); grid::grid.draw(combined)
grDevices::dev.off()
if (capabilities("cairo")) {
  grDevices::png(paste0(prefix,".png"),width=14,height=13,
                 units="in",res=220,type="cairo",bg="white")
  grid::grid.newpage(); grid::grid.draw(combined)
  grDevices::dev.off()
}

# Auditable plotting data: 259 condition x habitat x capability means.
write.table(plotdat[,c("GRESP","Capability","Display","Habitat", "Treatment", "Mean_pct")],
            paste0(prefix,"_data.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
write.table(meta[,c("Sample","Site","Treatment","Habitat","Detected_MAGs")],
            paste0(prefix,"_sample_counts.tsv"),sep="\t",row.names=FALSE,quote=FALSE)

cat("SUCCESS: Refined integrated BS/RZ heatmap V2\n")
cat("Filter:",FILTER,"\n")
cat("Capabilities:",length(unique(paste(plotdat$GRESP,plotdat$Capability))),"\n")
cat("Columns:",nrow(columns),"\n")
cat("Metagenomes by column:\n")
for (i in seq_len(nrow(columns))) {
  cat(columns$Habitat[i],columns$Treatment[i],": ",
      gsub("\n", ", ",column_labels[i]),"\n",sep="")
}
cat("SVG:",paste0(prefix,".svg"),"\n")
cat("PDF:",paste0(prefix,".pdf"),"\n")
