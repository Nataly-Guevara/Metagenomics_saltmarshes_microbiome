#!/usr/bin/env Rscript
# Figure: habitat-associated shared MAG abundance and capability enrichment
# Source: 999-permutation intact-profile results; NO additional analysis.
suppressPackageStartupMessages(library(ggplot2))
library(grid)

base <- "results/ecological_strategies/functional_carrier_taxonomy/carrier_turnover"
input <- file.path(base, "Shared_MAG_capability_count_null_pair_results.tsv")
stopifnot(file.exists(input))
d <- read.delim(input, check.names=FALSE, stringsAsFactors=FALSE)
d <- d[d$Filter == "Breadth10", , drop=FALSE]
stopifnot(nrow(d)==16, all(c("Pair_W_pct","Pair_C_pct","Pair_C_minus_W_pp",
  "Null_pair_2.5pct_pp","Null_pair_97.5pct_pp","Null_pair_expected_pp",
  "BH_q_within_breadth") %in% names(d)))

# Match BS and RZ within each site. TRUG BARE BS was not sampled;
# omit the unpaired TRUG RZ points for Natural-BARE and BESE-BARE.
d <- d[d$Contrast == "NAT_vs_BESE" | d$Site %in% c("ACP","BASS"),,drop=FALSE]
stopifnot(nrow(d)==14)
contrast_levels <- c("NAT_vs_BESE","NAT_vs_BARE","BESE_vs_BARE")
contrast_labels <- c("Natural – BESE","Natural – BARE","BESE – BARE")
d$Contrast <- factor(d$Contrast, levels=contrast_levels, labels=contrast_labels)
d$Habitat <- factor(d$Habitat, levels=c("BS","RZ"),
                    labels=c("Bulk soil","Rhizosphere"))
d$Site <- factor(d$Site, levels=c("ACP","BASS","TRUG"))
d$key <- interaction(d$Site, d$Habitat, drop=TRUE)
stopifnot(all(table(d$Contrast) == c(6,4,4)))

site_palette <- c(ACP="#39758E",BASS="#B37A35",TRUG="#8763A5")
base_theme <- theme_bw(base_size=11) +
  theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
        strip.background=element_rect(fill="#F1F3F5", colour="#DADDE1"),
        strip.text=element_text(face="bold",size=10),
        legend.position="bottom", legend.title=element_text(size=9),
        axis.text=element_text(colour="#26282C"),
        plot.title=element_text(face="bold",size=13),
        plot.subtitle=element_text(size=9,colour="#505866"),
        plot.margin=margin(7,10,7,7))

# Plot A: show the primary biologically meaningful result (population overlap).
# Each site is a paired comparison, with no invented matching where BS unavailable.
a <- ggplot(d, aes(x=Habitat, y=Pair_W_pct, group=Site,colour=Site)) +
  geom_line(linewidth=0.8, alpha=0.85) +
  geom_point(size=2.6) +
  facet_grid(Contrast ~ ., drop=FALSE) +
  scale_colour_manual(values=site_palette) +
  scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100),
                     expand=expansion(mult=c(0.02,0.03))) +
  labs(title="A  Population sharing",subtitle="Site-matched, 10% genome breadth",
       x=NULL,y="MAG abundance from shared populations (%)", colour="Site") +
  base_theme

# Plot B: observed capability enrichment and 95% permutation range.
# Display each site-habitat independently; 16 comparisons, no pseudo-replication.
null_limit <- max(abs(c(d$Pair_C_minus_W_pp, d$Null_pair_2.5pct_pp,
                        d$Null_pair_97.5pct_pp)),na.rm=TRUE)
null_limit <- ceiling(null_limit + 0.5)
d$SiteHabitat <- factor(paste(as.character(d$Site),
  ifelse(d$Habitat=="Bulk soil","BS","RZ"),sep=" · "),
  levels=as.vector(t(outer(c("ACP","BASS","TRUG"),c("BS","RZ"),paste,sep=" · "))))
b <- ggplot(d,aes(x=SiteHabitat,y=Pair_C_minus_W_pp)) +
  geom_hline(yintercept=0, linetype="dashed", colour="#9299A3",linewidth=0.5) +
  geom_linerange(aes(ymin=Null_pair_2.5pct_pp,ymax=Null_pair_97.5pct_pp),
                 linewidth=0.9,colour="#9CA9B3") +
  geom_point(aes(colour=Habitat),size=2.6) +
  geom_point(aes(y=Null_pair_expected_pp),shape=3,size=2.6,colour="#48515B") +
  facet_grid(Contrast ~ ., drop=FALSE) +
  scale_colour_manual(values=c("Bulk soil"="#596B85","Rhizosphere"="#AD693B")) +
  scale_y_continuous(expand=expansion(mult=c(0.04,0.04))) +
  coord_cartesian(ylim=c(-null_limit,null_limit)) +
  labs(title="B  Capability weighting",subtitle="Observed C − W and null distribution",
       x=NULL,y="C − W (percentage points)", colour="Habitat") +
  base_theme +
  theme(axis.text.x=element_text(angle=42,hjust=1,size=8),
        legend.position="bottom")

outpdf <- file.path(base,"Shared_MAG_habitat_and_capability_null_Breadth10.pdf")
outsvg <- file.path(base,"Shared_MAG_habitat_and_capability_null_Breadth10.svg")
render <- function() {
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(2,2,
    heights=unit(c(0.94,0.06),"null"),widths=unit(c(1,1),"null"))))
  print(a,vp=viewport(layout.pos.row=1,layout.pos.col=1),newpage=FALSE)
  print(b,vp=viewport(layout.pos.row=1,layout.pos.col=2),newpage=FALSE)
  grid.text("Each point compares two pooled metagenomes at one site and habitat. Only BS–RZ matched sites shown. Panel B: vertical bars = central 95% of 999 intact-profile permutations; short tick = null mean. No comparison was FDR-significant (q ≥ 0.16).",
            x=unit(0.025,"npc"),y=unit(0.02,"npc"),just=c("left","center"),
            gp=gpar(fontsize=8,col="#46505C"))
  upViewport()
}
pdf(outpdf,width=14,height=8.5,useDingbats=FALSE)
render()
dev.off()
svg(outsvg,width=14,height=8.5)
render()
dev.off()
cat("SUCCESS: Shared MAG and capability-null figure\n")
cat("PDF:",outpdf,"\nSVG:",outsvg,"\n")
cat("The vertical null bars are allocation-null ranges, NOT sampling confidence intervals.\n")
