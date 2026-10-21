#!/usr/bin/env Rscript
# Figure 2: GTDB taxonomy of the detected, reconstructed MAG community.
# Run from /path/to/metagenomics_project
# R >= 4.0; ggplot2 only. No annotation, mapping, tests or package installation.
# Primary: >=10% genome breadth. Sensitivity: >=50% genome breadth.
# All 373 catalogue MAGs, including Archaea, are eligible for detection.
# Relative abundance = retained MAG mean depth / sum(retained MAG mean depth).
# DO NOT use the assembly-origin treatment as the MAG's ecological occurrence.

suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(paste("ggplot2 is unavailable. On Gadi use R/4.4.2 and",
               "R_LIBS_USER=/path/to/local/Rlibs/4.4"), call. = FALSE)
  }
  library(ggplot2)
})

ROOT <- getwd()
OUT <- file.path(ROOT, "results/figures/Fig2_taxonomy")
BREADTH_FILTERS <- c(Breadth10 = 0.10, Breadth50 = 0.50)
TOP_PHYLA <- 12L
TOP_ORDERS <- 15L
CONDITIONS <- c("NAT", "BESE", "BARE", "DEGR")
SITES <- c("ACP", "BASS", "TRUG")
CONDITION_LABELS <- c(NAT = "Natural", BESE = "BESE", BARE = "BARE", DEGR = "Degraded")
CONDITION_COLOURS <- c(NAT = "#9BC1BC", BESE = "#FF5A5F", BARE = "#FFA630", DEGR = "#5D576B")

read_tsv <- function(relative, required) {
  path <- file.path(ROOT, relative)
  if (!file.exists(path)) stop("Missing input: ", path, call. = FALSE)
  x <- read.delim(path, sep = "\t", check.names = FALSE,
                  stringsAsFactors = FALSE, quote = "", comment.char = "",
                  na.strings = c("NA", "N/A"))
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("Missing columns in ", path, ": ", paste(missing, collapse = ", "))
  x
}
check_ids <- function(x, column) {
  if (anyNA(x[[column]]) || any(!nzchar(trimws(x[[column]]))) || anyDuplicated(x[[column]])) {
    stop("Missing or duplicated identifiers: ", column)
  }
}
write_tsv <- function(x, name) {
  write.table(x, file.path(OUT, name), sep = "\t", row.names = FALSE,
              quote = FALSE, na = "NA")
}

tax <- read_tsv("results/final_summaries/combined_binners/final_373_MAGs_GTDB_taxonomy.tsv",
                c("Final_MAG", "GTDB_classification"))
dep <- read_tsv("mapping/competitive_MAGs/MAG_abundance/final_373_MAGs_mean_depth_matrix.tsv", "MAG_ID")
brd <- read_tsv("mapping/competitive_MAGs/MAG_abundance/final_373_MAGs_breadth_matrix.tsv", "MAG_ID")
check_ids(tax, "Final_MAG"); check_ids(dep, "MAG_ID"); check_ids(brd, "MAG_ID")
stopifnot(nrow(dep) == 373L, nrow(brd) == 373L, nrow(tax) == 373L,
          setequal(dep$MAG_ID, brd$MAG_ID), setequal(dep$MAG_ID, tax$Final_MAG))
samples <- setdiff(names(dep), "MAG_ID")
stopifnot(length(samples) == 20L, setequal(samples, setdiff(names(brd), "MAG_ID")))
D <- as.matrix(dep[, samples, drop = FALSE]); storage.mode(D) <- "double"
B <- as.matrix(brd[match(dep$MAG_ID, brd$MAG_ID), samples, drop = FALSE]); storage.mode(B) <- "double"
rownames(D) <- rownames(B) <- dep$MAG_ID
stopifnot(all(is.finite(D)), all(D >= 0), all(is.finite(B)), all(B >= 0 & B <= 1))
tax <- tax[match(dep$MAG_ID, tax$Final_MAG), , drop = FALSE]

parts <- strsplit(samples, "_", fixed = TRUE)
if (any(lengths(parts) < 4L)) stop("Unexpected metagenome sample-name format.")
meta <- data.frame(Sample = samples,
                   Site = vapply(parts, `[`, character(1), 2),
                   Treatment = vapply(parts, `[`, character(1), 3),
                   Habitat = vapply(parts, `[`, character(1), 4))
meta$Treatment[meta$Treatment == "NATURAL"] <- "NAT"
meta$Treatment[meta$Treatment %in% c("DEG", "DEGRADED")] <- "DEGR"
cell <- function(x) paste(x$Site, x$Treatment, x$Habitat, sep = "|")
expected <- rbind(expand.grid(Site = SITES, Treatment = c("NAT", "BESE", "BARE"),
                             Habitat = c("BS", "RZ"), stringsAsFactors = FALSE),
                  data.frame(Site = SITES, Treatment = "DEGR", Habitat = "BS"))
expected <- expected[!(expected$Site == "TRUG" & expected$Treatment == "BARE" & expected$Habitat == "BS"), ]
stopifnot(!anyDuplicated(cell(meta)), setequal(cell(meta), cell(expected)))
meta <- meta[order(match(meta$Habitat, c("BS", "RZ")), match(meta$Treatment, CONDITIONS), match(meta$Site, SITES)), ]
D <- D[, meta$Sample, drop = FALSE]; B <- B[, meta$Sample, drop = FALSE]
meta$Habitat_label <- factor(meta$Habitat, c("BS", "RZ"), c("Bulk soil (BS)", "Rhizosphere (RZ)"))
meta$Cell <- paste(meta$Habitat, meta$Treatment, sep = "|")
meta$X <- NA_real_
headers <- list()
for (habitat in c("BS", "RZ")) {
  cursor <- 1
  for (condition in CONDITIONS) {
    ix <- which(meta$Habitat == habitat & meta$Treatment == condition)
    if (!length(ix)) next
    meta$X[ix] <- cursor + seq_along(ix) - 1
    headers[[length(headers) + 1L]] <- data.frame(
      Habitat = habitat, Treatment = condition,
      Left = min(meta$X[ix]) - 0.42, Right = max(meta$X[ix]) + 0.42,
      Mid = mean(range(meta$X[ix])), Label = CONDITION_LABELS[[condition]])
    cursor <- max(meta$X[ix]) + 1.65
  }
}
headers <- do.call(rbind, headers)
headers$Habitat_label <- factor(headers$Habitat, c("BS", "RZ"), c("Bulk soil (BS)", "Rhizosphere (RZ)"))

extract_rank <- function(classification, prefix) {
  if (is.na(classification) || !nzchar(classification)) return(NA_character_)
  tokens <- trimws(strsplit(classification, ";", fixed = TRUE)[[1]])
  tokens <- tokens[startsWith(tokens, prefix)]
  if (length(tokens) > 1L) stop("Duplicated rank in GTDB classification.")
  if (!length(tokens)) return(NA_character_)
  value <- substring(tokens, nchar(prefix) + 1L)
  if (nzchar(value)) value else NA_character_
}
tax$Domain <- vapply(tax$GTDB_classification, extract_rank, character(1), prefix = "d__")
tax$Phylum <- vapply(tax$GTDB_classification, extract_rank, character(1), prefix = "p__")
tax$Order <- vapply(tax$GTDB_classification, extract_rank, character(1), prefix = "o__")
rank_labels <- function(rank) {
  value <- tax[[rank]]
  unknown <- is.na(value) | !nzchar(value)
  value[unknown] <- "Unclassified"
  archaea <- !unknown & !is.na(tax$Domain) & tax$Domain == "Archaea"
  value[archaea] <- paste0(value[archaea], " [A]")
  value
}
relative_weights <- function(threshold) {
  retained <- D > 0 & B >= threshold
  x <- D; x[!retained] <- 0
  total <- colSums(x)
  if (any(total <= 0)) stop("No detected MAG abundance in a sample at breadth ", threshold)
  sweep(x, 2, total, "/") * 100
}
cell_means <- function(mat) {
  sapply(unique(meta$Cell), function(g) rowMeans(mat[, meta$Cell == g, drop = FALSE]))
}
primary <- relative_weights(BREADTH_FILTERS[["Breadth10"]])
selected <- list(); palettes <- list(); full_primary <- list()
for (rank in c("Phylum", "Order")) {
  full <- rowsum(primary, rank_labels(rank), reorder = FALSE)
  full_primary[[rank]] <- full
  # Equal-weight the seven observed condition x habitat cells for display ranking.
  score <- rowMeans(cell_means(full))
  candidate <- names(score)[names(score) != "Unclassified"]
  candidate <- candidate[order(-score[candidate], candidate)]
  selected[[rank]] <- head(candidate, if (rank == "Phylum") TOP_PHYLA else TOP_ORDERS)
  levels <- c(selected[[rank]], "Other", if ("Unclassified" %in% rownames(full)) "Unclassified")
  pal <- setNames(grDevices::hcl.colors(length(selected[[rank]]), "Dark 3"), selected[[rank]])
  if (rank == "Phylum") {
    stable <- c(Actinomycetota = "#CB8354", Pseudomonadota = "#467DAA", Bacteroidota = "#5EA59A",
      Gemmatimonadota = "#A791C5", Chloroflexota = "#AFB15E", "Halobacteriota [A]" = "#C66F94",
      Desulfobacterota = "#708C86", Acidobacteriota = "#9A755F", Cyanobacteriota = "#83B9D3",
      Myxococcota_A = "#9B809E", Planctomycetota = "#D9B661", Deinococcota = "#D69A93")
    common <- intersect(names(pal), names(stable)); pal[common] <- stable[common]
  }
  palettes[[rank]] <- c(pal, Other = "#DCE1E5", Unclassified = "#A9B1B9")[levels]
}

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
write_tsv(tax[, c("Final_MAG", "Domain", "Phylum", "Order", "GTDB_classification")], "Fig2_MAG_taxonomy_assignments.tsv")

plot_theme <- function() {
  theme_bw(base_size = 11) + theme(
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = "#E9ECEF", linewidth = 0.25),
    panel.border = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(),
    axis.text.y = element_text(colour = "#29313A"), axis.title.y = element_text(margin = margin(r = 8)),
    strip.background = element_rect(fill = "#EDF1F4", colour = NA),
    strip.text = element_text(face = "bold", size = 11),
    panel.spacing.x = grid::unit(0.6, "cm"),
    legend.position = "right", legend.key.size = grid::unit(0.38, "cm"),
    legend.text = element_text(size = 9), legend.title = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 10, margin = margin(b = 10)),
    plot.caption = element_text(size = 8, hjust = 0, margin = margin(t = 9)),
    plot.margin = margin(12, 14, 10, 10))
}
export_plot <- function(p, stem, width = 11.5, height = 6.3) {
  grDevices::pdf(file.path(OUT, paste0(stem, ".pdf")), width = width, height = height, useDingbats = FALSE)
  print(p); invisible(dev.off())
  if (capabilities("cairo")) {
    grDevices::svg(file.path(OUT, paste0(stem, ".svg")), width = width, height = height)
    print(p); invisible(dev.off())
    grDevices::png(file.path(OUT, paste0(stem, ".png")), width = width, height = height,
                   units = "in", res = 300, type = "cairo")
    print(p); invisible(dev.off())
  } else warning("Cairo is unavailable: PDF was saved, but SVG and PNG exports were skipped.")
}

for (filter_name in names(BREADTH_FILTERS)) {
  threshold <- BREADTH_FILTERS[[filter_name]]
  W <- relative_weights(threshold)
  sample_info <- meta
  sample_info$Detected_MAGs <- colSums(D > 0 & B >= threshold)
  sample_info$Retained_mean_depth <- colSums(D * (D > 0 & B >= threshold))
  sample_info$Retained_depth_pct <- 100 * sample_info$Retained_mean_depth / colSums(D)
  sample_info$Archaeal_relative_abundance_pct <- colSums(W[!is.na(tax$Domain) & tax$Domain == "Archaea", , drop = FALSE])
  write_tsv(sample_info, paste0("Fig2_sample_metadata_", filter_name, ".tsv"))
  for (rank in c("Phylum", "Order")) {
    raw <- rowsum(W, rank_labels(rank), reorder = FALSE)
    full_long <- data.frame(Sample = rep(colnames(raw), each = nrow(raw)),
                            Taxon = rep(rownames(raw), times = ncol(raw)),
                            Relative_abundance_pct = as.vector(raw))
    full_long <- cbind(full_long, meta[match(full_long$Sample, meta$Sample), c("Site", "Treatment", "Habitat")])
    write_tsv(full_long, paste0("Fig2_", rank, "_all_taxa_sample_values_", filter_name, ".tsv"))
    groups <- ifelse(rownames(raw) %in% selected[[rank]], rownames(raw),
                     ifelse(rownames(raw) == "Unclassified", "Unclassified", "Other"))
    mat <- rowsum(raw, groups, reorder = FALSE)
    stopifnot(max(abs(colSums(mat) - 100)) < 1e-7)
    levels <- names(palettes[[rank]])
    long <- data.frame(Sample = rep(colnames(mat), each = nrow(mat)),
                       Taxon = rep(rownames(mat), times = ncol(mat)),
                       Relative_abundance_pct = as.vector(mat))
    long <- cbind(long, meta[match(long$Sample, meta$Sample), setdiff(names(meta), "Sample"), drop = FALSE])
    long$Taxon <- factor(long$Taxon, levels = levels)
    caption <- paste0("Each bar represents one metagenome; abundances are normalised within the detected MAG catalogue.\n",
       "All 373 MAGs are eligible, including Archaea ([A]); non-retained MAGs have zero weight. ",
       "TRUG BARE BS is unavailable; Degraded has no RZ sample.")
    p <- ggplot(long, aes(X, Relative_abundance_pct, fill = Taxon)) +
      geom_col(width = 0.84, colour = "white", linewidth = 0.08, position = position_stack(reverse = TRUE)) +
      geom_segment(data = headers, aes(x = Left, xend = Right, y = 103.7, yend = 103.7, colour = Treatment),
                   inherit.aes = FALSE, linewidth = 1.7, show.legend = FALSE) +
      geom_text(data = headers, aes(x = Mid, y = 109, label = Label), inherit.aes = FALSE, size = 3.25, fontface = "bold") +
      geom_text(data = meta, aes(x = X, y = -4.2, label = Site), inherit.aes = FALSE, size = 2.9) +
      facet_grid(cols = vars(Habitat_label), scales = "free_x", space = "free_x") +
      scale_fill_manual(values = palettes[[rank]], limits = levels, drop = FALSE) +
      scale_colour_manual(values = CONDITION_COLOURS, guide = "none") +
      scale_y_continuous(breaks = seq(0, 100, 20), expand = expansion(mult = 0)) +
      scale_x_continuous(expand = expansion(add = 0.6)) +
      coord_cartesian(ylim = c(-7, 112), clip = "off") +
      labs(title = paste0("Taxonomic composition of saltmarsh MAG communities: ", tolower(rank)),
           subtitle = paste0("Site-level profiles by ecological condition and soil habitat | genome breadth >= ", threshold * 100, "%"),
           x = NULL, y = "Relative abundance within detected MAGs (%)", fill = rank, caption = caption) + plot_theme()
    stem <- if (rank == "Phylum") paste0("Fig2_Taxonomy_phylum_samples_", filter_name) else paste0("FigS_Taxonomy_order_samples_", filter_name)
    export_plot(p, stem, height = if (rank == "Phylum") 6.3 else 7.1)
    write_tsv(long, paste0(stem, "_plotted_values.tsv"))
    write_tsv(data.frame(Taxon = levels, Colour = unname(palettes[[rank]]),
                         Display_selection = "Top taxa by equal-weight mean across seven condition-habitat cells at Breadth10; fixed at Breadth50"),
              paste0("Fig2_", rank, "_palette.tsv"))

    if (rank == "Phylum") {
      means <- cell_means(mat)
      mean_long <- data.frame(Cell = rep(colnames(means), each = nrow(means)),
                              Taxon = rep(rownames(means), times = ncol(means)),
                              Relative_abundance_pct = as.vector(means))
      cells <- meta[!duplicated(meta$Cell), c("Cell", "Treatment", "Habitat", "Habitat_label")]
      cells$N_metagenomes <- vapply(cells$Cell, function(g) sum(meta$Cell == g), integer(1))
      cells$Mean_detected_MAGs <- vapply(cells$Cell, function(g) mean(sample_info$Detected_MAGs[meta$Cell == g]), numeric(1))
      cells$Label <- paste0(CONDITION_LABELS[cells$Treatment], "\n(n = ", cells$N_metagenomes, ")")
      cells$X <- ave(seq_len(nrow(cells)), cells$Habitat, FUN = seq_along)
      mean_long <- cbind(mean_long, cells[match(mean_long$Cell, cells$Cell), setdiff(names(cells), "Cell")])
      mean_long$Taxon <- factor(mean_long$Taxon, levels = levels)
      pm <- ggplot(mean_long, aes(X, Relative_abundance_pct, fill = Taxon)) +
        geom_col(width = 0.72, colour = "white", linewidth = 0.12, position = position_stack(reverse = TRUE)) +
        geom_text(data = cells, aes(X, -5.8, label = Label), inherit.aes = FALSE, size = 3.0, lineheight = 0.95) +
        facet_grid(cols = vars(Habitat_label), scales = "free_x", space = "free_x") +
        scale_fill_manual(values = palettes[[rank]], limits = levels, drop = FALSE) +
        scale_y_continuous(breaks = seq(0, 100, 20), expand = expansion(mult = 0)) +
        scale_x_continuous(expand = expansion(add = 0.6)) +
        coord_cartesian(ylim = c(-12, 102), clip = "off") +
        labs(title = "Taxonomic composition across ecological conditions and soil habitats",
          subtitle = paste0("Phylum-level means of sample-normalised MAG abundances | genome breadth >= ", threshold * 100, "%"),
          x = NULL, y = "Mean relative abundance within detected MAGs (%)", fill = rank,
          caption = paste0("Each bar is the arithmetic mean of the available site-level profiles; n = metagenomes, not MAGs.\n",
          "BARE bulk soil includes ACP and BASS only. Archaea are retained ([A]); the sample-level figure shows site variation.")) + plot_theme()
      export_plot(pm, paste0("Fig2_Taxonomy_phylum_condition_means_", filter_name), width = 10.5, height = 6.3)
      write_tsv(mean_long, paste0("Fig2_phylum_condition_means_", filter_name, ".tsv"))
      write_tsv(cells, paste0("Fig2_condition_sample_counts_", filter_name, ".tsv"))
    }
  }
}

caption <- c("Figure 2. Taxonomic composition of reconstructed saltmarsh MAG communities.",
  "GTDB phylum-level composition is shown separately for bulk soil and rhizosphere, with metagenomes grouped by ecological condition and labelled by site. Each primary-figure bar represents one metagenome.",
  "A MAG was retained when its mean depth was positive and at least 10% of its reconstructed genome length was covered. Relative mean-depth abundance was calculated within each sample by dividing each retained MAG's mean depth by the sum of mean depths of all retained MAGs, including bacterial and archaeal MAGs.",
  "The 12 highest-ranked phyla are selected using the equally weighted mean of the seven observed condition-by-habitat means. Remaining classified phyla are grouped as Other, and missing rank assignments are retained as Unclassified. Archaeal taxa are marked [A]. Identical display taxa and colours are used at the 50% breadth sensitivity threshold.",
  "The companion seven-bar summary uses arithmetic means of site-level relative abundances. All condition-by-habitat cells contain three metagenomes except BARE bulk soil (ACP and BASS; n = 2). No Degraded rhizosphere was sampled, and TRUG BARE bulk soil was unavailable.",
  "The proportions describe relative abundance within the detected MAG catalogue, not MAG assembly yield, absolute cell abundance, or the complete environmental microbiome. Mean detected-MAG counts and retained-depth summaries are provided separately; MAG counts are not biological replicate counts.")
writeLines(caption, file.path(OUT, "Fig2_caption_draft.txt"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "Fig2_R_sessionInfo.txt"))
cat("\nSUCCESS: Figure 2 taxonomy plots\n")
cat("Catalogue MAGs:", nrow(tax), "\n")
print(table(tax$Domain, useNA = "ifany"))
cat("Metagenomes:", nrow(meta), "(BS:", sum(meta$Habitat == "BS"), "; RZ:", sum(meta$Habitat == "RZ"), ")\n")
cat("Main: Fig2_Taxonomy_phylum_samples_Breadth10.pdf\n")
cat("Companion: Fig2_Taxonomy_phylum_condition_means_Breadth10.pdf\n")
cat("Supplement: FigS_Taxonomy_order_samples_Breadth10.pdf\n")
cat("Matching Breadth50 plots, source tables, caption and session info saved.\n")
cat("Output directory:", normalizePath(OUT), "\n")
