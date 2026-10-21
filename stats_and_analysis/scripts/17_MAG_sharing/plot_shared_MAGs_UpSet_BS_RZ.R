#!/usr/bin/env Rscript
# Shared-MAG UpSet-style plots: condition intersections WITHIN site and habitat.
# Uses only R's bundled grid, graphics-device and table-reading functions.
# No package installation, remapping, reannotation or statistical testing.
# Run from the metagenomics-saltmarshes project root. Optional first argument: root.
#
# Detection: positive mean depth AND breadth >= 0.10 / 0.50.
# Each displayed set is ONE metagenome, not the MAG's assembly-origin group.
# Top bars: EXCLUSIVE intersections. Inactive dots mean not detected under the
# selected rule, not proven biological absence. All 373 MAGs remain eligible.
# Sources also include INCLUSIVE intersections for comparison with zeta counts.
# Plot design reference: Lex et al. (2014), UpSet: Visualization of Intersecting Sets.
# https://doi.org/10.1109/TVCG.2014.2346248

suppressPackageStartupMessages(library(grid))
args <- commandArgs(trailingOnly = TRUE)
ROOT <- if (length(args)) normalizePath(args[1], mustWork = TRUE) else getwd()
OUT <- file.path(ROOT, "results/figures/Shared_MAGs_UpSet")
SITES <- c("ACP", "BASS", "TRUG")
HABITATS <- c("BS", "RZ")
CONDITIONS <- c("NAT", "BESE", "BARE", "DEGR")
FILTERS <- c(Breadth10 = 0.10, Breadth50 = 0.50)
LABELS <- c(NAT = "Natural", BESE = "BESE", BARE = "BARE", DEGR = "Degraded")
CONDITION_COLOURS <- c(NAT = "#9BC1BC", BESE = "#FF5A5F",
                      BARE = "#FFA630", DEGR = "#5D576B")
INK <- "#253440"
BAR_COLOURS <- c(One_condition = "#A9B1B9", Multiple_conditions = "#426780",
                All_sampled_conditions = "#198477")

read_tsv <- function(relative, required) {
  path <- file.path(ROOT, relative)
  if (!file.exists(path)) stop("Missing input: ", path, call. = FALSE)
  x <- read.delim(path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE,
                  quote = "", comment.char = "", na.strings = c("NA", "N/A", ""))
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "),
                            " in ", path, call. = FALSE)
  x
}
check_ids <- function(x) {
  if (anyNA(x$MAG_ID) || any(!nzchar(trimws(x$MAG_ID))) || anyDuplicated(x$MAG_ID))
    stop("Missing or duplicated MAG_ID values.", call. = FALSE)
}
write_tsv <- function(x, name) {
  write.table(x, file.path(OUT, name), sep = "\t", quote = FALSE,
              row.names = FALSE, na = "NA")
}

Dtab <- read_tsv("mapping/competitive_MAGs/MAG_abundance/final_373_MAGs_mean_depth_matrix.tsv", "MAG_ID")
Btab <- read_tsv("mapping/competitive_MAGs/MAG_abundance/final_373_MAGs_breadth_matrix.tsv", "MAG_ID")
check_ids(Dtab); check_ids(Btab)
stopifnot(nrow(Dtab) == 373L, nrow(Btab) == 373L,
          setequal(Dtab$MAG_ID, Btab$MAG_ID))
samples <- setdiff(names(Dtab), "MAG_ID")
stopifnot(length(samples) == 20L,
          setequal(samples, setdiff(names(Btab), "MAG_ID")))
D <- as.matrix(Dtab[, samples, drop = FALSE]); storage.mode(D) <- "double"
B <- as.matrix(Btab[match(Dtab$MAG_ID, Btab$MAG_ID), samples, drop = FALSE])
storage.mode(B) <- "double"
rownames(D) <- rownames(B) <- Dtab$MAG_ID
stopifnot(all(is.finite(D)), all(D >= 0), all(is.finite(B)), all(B >= 0 & B <= 1))

parts <- strsplit(samples, "_", fixed = TRUE)
if (any(lengths(parts) < 4L)) stop("Unexpected metagenome name format.")
meta <- data.frame(Sample = samples,
                   Site = vapply(parts, `[`, character(1), 2),
                   Treatment = vapply(parts, `[`, character(1), 3),
                   Habitat = vapply(parts, `[`, character(1), 4),
                   stringsAsFactors = FALSE)
meta$Treatment[meta$Treatment == "NATURAL"] <- "NAT"
meta$Treatment[meta$Treatment %in% c("DEG", "DEGRADED")] <- "DEGR"
cell_key <- function(x) paste(x$Site, x$Treatment, x$Habitat, sep = "|")
expected <- rbind(
  expand.grid(Site = SITES, Treatment = c("NAT", "BESE", "BARE"),
              Habitat = HABITATS, stringsAsFactors = FALSE),
  data.frame(Site = SITES, Treatment = "DEGR", Habitat = "BS"))
expected <- expected[!(expected$Site == "TRUG" & expected$Treatment == "BARE" &
                       expected$Habitat == "BS"), , drop = FALSE]
if (anyDuplicated(cell_key(meta)) || !setequal(cell_key(meta), cell_key(expected)))
  stop("Sample design differs from expected 20 cells. Do not fill missing samples with zeros.")
meta <- meta[order(match(meta$Habitat, HABITATS), match(meta$Site, SITES),
                   match(meta$Treatment, CONDITIONS)), , drop = FALSE]
D <- D[, meta$Sample, drop = FALSE]; B <- B[, meta$Sample, drop = FALSE]

# Enumerate all nonempty membership patterns; a MAG belongs to one exact pattern.
# Code 0 means not detected in ANY available set and is excluded from UpSet bars.
make_panel <- function(filter, threshold, site, habitat) {
  m <- meta[meta$Site == site & meta$Habitat == habitat, , drop = FALSE]
  m <- m[order(match(m$Treatment, CONDITIONS)), , drop = FALSE]
  P <- D[, m$Sample, drop = FALSE] > 0 & B[, m$Sample, drop = FALSE] >= threshold
  colnames(P) <- m$Treatment
  k <- ncol(P)
  powers <- 2^(seq_len(k) - 1L)
  codes <- as.integer(P %*% powers)
  maxcode <- 2^k - 1L
  pattern <- sapply(powers, function(bit) (seq_len(maxcode) %/% bit) %% 2L == 1L)
  colnames(pattern) <- m$Treatment
  count <- tabulate(codes, nbins = maxcode)
  inclusive <- vapply(seq_len(maxcode), function(j) {
    take <- which(pattern[j, ])
    sum(rowSums(P[, take, drop = FALSE]) == length(take))
  }, numeric(1))
  totals <- colSums(P)
  union_n <- sum(codes > 0L)
  all_n <- sum(rowSums(P) == k)
  stopifnot(sum(count) == union_n, count[maxcode] == all_n,
            all(colSums(pattern * count) == totals))
  list(Filter = filter, Site = site, Habitat = habitat, meta = m, P = P,
       codes = codes, pattern = pattern, count = count, inclusive = inclusive,
       totals = totals, union_n = union_n, all_n = all_n, maxcode = maxcode)
}

panels <- list()
for (f in names(FILTERS)) for (site in SITES) for (h in HABITATS)
  panels[[paste(f, site, h, sep = "|")]] <- make_panel(f, FILTERS[[f]], site, h)

# Fixed intersection ordering per site/habitat across both detection thresholds.
# Show patterns observed at EITHER threshold, plus all-sampled-condition intersection
# even if its count is zero. Nothing is truncated to a top-N list.
for (site in SITES) for (h in HABITATS) {
  p10 <- panels[[paste("Breadth10", site, h, sep = "|")]]
  p50 <- panels[[paste("Breadth50", site, h, sep = "|")]]
  show <- which(p10$count > 0L | p50$count > 0L | seq_along(p10$count) == p10$maxcode)
  show <- show[order(-p10$count[show], -rowSums(p10$pattern)[show], show)]
  for (f in names(FILTERS)) panels[[paste(f, site, h, sep = "|")]]$show <- show
}

# Cross-check previously computed whole-MAG pairwise sharing where available.
baseline_path <- "results/ecological_strategies/functional_carrier_taxonomy/carrier_turnover/Whole_MAG_community_turnover.tsv"
checks <- 0L
if (file.exists(file.path(ROOT, baseline_path))) {
  baseline <- read_tsv(baseline_path, c("Filter", "Site", "Habitat", "Contrast", "Shared_MAGs"))
  for (i in seq_len(nrow(baseline))) {
    r <- baseline[i, , drop = FALSE]
    key <- paste(r$Filter, r$Site, r$Habitat, sep = "|")
    if (!key %in% names(panels)) next
    p <- panels[[key]]
    ab <- strsplit(r$Contrast, "_vs_", fixed = TRUE)[[1]]
    ab[ab %in% c("DEG", "DEGRADED")] <- "DEGR"; ab[ab == "NATURAL"] <- "NAT"
    if (length(ab) != 2L || !all(ab %in% colnames(p$P))) stop("Unexpected contrast in baseline.")
    shared <- sum(p$P[, ab[1]] & p$P[, ab[2]])
    union_n <- sum(p$P[, ab[1]] | p$P[, ab[2]])
    stopifnot(shared == as.integer(r$Shared_MAGs))
    if ("Community_Jaccard" %in% names(r) && union_n > 0L)
      stopifnot(abs(shared / union_n - as.numeric(r$Community_Jaccard)) < 1e-6)
    checks <- checks + 1L
  }
}

# Exact source tables, including membership of MAGs absent from an entire panel.
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
intersection_rows <- list(); membership_rows <- list(); set_rows <- list(); summary_rows <- list()
for (p in panels) {
  k <- ncol(p$P)
  for (j in seq_len(p$maxcode)) {
    included <- colnames(p$P)[p$pattern[j, ]]
    excluded <- colnames(p$P)[!p$pattern[j, ]]
    intersection_rows[[length(intersection_rows) + 1L]] <- data.frame(
      Filter = p$Filter, Site = p$Site, Habitat = p$Habitat,
      Available_conditions = paste(colnames(p$P), collapse = ","),
      Conditions_detected = paste(included, collapse = ","),
      Conditions_not_detected = paste(excluded, collapse = ","),
      N_conditions_detected = length(included), Membership_code = j,
      Exclusive_MAG_count = p$count[j], Inclusive_MAG_count = p$inclusive[j],
      Displayed = as.integer(j %in% p$show), Plot_column = match(j, p$show),
      stringsAsFactors = FALSE)
  }
  membership <- data.frame(Filter = p$Filter, Site = p$Site, Habitat = p$Habitat,
    MAG_ID = rownames(p$P), Membership_code = p$codes,
    N_conditions_detected = rowSums(p$P), Included_in_union = as.integer(p$codes > 0L),
    stringsAsFactors = FALSE)
  for (condition in CONDITIONS) {
    membership[[condition]] <- if (condition %in% colnames(p$P))
      as.integer(p$P[, condition]) else NA_integer_
  }
  membership_rows[[length(membership_rows) + 1L]] <- membership
  set_rows[[length(set_rows) + 1L]] <- data.frame(
    Filter = p$Filter, p$meta, Detected_MAGs = unname(p$totals),
    stringsAsFactors = FALSE)
  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    Filter = p$Filter, Site = p$Site, Habitat = p$Habitat, N_metagenomes = k,
    Available_conditions = paste(colnames(p$P), collapse = ","),
    MAG_union = p$union_n, MAGs_shared_all_sampled_conditions = p$all_n,
    MAGs_not_detected_in_any_set = nrow(p$P) - p$union_n,
    Nonzero_exclusive_intersections = sum(p$count > 0L), stringsAsFactors = FALSE)
}
write_tsv(do.call(rbind, intersection_rows), "Shared_MAGs_intersection_counts.tsv")
write_tsv(do.call(rbind, membership_rows), "Shared_MAGs_membership_by_site_habitat.tsv")
write_tsv(do.call(rbind, set_rows), "Shared_MAGs_sample_set_sizes.tsv")
summary_table <- do.call(rbind, summary_rows)
write_tsv(summary_table, "Shared_MAGs_panel_summary.tsv")

# Shared count axes across all sites, habitats and detection thresholds.
count_max <- max(vapply(panels, function(p) max(p$count), numeric(1)))
set_max <- max(vapply(panels, function(p) max(p$totals), numeric(1)))
YMAX <- max(10, tail(pretty(c(0, count_max), n = 4), 1))
SMAX <- max(10, tail(pretty(c(0, set_max), n = 3), 1))

text_at <- function(label, x, y, size = 9, just = "left", colour = INK,
                    fontface = "plain", rot = 0) {
  grid.text(label, x = x, y = y, just = just, rot = rot,
            gp = gpar(fontsize = size, col = colour, fontface = fontface, fontfamily = "sans"))
}

# Draw a classic UpSet layout inside the current viewport.
draw_panel <- function(p, letter) {
  show <- p$show; k <- ncol(p$P); n <- length(show)
  text_at(paste0(letter, "  ", p$Site), 0.015, 0.968, 13, fontface = "bold")
  text_at(sprintf("Union: %d MAGs  |  Shared by all %d sampled conditions: %d",
                  p$union_n, k, p$all_n), 0.015, 0.905, 8.5)
  if (p$Site == "TRUG" && p$Habitat == "BS")
    text_at("BARE BS unavailable.",
            0.015, 0.85, 8, colour = "#6A7178")

  # Top intersection bars, native coordinates mapped into this panel.
  left <- 0.365; right <- 0.98; y0 <- 0.405; ytop <- 0.805
  step <- (right - left) / n
  xs <- left + (seq_len(n) - 0.5) * step
  axis_y <- function(value) y0 + (ytop - y0) * value / YMAX
  ticks <- pretty(c(0, YMAX), n = 4); ticks <- ticks[ticks >= 0 & ticks <= YMAX]
  for (v in ticks) {
    yy <- axis_y(v)
    grid.lines(x = c(left, right), y = c(yy, yy), gp = gpar(col = "#E6E9EC", lwd = 0.7))
    text_at(format(v, trim = TRUE), left - 0.012, yy, 8, just = "right")
  }
  text_at("Exclusive MAG count", left, 0.86, 8.5, fontface = "bold")
  degree <- rowSums(p$pattern)[show]
  category <- ifelse(degree == k, "All_sampled_conditions",
                     ifelse(degree == 1L, "One_condition", "Multiple_conditions"))
  heights <- (ytop - y0) * p$count[show] / YMAX
  grid.rect(x = xs, y = y0 + heights / 2, width = step * 0.7, height = heights,
             gp = gpar(fill = unname(BAR_COLOURS[category]), col = NA))
  for (j in seq_len(n)) text_at(as.character(p$count[show[j]]), xs[j],
      axis_y(p$count[show[j]]) + 0.020, 8.2, just = "centre", fontface = "bold")

  # Set sizes and the membership matrix share row coordinates.
  ys <- seq(0.305, 0.115, length.out = k)
  row_h <- 0.046
  text_at("Detected MAGs", 0.16, 0.357, 8, fontface = "bold")
  for (i in seq_len(k)) {
    if (i %% 2L == 1L) grid.rect(x = 0.496, y = ys[i], width = 0.968, height = row_h,
                                gp = gpar(fill = "#F3F5F6", col = NA))
    condition <- colnames(p$P)[i]
    text_at(LABELS[[condition]], 0.145, ys[i], 9, just = "right")
    width <- 0.13 * p$totals[i] / SMAX
    grid.rect(x = 0.16 + width / 2, y = ys[i], width = width, height = 0.030,
               gp = gpar(fill = CONDITION_COLOURS[[condition]], col = NA))
    text_at(as.character(p$totals[i]), 0.301, ys[i], 8.2)
    grid.points(x = xs, y = rep(ys[i], n), pch = 16, size = unit(1.65, "mm"),
                 gp = gpar(col = "#D4D9DE"))
  }
  # A common scale is used for the set-size bars too.
  for (v in c(0, SMAX)) text_at(format(v, trim = TRUE), 0.16 + 0.13 * v / SMAX,
                              0.059, 7.5, just = "centre", colour = "#68717A")
  grid.lines(x = c(0.16, 0.29), y = c(0.08, 0.08), gp = gpar(col = "#B0B8BF", lwd = 0.7))
  for (j in seq_len(n)) {
    active <- which(p$pattern[show[j], ])
    if (length(active) > 1L)
      grid.lines(x = rep(xs[j], 2), y = range(ys[active]), gp = gpar(col = INK, lwd = 1.15))
    grid.points(x = rep(xs[j], length(active)), y = ys[active], pch = 16,
                 size = unit(2.25, "mm"), gp = gpar(col = INK))
  }
  grid.lines(x = c(0.01, 0.985), y = c(0.012, 0.012), gp = gpar(col = "#DEE3E7", lwd = 0.8))
}

draw_page <- function(filter, habitats = HABITATS) {
  grid.newpage()
  title <- if (length(habitats) == 2L) "Shared MAGs across saltmarsh conditions" else
    paste0("Shared MAGs across conditions: ", if (habitats == "BS") "bulk soil" else "rhizosphere")
  text_at(title, 0.03, 0.972, 17, fontface = "bold")
  text_at(sprintf("Site-resolved UpSet plots | Genome breadth >= %d%% and positive mean depth | All 373 catalogue MAGs eligible",
                  round(100 * FILTERS[[filter]])), 0.03, 0.943, 9.5)
  ncol <- length(habitats); gap <- 0.025
  panel_width <- (0.94 - (ncol - 1L) * gap) / ncol
  for (cc in seq_along(habitats)) {
    xleft <- 0.03 + (cc - 1L) * (panel_width + gap)
    grid.rect(x = xleft + panel_width / 2, y = 0.910, width = panel_width, height = 0.029,
               gp = gpar(fill = "#E9EEF1", col = NA))
    text_at(if (habitats[cc] == "BS") "Bulk soil (BS)" else "Rhizosphere (RZ)",
            xleft + panel_width / 2, 0.910, 11, just = "centre", fontface = "bold")
    for (rr in seq_along(SITES)) {
      ytop <- 0.883 - (rr - 1L) * 0.257
      pushViewport(viewport(x = xleft, y = ytop, width = panel_width, height = 0.250,
                            just = c("left", "top"), clip = "off"))
      index <- (rr - 1L) * 2L + match(habitats[cc], HABITATS)
      draw_panel(panels[[paste(filter, SITES[rr], habitats[cc], sep = "|")]], letters[index])
      popViewport()
    }
  }
  legend_x <- c(0.04, 0.29, 0.54)
  legend_text <- c("One condition only", "Multiple, but not all conditions", "All sampled conditions")
  for (j in 1:3) {
    grid.rect(x = legend_x[j], y = 0.080, width = 0.012, height = 0.010,
               gp = gpar(fill = BAR_COLOURS[j], col = NA))
    text_at(legend_text[j], legend_x[j] + 0.011, 0.080, 8.5)
  }
  footer <- paste(
    "Each bar counts MAGs detected in exactly the conditions marked by filled dots. Counts are not abundance weighted.",
    "Each set is one pooled-DNA metagenome; intersections are within sites. Degraded RZ was not sampled; TRUG BARE BS is unavailable.",
    "Detection-only exclusivity is not proof of habitat restriction. Four-condition BS and three-condition RZ intersections are different comparisons.",
    sep = "\n")
  text_at(footer, 0.03, 0.051, 8, just = c("left", "top"))
}

export_page <- function(stem, draw, width, height = 11.6) {
  grDevices::pdf(file.path(OUT, paste0(stem, ".pdf")), width = width, height = height,
                 useDingbats = FALSE, onefile = TRUE)
  tryCatch(draw(), finally = invisible(grDevices::dev.off()))
  if (capabilities("cairo")) {
    grDevices::svg(file.path(OUT, paste0(stem, ".svg")), width = width, height = height)
    tryCatch(draw(), finally = invisible(grDevices::dev.off()))
  }
}
for (f in names(FILTERS)) {
  export_page(paste0("Shared_MAGs_UpSet_BS_RZ_", f), function() draw_page(f), width = 17)
  for (h in HABITATS)
    export_page(paste0("Shared_MAGs_UpSet_", h, "_", f), function() draw_page(f, h), width = 10.5)
}

caption <- c(
  "Shared MAG identities across ecological conditions within each site and soil habitat.",
  "UpSet-style panels show exact intersections among Natural, BESE, BARE and Degraded metagenomes where sampled, separately for ACP, BASS and TRUG and for bulk soil (BS) and rhizosphere (RZ). Each set contains MAGs with positive mean depth and genome breadth at or above the stated threshold. The full 373-MAG catalogue, including bacterial and archaeal MAGs, was eligible. Each set corresponds to one metagenome produced by pooling aliquots from five replicate DNA extracts; individual MAGs and intersections are not biological replicates.",
  "Vertical bars give exclusive intersection counts: filled dots identify the conditions in which the counted MAGs were detected, and light dots the other sampled conditions in which they did not meet the detection criterion. Horizontal bars give the total detected-MAG count for each metagenome. Colours distinguish single-condition, multi-condition and all-sampled-condition intersections. Patterns observed at either 10% or 50% breadth, plus the all-sampled-condition pattern, are displayed in the same order at both thresholds; zero counts are retained where relevant. Both count axes have common scales across panels and thresholds.",
  "Degraded rhizosphere was not sampled, and the TRUG BARE bulk-soil metagenome was unavailable; neither is inserted as an empty set. The all-condition BS intersection involves four conditions at ACP and BASS, but three at TRUG; RZ intersections involve Natural, BESE and BARE. Consequently, raw all-condition counts are not a matched habitat-effect test. No sites are pooled into the plotted condition sets. Counts describe detection of reconstructed MAG identities, not absolute biomass, functional equivalence, habitat specificity or longitudinal recovery.",
  "The intersection source table also includes inclusive counts. For example, an inclusive Natural-BESE count includes MAGs also detected in BARE or Degraded, whereas the Natural-BESE exclusive bar excludes them. The complete all-set intersection corresponds to the highest available within-panel zeta order.",
  "Reference: Lex A, Gehlenborg N, Strobelt H, Vuillemot R, Pfister H. 2014. UpSet: Visualization of Intersecting Sets. IEEE Transactions on Visualization and Computer Graphics 20:1983-1992. doi:10.1109/TVCG.2014.2346248.",
  "Implementation: R grid graphics; no external UpSet package required. No statistical tests were added.")
writeLines(caption, file.path(OUT, "Shared_MAGs_UpSet_caption_and_definitions.txt"))
capture.output(sessionInfo(), file = file.path(OUT, "Shared_MAGs_UpSet_sessionInfo.txt"))

cat("\nSUCCESS: Site-resolved shared-MAG UpSet plots\n")
cat("Catalogue MAGs:", nrow(D), "\n")
cat("Metagenomes:", ncol(D), "(BS:", sum(meta$Habitat == "BS"), "; RZ:", sum(meta$Habitat == "RZ"), ")\n")
if (checks) cat("Existing whole-MAG pairwise counts verified:", checks, "\n")
cat("\nPRIMARY 10% BREADTH\n")
print(subset(summary_table, Filter == "Breadth10",
             select = c(Site, Habitat, N_metagenomes, MAG_union, MAGs_shared_all_sampled_conditions)), row.names = FALSE)
cat("\nCombined: Shared_MAGs_UpSet_BS_RZ_Breadth10.pdf\n")
cat("Separate: Shared_MAGs_UpSet_BS_Breadth10.pdf\n")
cat("Separate: Shared_MAGs_UpSet_RZ_Breadth10.pdf\n")
cat("Matching Breadth50 PDFs, SVGs (if Cairo available), count tables, MAG memberships and caption saved.\n")
cat("Output:", OUT, "\n")
cat("Top bars are exclusive intersections, not inclusive pairwise shared counts.\n")
