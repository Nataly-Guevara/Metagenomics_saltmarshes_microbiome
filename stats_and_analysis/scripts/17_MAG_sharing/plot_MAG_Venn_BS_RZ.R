#!/usr/bin/env Rscript
# Two condition-level Venn diagrams, using identical site coverage in BS and RZ.
# Run from /path/to/metagenomics_project
# Only R's bundled grid package is required.
#
# Each condition set = the UNION of MAGs detected in its selected-site samples.
# A MAG detected at both selected sites is counted once in that condition.
# Intersections show regional co-detection, NOT necessarily within-site sharing.
# All eligible MAGs, including Archaea, remain in the analysis.
# Shapes are schematic: areas are NOT proportional to counts.

suppressPackageStartupMessages(library(grid))

sites_use <- c("ACP", "BASS")
filters <- c(Breadth10 = 0.10, Breadth50 = 0.50)
AB <- "mapping/competitive_MAGs/MAG_abundance"
OUT <- "results/figures/MAG_Venn_BS_RZ"
condition_names <- c(NAT="Natural", BESE="BESE", BARE="BARE", DEGR="Degraded")
condition_colours <- c(NAT="#9BC1BC", BESE="#FF5A5F", BARE="#FFA630", DEGR="#5D576B")
conditions <- list(BS=c("NAT", "BESE", "BARE", "DEGR"), RZ=c("NAT", "BESE", "BARE"))

read_matrix <- function(path) {
  if (!file.exists(path)) stop("Input not found: ", path,
    "\nRun this script from the metagenomics-saltmarshes project directory.")
  d <- read.delim(path, check.names=FALSE, stringsAsFactors=FALSE,
                  quote="", comment.char="", fileEncoding="UTF-8-BOM")
  if (!("MAG_ID" %in% names(d))) stop("Missing MAG_ID in ", path)
  if (anyNA(d$MAG_ID) || any(!nzchar(d$MAG_ID)) || anyDuplicated(d$MAG_ID))
    stop("Missing or duplicated MAG IDs in ", path)
  x <- as.matrix(d[, setdiff(names(d), "MAG_ID"), drop=FALSE])
  suppressWarnings(storage.mode(x) <- "double")
  if (anyNA(x) || any(!is.finite(x)) || any(x < 0))
    stop("Non-numeric, missing, nonfinite or negative values in ", path)
  rownames(x) <- d$MAG_ID
  x
}

depth <- read_matrix(file.path(AB, "final_373_MAGs_mean_depth_matrix.tsv"))
breadth <- read_matrix(file.path(AB, "final_373_MAGs_breadth_matrix.tsv"))
if (!setequal(rownames(depth), rownames(breadth)) ||
    !setequal(colnames(depth), colnames(breadth))) stop("Depth/breadth identifiers differ.")
breadth <- breadth[rownames(depth), colnames(depth), drop=FALSE]
if (nrow(depth) != 373L || ncol(depth) != 20L) stop("Expected 373 MAGs and 20 samples.")
if (any(breadth > 1 + 1e-8)) stop("Genome breadth must be a fraction between 0 and 1.")

parts <- strsplit(colnames(depth), "_", fixed=TRUE)
if (any(lengths(parts) < 4L)) stop("Unexpected sample-name format.")
metadata <- data.frame(
  Sample=colnames(depth),
  Site=vapply(parts, `[`, character(1), 2L),
  Condition=vapply(parts, `[`, character(1), 3L),
  Habitat=vapply(parts, `[`, character(1), 4L), stringsAsFactors=FALSE
)
metadata$Condition[metadata$Condition == "NATURAL"] <- "NAT"
metadata$Condition[metadata$Condition %in% c("DEG", "DEGRADED")] <- "DEGR"
if (anyDuplicated(metadata[, c("Site", "Condition", "Habitat")]))
  stop("More than one sample per site x condition x habitat.")
for (h in names(conditions)) for (t in conditions[[h]]) for (s in sites_use) {
  found <- with(metadata, Site == s & Condition == t & Habitat == h)
  if (sum(found) != 1L) stop("Incomplete matched design: ", paste(s, t, h))
}

dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
scope <- paste(sites_use, collapse="_")
write_tsv <- function(x, name) write.table(x, file.path(OUT, name),
  sep="\t", quote=FALSE, row.names=FALSE, na="NA", fileEncoding="UTF-8")

# Schematic geometry; 4-set diagram uses ellipses (15 distinct intersections).
# Anchors below are validated against analytic ellipse membership, not inferred
# from a plotting package's area-vector order.
shape_geometry <- function(n) {
  if (n == 4L) {
    shapes <- rbind(c(.37,.46,.35,.20,135), c(.63,.46,.35,.20,45),
                    c(.49,.58,.34,.15,135), c(.51,.58,.34,.15,45))
    anchors <- rbind(
      c(.190,.562), c(.810,.562), c(.500,.270), c(.394,.760),
      c(.310,.668), c(.682,.458), c(.592,.386), c(.606,.760),
      c(.318,.458), c(.690,.668), c(.408,.386), c(.500,.692),
      c(.382,.570), c(.618,.570), c(.500,.484))
    labels <- rbind(c(.15,.15), c(.85,.15), c(.20,.91), c(.80,.91))
    leader_start <- rbind(c(.15,.205),c(.85,.205),c(.20,.858),c(.80,.858))
    leader_end <- rbind(c(.21,.32),c(.79,.32),c(.29,.79),c(.71,.79))
  } else if (n == 3L) {
    shapes <- rbind(c(.385,.57,.25,.25,0), c(.615,.57,.25,.25,0),
                    c(.50,.37,.25,.25,0))
    anchors <- rbind(c(.260,.636),c(.740,.636),c(.500,.698),
                     c(.500,.230),c(.332,.406),c(.668,.406),c(.500,.502))
    labels <- rbind(c(.21,.88),c(.79,.88),c(.50,.070))
    leader_start <- rbind(c(.21,.827),c(.79,.827),c(.50,.123))
    leader_end <- rbind(c(.29,.77),c(.71,.77),c(.50,.145))
  } else stop("Only 3- and 4-set diagrams are supported.")
  colnames(shapes) <- c("x", "y", "a", "b", "rotation")
  anchor_masks <- integer(nrow(anchors))
  for (j in seq_len(n)) {
    sh <- shapes[j,]; angle <- sh["rotation"] * pi/180
    u <- (anchors[,1]-sh["x"])*cos(angle) + (anchors[,2]-sh["y"])*sin(angle)
    v <- -(anchors[,1]-sh["x"])*sin(angle) + (anchors[,2]-sh["y"])*cos(angle)
    anchor_masks <- anchor_masks + as.integer((u/sh["a"])^2 + (v/sh["b"])^2 < 1)*2^(j-1)
  }
  if (!identical(as.integer(anchor_masks), seq_len(2^n-1L)))
    stop("Venn geometry/label membership validation failed.")
  list(shapes=shapes, anchors=anchors, labels=labels,
       leader_start=leader_start, leader_end=leader_end)
}

make_sets <- function(h, filt, present) {
  sets <- list(); sample_rows <- list()
  for (t in conditions[[h]]) {
    ix <- which(metadata$Site %in% sites_use & metadata$Condition == t & metadata$Habitat == h)
    ix <- ix[match(sites_use, metadata$Site[ix])]
    x <- present[, ix, drop=FALSE]
    sets[[t]] <- rownames(present)[rowSums(x) > 0L]
    sample_rows[[t]] <- data.frame(Filter=filt, Habitat=h,
      Condition=t, Sites=paste(sites_use,collapse=","), N_metagenomes=length(ix),
      MAGs_detected=length(sets[[t]]), Samples=paste(metadata$Sample[ix],collapse=";"))
  }
  presence <- do.call(cbind, lapply(sets, function(x) rownames(present) %in% x))
  dimnames(presence) <- list(rownames(present), names(sets))
  n <- length(sets)
  masks <- as.integer(presence %*% (2^(seq_len(n)-1L)))
  counts <- tabulate(masks, nbins=2^n-1L) # Mask zero is outside all selected sets.
  stopifnot(sum(counts) == sum(masks > 0L))
  for (j in seq_len(n)) {
    includes <- bitwAnd(seq_along(counts), 2^(j-1L)) > 0L
    stopifnot(sum(counts[includes]) == length(sets[[j]]))
  }
  combos <- vapply(seq_along(counts), function(m) {
    paste(names(sets)[bitwAnd(m, 2^(seq_len(n)-1L)) > 0L], collapse=" + ")
  }, character(1))
  membership <- data.frame(Filter=filt, Habitat=h, Sites=paste(sites_use,collapse=","),
    MAG_ID=rownames(present), Mask=masks, presence*1L, check.names=FALSE)
  membership$Exclusive_conditions <- ifelse(masks == 0L, "Not detected in selected samples", combos[pmax(masks,1L)])
  exclusive <- data.frame(Filter=filt, Habitat=h, Sites=paste(sites_use,collapse=","),
    Mask=seq_along(counts), Exclusive_conditions=combos, MAG_count=counts)
  inclusive <- do.call(rbind, lapply(combn(names(sets), 2L, simplify=FALSE), function(pair) {
    a <- sets[[pair[1]]]; b <- sets[[pair[2]]]
    shared <- intersect(a,b); united <- union(a,b)
    data.frame(Filter=filt, Habitat=h, Sites=paste(sites_use,collapse=","),
      Condition_A=pair[1], Condition_B=pair[2],
      MAGs_A=length(a), MAGs_B=length(b),
      Inclusive_shared_MAGs=length(shared),
      Jaccard=if (length(united)) length(shared)/length(united) else NA_real_)
  }))
  list(habitat=h, filter=filt, sets=sets, counts=counts, union=sum(counts),
       shared_all=tail(counts,1L), set_counts=do.call(rbind,sample_rows),
       membership=membership, exclusive=exclusive, pairwise=inclusive)
}

ellipse_vertices <- function(sh) {
  t <- seq(0, 2*pi, length.out=721L); angle <- sh["rotation"]*pi/180
  cbind(sh["x"] + sh["a"]*cos(t)*cos(angle) - sh["b"]*sin(t)*sin(angle),
        sh["y"] + sh["a"]*cos(t)*sin(angle) + sh["b"]*sin(t)*cos(angle))
}

draw_panel <- function(result, letter="") {
  n <- length(result$sets); geometry <- shape_geometry(n)
  title <- if (result$habitat == "BS") "Bulk soil (BS)" else "Rhizosphere (RZ)"
  title <- paste0(if (nzchar(letter)) paste0(letter, "  ") else "", title)
  grid.text(title, x=.5,y=.965, gp=gpar(fontfamily="sans",fontsize=17,fontface="bold"))
  grid.text(paste0(paste(sites_use,collapse=" + "), " | ", n,
                   " conditions | detection breadth >= ", round(100*filters[result$filter]), "%"),
            x=.5,y=.923,gp=gpar(fontfamily="sans",fontsize=10.5,col="#454545"))
  pushViewport(viewport(x=.5,y=.535,width=unit(.92,"snpc"),height=unit(.92,"snpc")))
  # Layer translucent fills first, then borders, to retain each set outline.
  for (j in seq_len(n)) {
    xy <- ellipse_vertices(geometry$shapes[j,])
    grid.polygon(xy[,1],xy[,2],gp=gpar(col=NA,
      fill=adjustcolor(condition_colours[names(result$sets)[j]], alpha.f=.19)))
  }
  for (j in seq_len(n)) {
    xy <- ellipse_vertices(geometry$shapes[j,])
    grid.polygon(xy[,1],xy[,2],gp=gpar(fill=NA,
      col=condition_colours[names(result$sets)[j]],lwd=2))
  }
  for (mask in seq_along(result$counts)) {
    ax <- geometry$anchors[mask,1]; ay <- geometry$anchors[mask,2]
    all_shared <- mask == length(result$counts)
    grid.roundrect(x=ax,y=ay,width=unit(if (all_shared) 11 else 9,"mm"),
      height=unit(5.8,"mm"),r=unit(1.1,"mm"),
      gp=gpar(fill=adjustcolor("white",alpha.f=.94),
              col=if (all_shared) "#37474F" else NA, lwd=.8))
    grid.text(as.character(result$counts[mask]),ax,ay,
      gp=gpar(fontfamily="sans",fontsize=12.5,
              fontface=if (all_shared) "bold" else "plain",
              col=if (result$counts[mask] == 0L) "#6B7280" else "#171717"))
  }
  for (j in seq_len(n)) {
    condition <- names(result$sets)[j]
    grid.lines(x=c(geometry$leader_start[j,1],geometry$leader_end[j,1]),
               y=c(geometry$leader_start[j,2],geometry$leader_end[j,2]),
               gp=gpar(col=condition_colours[condition],lwd=1.6))
    grid.text(paste0(condition_names[condition], "\nn = ", length(result$sets[[j]])),
      x=geometry$labels[j,1],y=geometry$labels[j,2],
      gp=gpar(fontfamily="sans",fontsize=12.3,fontface="bold",lineheight=1.12))
  }
  popViewport()
  grid.text(paste0("Detected union: ", result$union, " MAGs   |   Shared by all ",
                    n, ": ", result$shared_all), x=.5,y=.105,
            gp=gpar(fontfamily="sans",fontsize=10.8,fontface="bold"))
  grid.text("One condition set = detection at either selected site; MAG IDs counted once.",
            x=.5,y=.071,gp=gpar(fontfamily="sans",fontsize=8.5,col="#454545"))
  grid.text("Region numbers are exclusive counts. Shape areas are not proportional.",
            x=.5,y=.044,gp=gpar(fontfamily="sans",fontsize=8.5,col="#454545"))
}

save_figure <- function(stem, result_list) {
  combined <- length(result_list) == 2L
  width <- if (combined) 14 else 7
  height <- 7.7
  draw <- function() {
    grid.newpage()
    if (combined) {
      for (j in 1:2) {
        pushViewport(viewport(x=(j-.5)/2,y=.5,width=.5,height=1))
        draw_panel(result_list[[j]], letters[j]); popViewport()
      }
    } else draw_panel(result_list[[1]])
  }
  file <- file.path(OUT,paste0(stem,".pdf"))
  grDevices::pdf(file,width=width,height=height,useDingbats=FALSE,
                family="Helvetica",onefile=TRUE)
  tryCatch(draw(), finally=invisible(grDevices::dev.off()))
  if (capabilities("cairo")) {
    grDevices::svg(file.path(OUT,paste0(stem,".svg")),width=width,height=height,
                  family="sans",onefile=FALSE)
    tryCatch(draw(), finally=invisible(grDevices::dev.off()))
  }
}

for (filt in names(filters)) {
  detected <- depth > 0 & breadth >= filters[filt]
  result <- lapply(names(conditions), make_sets, filt=filt, present=detected)
  names(result) <- names(conditions)
  for (h in names(result)) {
    r <- result[[h]]
    tag <- paste(h,scope,filt,sep="_")
    write_tsv(r$set_counts,paste0("Venn_set_sizes_",tag,".tsv"))
    write_tsv(r$exclusive,paste0("Venn_exclusive_intersections_",tag,".tsv"))
    write_tsv(r$pairwise,paste0("Venn_inclusive_pairwise_counts_",tag,".tsv"))
    write_tsv(r$membership,paste0("Venn_MAG_memberships_",tag,".tsv"))
    save_figure(paste0("MAG_Venn_",tag),list(r))
    cat("\n",h,"|",filt,"| sites:",paste(sites_use,collapse=", "),"\n")
    print(r$set_counts[,c("Condition","N_metagenomes","MAGs_detected")],row.names=FALSE)
    cat("Union:",r$union,"| Shared by all",length(r$sets),"conditions:",r$shared_all,"\n")
  }
  save_figure(paste("MAG_Venn_BS_RZ",scope,filt,sep="_"),result)
}

selected <- metadata$Site %in% sites_use
write_tsv(metadata[selected,],paste0("Venn_selected_samples_",scope,".tsv"))
caption <- c(
  "Detected MAG overlap across ecological conditions in bulk soil and rhizosphere.",
  paste0("Both habitats are restricted to ",paste(sites_use,collapse=" and "),
    ", with one pooled-DNA metagenome per site x condition x habitat."),
  "Each condition set is the union of MAG IDs detected in either selected-site metagenome.",
  "A MAG is counted once per condition, even when detected at both sites.",
  "Detection requires positive mean depth and the indicated genome-breadth threshold.",
  "All 373 bacterial and archaeal catalogue MAGs are eligible; counts are not abundance weighted.",
  "Labels give condition totals; numbers inside regions give exclusive intersection counts.",
  "Pairwise-inclusive shared counts are also supplied separately, because they include multi-condition overlaps.",
  "Overlap denotes occurrence across these sites, not necessarily co-detection within the same site.",
  "TRUG is excluded from both panels because its BARE bulk-soil metagenome is unavailable.",
  "Degraded rhizosphere was not sampled. The four-condition BS and three-condition RZ cores are not equivalent comparisons.",
  "Shape areas are schematic and are not proportional to MAG counts.",
  "Condition-only detection does not establish ecological exclusivity or habitat specialisation.",
  "All original input data and earlier UpSet outputs are left unchanged."
)
writeLines(caption,file.path(OUT,"MAG_Venn_caption_and_definitions.txt"))
writeLines(capture.output(sessionInfo()),file.path(OUT,"R_sessionInfo.txt"))
cat("\nSUCCESS: Two habitat Venn diagrams, with matched site coverage\n")
cat("Primary combined: MAG_Venn_BS_RZ_",scope,"_Breadth10.pdf\n",sep="")
cat("Separate BS and RZ PDFs, Breadth50 sensitivity PDFs, count tables and MAG memberships saved.\n")
cat("Each condition set combines detections across sites; this is not a within-site overlap test.\n")
cat("Output:",normalizePath(OUT),"\n")
