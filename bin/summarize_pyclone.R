#!/usr/bin/env Rscript
#
# summarise_pyclone_vi.R
#
# Read a PyClone-VI results file and produce the cluster-level view:
#   <case>.clusters.tsv  one row per cluster, cellular prevalence per sample
#   <case>.ccf.pdf       prevalence at timepoint 1 vs timepoint 2
#
# Usage:
#   Rscript summarise_pyclone_vi.R CASE.results.tsv outdir [ref_sample]
#
# ref_sample is the x-axis sample (default: first alphabetically, i.e. BL).

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: Rscript summarise_pyclone_vi.R <results.tsv> <outdir> [ref_sample]")
}
res_file <- args[1]
outdir   <- args[2]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
case <- sub("\\.results\\.tsv$", "", basename(res_file))

res <- fread(res_file)
samples <- sort(unique(res$sample_id))
ref <- if (length(args) >= 3L) args[3] else samples[1]
alt <- setdiff(samples, ref)
if (!ref %in% samples) stop("ref_sample not in results: ", ref)

## --- assignment confidence --------------------------------------------------
# cluster_assignment_prob is per mutation; low-confidence mutations make a
# cluster look larger than the evidence supports.
if ("cluster_assignment_prob" %in% names(res)) {
  conf <- unique(res[, .(mutation_id, cluster_id, cluster_assignment_prob)])
  message(sprintf("Mutations with assignment prob < 0.8: %d / %d",
                  conf[cluster_assignment_prob < 0.8, .N], nrow(conf)))
}

## --- cluster table ----------------------------------------------------------

sizes <- unique(res[, .(mutation_id, cluster_id)])[, .(n_mut = .N), by = cluster_id]

wide <- dcast(unique(res[, .(cluster_id, sample_id, cellular_prevalence)]),
              cluster_id ~ sample_id, value.var = "cellular_prevalence")

clusters <- merge(wide, sizes, by = "cluster_id")
setorder(clusters, -n_mut)

# classify only for the two-timepoint case
if (length(samples) == 2L) {
  a <- clusters[[ref]]; b <- clusters[[alt]]
  clusters[, delta := b - a]
  clusters[, call := fifelse(a > 0.85 & b > 0.85, "clonal",
                     fifelse(b - a >  0.15,       "expanding",
                     fifelse(a - b >  0.15,       "contracting", "stable")))]
  clusters[a < 0.05 & b > 0.20, call := "emergent"]
  clusters[a > 0.20 & b < 0.05, call := "extinct"]
}

f <- file.path(outdir, paste0(case, ".clusters.tsv"))
fwrite(clusters, f, sep = "\t", quote = FALSE)
message("Wrote ", f)
print(clusters)

## --- plot -------------------------------------------------------------------

if (length(samples) == 2L) {
  # clusters with very few mutations are not interpretable regardless of shift
  clusters[, small := n_mut < 5]

  p <- ggplot(clusters, aes(x = .data[[ref]], y = .data[[alt]])) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
    geom_point(aes(size = n_mut, colour = call, alpha = !small)) +
    ggrepel::geom_text_repel(aes(label = paste0("C", cluster_id, " (n=", n_mut, ")")),
                             size = 3, min.segment.length = 0) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.35), guide = "none") +
    scale_size_continuous(range = c(2, 9)) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = paste("Cellular prevalence -", ref),
         y = paste("Cellular prevalence -", alt),
         title = case,
         subtitle = "above the diagonal = expanded at progression; faded = n<5 mutations") +
    theme_bw(base_size = 11)

  ggsave(file.path(outdir, paste0(case, ".ccf.pdf")), p, width = 7, height = 6)
  message("Wrote ", file.path(outdir, paste0(case, ".ccf.pdf")))
}

## --- mutation-to-cluster map, for linking to single-cell data ---------------

map <- unique(res[, .(mutation_id, cluster_id)])
fwrite(map, file.path(outdir, paste0(case, ".mutation_clusters.tsv")),
       sep = "\t", quote = FALSE)
message("Wrote mutation-to-cluster map (", nrow(map), " mutations)")