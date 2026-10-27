#!/usr/bin/env Rscript
#
# build_pyclone_input.R  (union strategy)
#
# Strategy:
#   1. Call SNVs per sample, keep PASS.
#   2. Take the UNION of PASS SNVs across all samples of a patient.
#   3. Extract allele-specific read counts for every union site in EVERY
#      sample, whether or not it was called in that sample.
#   4. Assemble the PyClone table from those counts + per-sample ASCAT.
#
# This script covers step 2 (union construction) and step 4 (table assembly).
# Step 3 needs BAM/CRAM access and is done by an allele counter -- use
# --emit-sites to write the site list to feed it.
#
# Usage:
#   # (a) build the union site list from per-sample PASS VCFs
#   Rscript build_pyclone_input.R --emit-sites manifest.tsv outdir
#
#   # ... run the allele counter at those sites in every sample ...
#
#   # (b) assemble the PyClone input from the force-called counts
#   Rscript build_pyclone_input.R manifest.tsv outdir [classic|vi]
#
# Manifest (tab-separated), one row per TUMOUR SAMPLE:
#   CASE_ID     patient identifier
#   SAMPLE      sample identifier, unique within patient (e.g. timepoint)
#   counts      force-called allele counts at union sites:
#               CHROM POS REF ALT REF_COUNT VAR_COUNT   (no header)
#   ascat_cnv   this sample's ASCAT *.cnvs.txt
#   ascat_pp    this sample's ASCAT *.purityploidy.txt
#   pass_vcf    (only for --emit-sites) per-sample filtered VCF
#   sex         optional: "XX" or "XY"

suppressPackageStartupMessages(library(data.table))

## ---------------------------------------------------------------- options ---

KEEP_CHRS   <- as.character(1:22)  # add "X" only with sex set in the manifest
DROP_MAJOR0 <- TRUE                # PyClone requires major_cn >= 1
MIN_DEPTH   <- 20                  # per-sample depth floor at a union site.
                                   # A site below this in ANY sample is dropped:
                                   # an undersampled zero is not evidence of
                                   # absence, and PyClone would read it as a
                                   # confident CCF of 0.

## ------------------------------------------------------------- arguments ---

args <- commandArgs(trailingOnly = TRUE)
emit_sites <- "--emit-sites" %in% args
args <- setdiff(args, "--emit-sites")

if (length(args) < 2L) {
  stop("Usage:\n",
       "  Rscript build_pyclone_input.R --emit-sites <manifest.tsv> <outdir>\n",
       "  Rscript build_pyclone_input.R <manifest.tsv> <outdir> [classic|vi]")
}
manifest_file <- args[1]
outdir        <- args[2]
fmt           <- if (length(args) >= 3L) tolower(args[3]) else "vi"
if (!fmt %in% c("classic", "vi")) stop("Format must be 'classic' or 'vi'")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

manifest <- fread(manifest_file, colClasses = "character")
need <- if (emit_sites) c("CASE_ID", "SAMPLE", "pass_vcf") else
                        c("CASE_ID", "SAMPLE", "counts", "ascat_cnv", "ascat_pp")
missing <- setdiff(need, names(manifest))
if (length(missing)) {
  stop("Manifest is missing column(s): ", paste(missing, collapse = ", "))
}
if (!"sex" %in% names(manifest)) manifest[, sex := NA_character_]

## ------------------------------------------------------------- utilities ---

norm_chr <- function(x) sub("^chr", "", as.character(x), ignore.case = TRUE)

pick_col <- function(dt, patterns, what) {
  for (p in patterns) {
    hit <- grep(p, names(dt), ignore.case = TRUE, value = TRUE)
    if (length(hit)) return(hit[1])
  }
  stop("Could not find a column for '", what, "' in: ",
       paste(names(dt), collapse = ", "))
}

read_segments <- function(path) {
  seg <- fread(path)
  out <- seg[, .(CHR      = norm_chr(get(pick_col(seg, c("^chr$", "^chrom"), "chr"))),
                 start    = as.integer(get(pick_col(seg, c("^startpos$", "^start$"), "start"))),
                 end      = as.integer(get(pick_col(seg, c("^endpos$", "^end$"), "end"))),
                 major_cn = as.integer(get(pick_col(seg, c("^nmajor$", "major"), "nMajor"))),
                 minor_cn = as.integer(get(pick_col(seg, c("^nminor$", "minor"), "nMinor"))))]
  out[!is.na(start) & !is.na(end)]
}

read_purity <- function(path) {
  pp <- fread(path)
  if (ncol(pp) == 2L && nrow(pp) > 1L && is.character(pp[[1]]) &&
      !any(grepl("ploidy", names(pp), ignore.case = TRUE))) {
    acf <- as.numeric(pp[[2]])[grep("aberrant|purity|rho", tolower(pp[[1]]))][1]
  } else {
    acf <- as.numeric(pp[[pick_col(pp, c("aberrantcellfraction", "purity", "^rho$"),
                                   "ACF")]][1])
  }
  if (is.na(acf)) stop("Failed to parse purity from ", path)
  acf
}

read_counts <- function(path) {
  ct <- fread(path, header = FALSE,
              col.names = c("CHR", "POS", "REF", "ALT", "ref_counts", "var_counts"))
  ct[, `:=`(CHR = norm_chr(CHR), POS = as.integer(POS),
            ref_counts = as.integer(ref_counts),
            var_counts = as.integer(var_counts))]
  ct[!is.na(POS)]
}

normal_cn_for <- function(chr, sex) {
  ifelse(chr %in% c("X", "Y"),
         ifelse(is.na(sex) || sex == "XX", ifelse(chr == "X", 2L, 0L), 1L),
         2L)
}

keep_chrs_for <- function(sex) {
  if (is.na(sex)) KEEP_CHRS else c(KEEP_CHRS, "X", if (sex == "XY") "Y")
}

## =========================================================== STEPS 1-2 =======
## Union of PASS SNVs across the samples of each patient.

if (emit_sites) {
  if (nchar(Sys.which("bcftools")) == 0L) stop("bcftools not found on PATH")

  for (case in unique(manifest$CASE_ID)) {
    rows <- manifest[CASE_ID == case]
    message(case, ":")
    per <- lapply(seq_len(nrow(rows)), function(i) {
      v <- rows$pass_vcf[i]
      if (!file.exists(v)) stop("No such VCF: ", v)
      # split multiallelics so REF/ALT are unambiguous, then PASS SNVs only
      cmd <- sprintf(
        "bcftools norm -m -any %s 2>/dev/null | bcftools view -f PASS -v snps -H",
        shQuote(v))
      txt <- system(cmd, intern = TRUE)
      if (!length(txt)) {
        warning("No PASS SNVs in ", v,
                " -- was FilterMutectCalls run with --ob-priors?")
        return(NULL)
      }
      f <- tstrsplit(txt, "\t", fixed = TRUE)
      d <- data.table(CHR = norm_chr(f[[1]]), POS = as.integer(f[[2]]),
                      REF = f[[4]], ALT = f[[5]])
      message(sprintf("  %-24s %d PASS SNVs", rows$SAMPLE[i], nrow(d)))
      d
    })
    u <- unique(rbindlist(Filter(Negate(is.null), per)))
    u <- u[CHR %in% keep_chrs_for(rows$sex[1])]
    setorder(u, CHR, POS)

    f <- file.path(outdir, paste0(case, ".union_sites.tsv"))
    fwrite(u[, .(paste0("chr", CHR), POS, REF, ALT)], f,
           sep = "\t", col.names = FALSE, quote = FALSE)
    message(sprintf("  union: %d sites across %d sample(s) -> %s",
                    nrow(u), nrow(rows), f))
  }

  message("\nNow count alleles at these sites in EVERY sample, e.g.\n",
          "  cut -f1,2 CASE.union_sites.tsv > regions.txt\n",
          "  awk -v OFS='\\t' '{print $1,$2,$3\",\"$4}' CASE.union_sites.tsv \\\n",
          "    | bgzip > alleles.tsv.gz && tabix -s1 -b2 -e2 alleles.tsv.gz\n",
          "  bcftools mpileup -f $REF -R regions.txt -a AD -d 8000 -q 20 -Q 20 -Ou s.cram \\\n",
          "    | bcftools call -m -C alleles -T alleles.tsv.gz \\\n",
          "    | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%AD{0}\\t%AD{1}]\\n' > s.counts.tsv\n",
          "then re-run this script without --emit-sites.")
  quit(status = 0)
}

## =========================================================== STEP 4 ==========

message("Assembling PyClone input (", fmt, ") from ", nrow(manifest), " sample(s)")

build_sample <- function(row) {
  ct  <- read_counts(row$counts)
  seg <- read_segments(row$ascat_cnv)
  acf <- read_purity(row$ascat_pp)

  keep <- keep_chrs_for(row$sex)
  ct <- ct[CHR %in% keep]; seg <- seg[CHR %in% keep]
  if (!nrow(ct)) { warning("No counts for ", row$SAMPLE); return(NULL) }

  ct[, `:=`(mstart = POS, mend = POS)]
  setkey(seg, CHR, start, end)
  j <- foverlaps(ct, seg, by.x = c("CHR", "mstart", "mend"),
                 by.y = c("CHR", "start", "end"), type = "within", nomatch = NA)

  j[, `:=`(depth = ref_counts + var_counts, no_cn = is.na(major_cn))]
  # mutation_id must be IDENTICAL across samples -- never embed the sample name
  j[, mutation_id := paste(paste0("chr", CHR), POS, REF, ALT, sep = ":")]
  j[, `:=`(CASE_ID = row$CASE_ID, SAMPLE = row$SAMPLE,
           normal_cn = normal_cn_for(CHR, row$sex), tumour_content = acf)]

  message(sprintf("  %-24s %6d union sites counted  (ACF %.2f, clonal VAF ~%.2f)",
                  row$SAMPLE, nrow(j), acf, acf / 2))

  j[, .(CASE_ID, SAMPLE, mutation_id, CHR, POS, ref_counts, var_counts,
        depth, normal_cn, minor_cn, major_cn, no_cn, tumour_content)]
}

tbl <- rbindlist(Filter(Negate(is.null),
                        lapply(seq_len(nrow(manifest)),
                               function(i) build_sample(as.list(manifest[i])))))
if (!nrow(tbl)) stop("No rows produced - check inputs.")
n0 <- uniqueN(tbl$mutation_id)

## --- QC gate -----------------------------------------------------------------
# At purity p on diploid ground a clonal mutation sits at VAF ~ p/2, and that
# mode must be visible. A monotone decline from the first bin means the call
# set is artefact-dominated: fix the filtering rather than cluster it.

message("\nVAF distribution (union sites, depth >= ", MIN_DEPTH, "):")
for (s in unique(tbl$SAMPLE)) {
  v <- tbl[SAMPLE == s & depth >= MIN_DEPTH, var_counts / depth]
  if (!length(v)) next
  h <- table(cut(v, breaks = seq(0, 1, 0.05), include.lowest = TRUE))
  message(sprintf("  %s  (expect clonal mode near %.2f)",
                  s, tbl[SAMPLE == s, tumour_content[1]] / 2))
  for (i in seq_along(h)) {
    if (h[i] > 0)
      message(sprintf("    %-12s %-40s %d", names(h)[i],
                      strrep("#", min(40, as.integer(h[i]) %/% 5 + 1)),
                      as.integer(h[i])))
  }
}

## --- drop sites unusable for the model ---------------------------------------

bad_cn <- tbl[no_cn | (DROP_MAJOR0 & major_cn == 0L), unique(mutation_id)]
if (length(bad_cn)) {
  message(sprintf("\nDropping %d site(s): no ASCAT segment or major_cn == 0 in >=1 sample",
                  length(bad_cn)))
  tbl <- tbl[!mutation_id %in% bad_cn]
}

low <- tbl[depth < MIN_DEPTH, unique(mutation_id)]
if (length(low)) {
  message(sprintf("Dropping %d site(s): depth < %d in >=1 sample", length(low), MIN_DEPTH))
  tbl <- tbl[!mutation_id %in% low]
}

## --- completeness ------------------------------------------------------------
# After force-calling every union site should appear in every sample. Any that
# do not were skipped by the counter (usually zero coverage); drop them rather
# than hand PyClone a ragged matrix.

for (case in unique(tbl$CASE_ID)) {
  sub <- tbl[CASE_ID == case]
  ns  <- uniqueN(sub$SAMPLE)
  ragged <- sub[, .N, by = mutation_id][N < ns, mutation_id]
  if (length(ragged)) {
    message(sprintf("%s: %d site(s) absent from >=1 sample's counts - dropped.",
                    case, length(ragged)))
    if (length(ragged) > 0.05 * n0)
      message("      That is a large fraction: check the allele counter ran over\n",
              "      the full union region list in every sample.")
    tbl <- tbl[!(CASE_ID == case & mutation_id %in% ragged)]
  }
}
if (!nrow(tbl)) stop("Nothing left after filtering.")
message(sprintf("\n%d of %d union sites retained", uniqueN(tbl$mutation_id), n0))

## --- write -------------------------------------------------------------------

setorder(tbl, CASE_ID, SAMPLE, CHR, POS)

if (fmt == "vi") {
  for (case in unique(tbl$CASE_ID)) {
    out <- tbl[CASE_ID == case,
               .(mutation_id, sample_id = SAMPLE, ref_counts,
                 alt_counts = var_counts, major_cn, minor_cn, normal_cn,
                 tumour_content)]
    f <- file.path(outdir, paste0(case, ".pyclone_vi.tsv"))
    fwrite(out, f, sep = "\t", quote = FALSE)
    message("Wrote ", f, " (", uniqueN(out$mutation_id), " mutations x ",
            uniqueN(out$sample_id), " samples)")
  }
} else {
  for (case in unique(tbl$CASE_ID)) {
    sub <- tbl[CASE_ID == case]
    for (s in unique(sub$SAMPLE)) {
      out <- sub[SAMPLE == s,
                 .(mutation_id, ref_counts, var_counts, normal_cn, minor_cn,
                   major_cn, variant_case = SAMPLE,
                   variant_freq = fifelse(depth > 0, var_counts / depth, NA_real_),
                   genotype = "NA")]
      f <- file.path(outdir, paste0(case, ".", s, ".tsv"))
      fwrite(out, f, sep = "\t", quote = FALSE)
      message("Wrote ", f, " (", nrow(out), " mutations)")
    }
    tc <- unique(sub[, .(SAMPLE, tumour_content)]); setorder(tc, SAMPLE)
    message("  tumour contents: ",
            paste(sprintf("%s=%.3f", tc$SAMPLE, tc$tumour_content), collapse = " "))
  }
}