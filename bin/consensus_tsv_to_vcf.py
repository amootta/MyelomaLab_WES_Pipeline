#!/usr/bin/env python3
"""
consensus_tsv_to_vcf.py

Convert an mgp1000 union/consensus somatic SNV TSV (post filter_consensus_calls.py,
i.e. already carrying n_callers_pass / callers_pass) into a minimal, valid VCFv4.2
file suitable for sorting/bgzip/tabix -> VEP annotation -> pVACtools (pvacseq).

Why this exists: pVACtools requires a VEP-annotated VCF as input, not a flat table.
This script produces the *pre-VEP* VCF. It does NOT run VEP, bgzip, tabix, or
pvacseq itself -- see the accompanying shell commands for those steps.

Input assumptions (matches Amatta's real 117+2-column header exactly):
  seqnames, start, end, strand, width, CALLER, TUMOR, NORMAL, SAMPLE, PATIENT,
  REF, ALT, n_callers_pass, callers_pass, then MUTECT_*, STRELKA_*, VARSCAN_*
  blocks, then nearest_gene, alt_read_depth_combo, alt_read_depth_mean,
  total_depth_combo, total_depth_mean, vaf_combo, vaf_mean.
  Missing values are the literal string "NA". SNVs only (single-base REF/ALT).

Genotype/AD/DP source per row: picks the first caller (in --caller-priority
order, default mutect,varscan,strelka) that is listed in that row's
callers_pass, and pulls that caller's own TUMOR/NORMAL AD/DP fields. This
keeps FORMAT values internally consistent (AD+other_AD == DP from one caller's
own accounting) rather than blending stats across callers with different
counting conventions.

Usage (single file):
    python3 consensus_tsv_to_vcf.py \
        --input COS_001_018_BL_vs_COS_001_018_G.hq.union.consensus.somatic.snv.callers_ge2.txt.gz \
        --out COS_001_018_BL_vs_COS_001_018_G.consensus.vcf

Usage (cohort / glob mode, one output VCF per input file, mirrors
filter_consensus_calls.py's cohort mode):
    python3 consensus_tsv_to_vcf.py \
        --input "*.callers_ge2.txt.gz" \
        --outdir vcfs/

Optional:
    --reference-fai /path/to/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.fai
        Adds proper ##contig lines (recommended before bcftools sort / VEP) and
        fixes chromosome output order to match the reference. Without it,
        contigs are emitted in first-seen order and ##contig lines are omitted
        (VEP/bcftools will generally still work, just less strict).
    --min-vaf-for-hom-alt FLOAT (default 0.9)
        Tumor VAF at/above this is called 1/1 instead of 0/1.
    --caller-priority mutect,varscan,strelka
        Order to try callers in for FORMAT fields.
"""

import argparse
import glob
import gzip
import os
import sys
from datetime import date

import pandas as pd

NA = "NA"

BASE_TO_STRELKA_COL = {"A": "AU", "C": "CU", "G": "GU", "T": "TU"}

REQUIRED_COLS = [
    "seqnames", "start", "end", "REF", "ALT", "TUMOR", "NORMAL",
    "SAMPLE", "PATIENT", "n_callers_pass", "callers_pass", "nearest_gene",
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, help="Input TSV(.gz) path, or a quoted glob for cohort mode")
    p.add_argument("--out", help="Output VCF path (single-file mode)")
    p.add_argument("--outdir", help="Output directory (cohort/glob mode) -- one <input_stem>.vcf per input file")
    p.add_argument("--reference-fai", help="Path to the reference .fai (e.g. GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.fai) for ##contig lines and chrom ordering")
    p.add_argument("--min-vaf-for-hom-alt", type=float, default=0.9)
    p.add_argument("--caller-priority", default="mutect,varscan,strelka")
    return p.parse_args()


def read_fai(path):
    """Return list of (chrom, length) in file order."""
    contigs = []
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                contigs.append((parts[0], int(parts[1])))
    return contigs


def na_to_none(v):
    if v is None:
        return None
    v = v.strip()
    if v == "" or v == NA or v.upper() == "NA":
        return None
    return v


def get_mutect_ad_dp(row, group):
    ref = na_to_none(row.get(f"MUTECT_AD_REF_{group}"))
    alt = na_to_none(row.get(f"MUTECT_AD_ALT_{group}"))
    dp = na_to_none(row.get(f"MUTECT_DP_{group}"))
    gt = na_to_none(row.get(f"MUTECT_GT_{group}"))
    if ref is None or alt is None:
        return None
    try:
        ref_i, alt_i = int(float(ref)), int(float(alt))
    except ValueError:
        return None
    dp_i = int(float(dp)) if dp is not None else ref_i + alt_i
    return ref_i, alt_i, dp_i, gt


def get_varscan_ad_dp(row, group):
    rd = na_to_none(row.get(f"VARSCAN_RD_{group}"))
    ad = na_to_none(row.get(f"VARSCAN_AD_{group}"))
    dp = na_to_none(row.get(f"VARSCAN_DP_{group}"))
    gt = na_to_none(row.get(f"VARSCAN_GT_{group}"))
    if rd is None or ad is None:
        return None
    try:
        ref_i, alt_i = int(float(rd)), int(float(ad))
    except ValueError:
        return None
    dp_i = int(float(dp)) if dp is not None else ref_i + alt_i
    return ref_i, alt_i, dp_i, gt


def get_strelka_ad_dp(row, group, ref_base, alt_base):
    ref_col = BASE_TO_STRELKA_COL.get(ref_base.upper())
    alt_col = BASE_TO_STRELKA_COL.get(alt_base.upper())
    if ref_col is None or alt_col is None:
        return None  # not a simple single-base SNV
    ref_v = na_to_none(row.get(f"STRELKA_{ref_col}_TIER1_{group}"))
    alt_v = na_to_none(row.get(f"STRELKA_{alt_col}_TIER1_{group}"))
    dp = na_to_none(row.get(f"STRELKA_DP_{group}"))
    if ref_v is None or alt_v is None:
        return None
    try:
        ref_i, alt_i = int(float(ref_v)), int(float(alt_v))
    except ValueError:
        return None
    dp_i = int(float(dp)) if dp is not None else ref_i + alt_i
    return ref_i, alt_i, dp_i, None  # strelka has no NORMAL/TUMOR GT field here


GETTERS = {
    "mutect": lambda row, group, ref_b, alt_b: get_mutect_ad_dp(row, group),
    "varscan": lambda row, group, ref_b, alt_b: get_varscan_ad_dp(row, group),
    "strelka": lambda row, group, ref_b, alt_b: get_strelka_ad_dp(row, group, ref_b, alt_b),
}


def build_format_sample(ref_i, alt_i, dp_i, gt, is_tumor, min_vaf_hom):
    vaf = alt_i / dp_i if dp_i else 0.0
    if gt is not None:
        gt_out = gt.replace("|", "/")
    else:
        if is_tumor:
            gt_out = "1/1" if vaf >= min_vaf_hom else "0/1"
        else:
            gt_out = "0/0"
    return f"{gt_out}:{ref_i},{alt_i}:{dp_i}:{vaf:.4f}"


def convert_file(in_path, out_path, contigs, caller_priority, min_vaf_hom):
    df = pd.read_csv(in_path, sep="\t", dtype=str, keep_default_na=False, na_filter=False)

    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        sys.exit(f"[{in_path}] missing expected columns: {missing}")

    tumor_names = df["TUMOR"].unique()
    normal_names = df["NORMAL"].unique()
    if len(tumor_names) != 1 or len(normal_names) != 1:
        print(f"WARNING [{in_path}]: multiple distinct TUMOR/NORMAL values found "
              f"(TUMOR={list(tumor_names)}, NORMAL={list(normal_names)}); "
              f"using the first row's values as the VCF sample names.", file=sys.stderr)
    tumor_sample = df["TUMOR"].iloc[0]
    normal_sample = df["NORMAL"].iloc[0]

    contig_order = {c: i for i, (c, _) in enumerate(contigs)} if contigs else None

    records = []
    n_skipped_no_caller_data = 0

    for _, row in df.iterrows():
        chrom = row["seqnames"]
        pos = row["start"]
        ref = row["REF"]
        alt = row["ALT"]

        if len(ref) != 1 or len(alt) != 1:
            # This script only handles simple SNVs (matches the .snv. file); skip
            # anything else (shouldn't occur in an snv-only union file, but guard anyway).
            n_skipped_no_caller_data += 1
            continue

        callers_pass = [c.strip().lower() for c in row["callers_pass"].split(",") if c.strip()]

        chosen = None
        tumor_vals = normal_vals = None
        for caller in caller_priority:
            if caller not in callers_pass:
                continue
            t = GETTERS[caller](row, "TUMOR", ref, alt)
            if t is None:
                continue
            n = GETTERS[caller](row, "NORMAL", ref, alt)
            chosen = caller
            tumor_vals = t
            normal_vals = n
            break

        if chosen is None:
            n_skipped_no_caller_data += 1
            continue

        t_ref, t_alt, t_dp, t_gt = tumor_vals
        tumor_fmt = build_format_sample(t_ref, t_alt, t_dp, t_gt, is_tumor=True, min_vaf_hom=min_vaf_hom)

        if normal_vals is not None:
            n_ref, n_alt, n_dp, n_gt = normal_vals
            normal_fmt = build_format_sample(n_ref, n_alt, n_dp, n_gt, is_tumor=False, min_vaf_hom=min_vaf_hom)
        else:
            normal_fmt = "0/0:.,.:.:."

        gene = row["nearest_gene"] if row["nearest_gene"] not in ("", NA) else "."
        info = (f"NCALLERS={row['n_callers_pass']};"
                f"CALLERS={row['callers_pass']};"
                f"SRC_CALLER={chosen};"
                f"GENE={gene}")

        sort_key = (contig_order.get(chrom, len(contig_order) if contig_order else 0) if contig_order else 0,
                    chrom, int(pos))
        records.append((sort_key, chrom, pos, ref, alt, info, tumor_fmt, normal_fmt))

    records.sort(key=lambda r: (r[0][0], r[0][1], r[0][2]))

    with open(out_path, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write(f"##fileDate={date.today().strftime('%Y%m%d')}\n")
        out.write(f"##source=consensus_tsv_to_vcf.py (from {os.path.basename(in_path)})\n")
        if contigs:
            for c, length in contigs:
                out.write(f"##contig=<ID={c},length={length}>\n")
        out.write('##INFO=<ID=NCALLERS,Number=1,Type=Integer,Description="Number of callers with FILTER=PASS supporting this variant">\n')
        out.write('##INFO=<ID=CALLERS,Number=1,Type=String,Description="Comma-separated list of callers with FILTER=PASS">\n')
        out.write('##INFO=<ID=SRC_CALLER,Number=1,Type=String,Description="Caller whose AD/DP/GT values were used to populate FORMAT fields for this record">\n')
        out.write('##INFO=<ID=GENE,Number=1,Type=String,Description="Nearest gene from the mgp1000 consensus TSV">\n')
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        out.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles">\n')
        out.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n')
        out.write('##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele fraction (alt/DP) from SRC_CALLER">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + tumor_sample + "\t" + normal_sample + "\n")
        for _, chrom, pos, ref, alt, info, tumor_fmt, normal_fmt in records:
            out.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info}\tGT:AD:DP:AF\t{tumor_fmt}\t{normal_fmt}\n")

    n_total = len(df)
    n_written = len(records)
    print(f"[{in_path}] {n_written}/{n_total} rows written to {out_path} "
          f"({n_skipped_no_caller_data} skipped: non-SNV or no usable caller AD/DP data)")
    print(f"    samples: TUMOR={tumor_sample}  NORMAL={normal_sample}")


def main():
    args = parse_args()
    caller_priority = [c.strip().lower() for c in args.caller_priority.split(",") if c.strip()]
    for c in caller_priority:
        if c not in GETTERS:
            sys.exit(f"Unknown caller in --caller-priority: {c} (must be one of {list(GETTERS)})")

    contigs = read_fai(args.reference_fai) if args.reference_fai else []

    is_glob = any(ch in args.input for ch in "*?[")
    if is_glob:
        if not args.outdir:
            sys.exit("--outdir is required in cohort/glob mode")
        os.makedirs(args.outdir, exist_ok=True)
        paths = sorted(glob.glob(args.input))
        if not paths:
            sys.exit(f"No files matched glob: {args.input}")
        for p in paths:
            stem = os.path.basename(p)
            for suf in (".txt.gz", ".tsv.gz", ".txt", ".tsv", ".gz"):
                if stem.endswith(suf):
                    stem = stem[: -len(suf)]
                    break
            out_path = os.path.join(args.outdir, stem + ".vcf")
            convert_file(p, out_path, contigs, caller_priority, args.min_vaf_for_hom_alt)
    else:
        if not args.out:
            sys.exit("--out is required in single-file mode")
        convert_file(args.input, args.out, contigs, caller_priority, args.min_vaf_for_hom_alt)


if __name__ == "__main__":
    main()