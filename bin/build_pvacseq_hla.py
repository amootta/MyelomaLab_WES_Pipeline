#!/usr/bin/env python3
"""
build_pvacseq_hla.py

Write the per-patient --hla_csv that nextflow-pvacseq needs, instead of typing
it by hand (HLA_pvacseq_COS_001_063.csv used to be written manually).

  * Sample_ID rows come from the tumour sample column of every VCF in
    <vcf-dir>/<patient>/ (stage 03 output), so they always match the VCFs.
  * HLA_Types come from one caller's class I (A, B, C) calls for that patient
    (stage 04 output, linked into <calls-dir>/<caller>/ by `wes.sh hla collect`),
    parsed with the same functions hla_concordance.py uses.

Output: <outdir>/<patient>.csv
    Sample_ID,HLA_Types
    COS_001_063_PRG_merged,HLA-A02:01;HLA-A25:01;HLA-B50:01;...

Usage:
    python3 build_pvacseq_hla.py --calls-dir HLA/calls --source optitype \
        --vcf-dir consensus/pvacseq_input --outdir pvacseq/hla
"""

import argparse
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hla_concordance import parse_hlahd_file, parse_optitype_file, parse_t1k_file  # noqa: E402

CLASS_I = ["A", "B", "C"]
FILE_GLOB = {
    "optitype": "{p}*_result.tsv",
    "hlahd": "{p}*_hlahd_final.txt",
    "t1k": "{p}*_HLA_genotype.tsv",
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--calls-dir", required=True, help="Folder with optitype/ hlahd/ t1k/ subfolders (wes.sh hla collect)")
    p.add_argument("--source", choices=sorted(FILE_GLOB), default="optitype", help="Which caller's typing to use")
    p.add_argument("--vcf-dir", required=True, help="Folder of <patient>/ subfolders holding that patient's tumour VCFs")
    p.add_argument("--outdir", required=True)
    p.add_argument("--resolution", type=int, default=2, help="Allele fields to keep (default 2, e.g. A*02:01)")
    return p.parse_args()


def tumour_names(patient_dir):
    names = []
    for vcf in sorted(glob.glob(os.path.join(patient_dir, "*.vcf"))):
        with open(vcf) as fh:
            for line in fh:
                if line.startswith("#CHROM"):
                    names.append(line.rstrip("\n").split("\t")[9])
                    break
    return names


def class_i_alleles(path, source, resolution):
    if source == "optitype":
        calls = parse_optitype_file(path, resolution)
    elif source == "hlahd":
        calls = parse_hlahd_file(path, resolution)
    else:
        calls, _ = parse_t1k_file(path, resolution, quality_threshold=0)
    alleles = []
    for gene in CLASS_I:
        for a in calls.get(gene, (None, None)):
            # "A*02:01" -> "HLA-A02:01": the format nextflow-pvacseq ran
            # successfully with (HLA_pvacseq_COS_001_063.csv), so kept as is.
            if a:
                fmt = "HLA-" + a.replace("*", "")
                if fmt not in alleles:
                    alleles.append(fmt)
    return alleles


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    calls_dir = os.path.join(args.calls_dir, args.source)
    written = 0
    for patient_dir in sorted(glob.glob(os.path.join(args.vcf_dir, "*"))):
        if not os.path.isdir(patient_dir):
            continue
        patient = os.path.basename(patient_dir)
        samples = tumour_names(patient_dir)
        hits = sorted(glob.glob(os.path.join(calls_dir, FILE_GLOB[args.source].format(p=patient))))
        if not samples:
            print(f"  [skip] {patient}: no VCFs in {patient_dir}", file=sys.stderr)
            continue
        if not hits:
            print(f"  [skip] {patient}: no {args.source} call file in {calls_dir}", file=sys.stderr)
            continue
        alleles = class_i_alleles(hits[0], args.source, args.resolution)
        if len(alleles) < 3:
            print(f"  [warn] {patient}: only {len(alleles)} class I alleles from {hits[0]}", file=sys.stderr)
        out = os.path.join(args.outdir, f"{patient}.csv")
        with open(out, "w", newline="\n") as fh:
            fh.write("Sample_ID,HLA_Types\n")
            for s in samples:
                fh.write(f"{s},{';'.join(alleles)}\n")
        print(f"  {patient}: {len(samples)} tumour sample(s), {len(alleles)} alleles -> {out}")
        written += 1
    if not written:
        sys.exit("No HLA csv written -- check --vcf-dir and --calls-dir")


if __name__ == "__main__":
    main()
