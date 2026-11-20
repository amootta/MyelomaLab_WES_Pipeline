#!/usr/bin/env python3
"""
build_lohhla_inputs.py

Lay out one LOHHLA input folder per patient from what the pipeline already has:

  <outdir>/<P>/input/bams/<S>.bam (+ .bam.bai)   links to the recal bams.
        LOHHLA (bitbucket original) takes every *.bam in --BAMDir as a sample,
        names each one after its file name minus ".bam", and keys read counts on
        the name up to the first dot -- so "<S>.recal.bam" would break it. The
        normal is the one whose path equals --normalBAMfile; the rest are tumours.
  <outdir>/<P>/input/hla_alleles.txt              one allele per line, e.g.
        hla_a_02_01_01_01, picked from --hla-fasta for each germline call.
  <outdir>/<P>/input/solutions.txt                purity/ploidy per tumour:
        tumorPurity<TAB>tumorPloidy          (header one field short, so R's
        <T><TAB>0.62<TAB>2.94                 read.table uses column 1 as row names)
        Tumours without an ASCAT result are left out; LOHHLA then reports
        allelic imbalance for them but no integer copy number.

The manifest (written by stages/11_lohhla/run.sh) is tab-separated, no header:
    patient  sample  role(normal|tumour)  bam  ascat_purityploidy_or_empty
"""

import argparse
import glob
import os
import re
import sys
from collections import OrderedDict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "bin"))
from hla_concordance import parse_hlahd_file, parse_optitype_file, parse_t1k_file  # noqa: E402

CLASS_I = ["A", "B", "C"]
FILE_GLOB = {
    "optitype": "{p}*_result.tsv",
    "hlahd": "{p}*_hlahd_final.txt",
    "t1k": "{p}*_HLA_genotype.tsv",
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--manifest", required=True)
    p.add_argument("--calls-dir", required=True, help="HLA/calls (wes.sh hla collect)")
    p.add_argument("--source", choices=sorted(FILE_GLOB), default="optitype")
    p.add_argument("--hla-fasta", required=True, help="fasta named hla_a_01_01_01_01 etc.")
    p.add_argument("--outdir", required=True)
    p.add_argument("--resolution", type=int, default=2, help="fields of the call to match on (default 2)")
    return p.parse_args()


def fasta_names(path):
    with open(path) as fh:
        return [line[1:].split()[0] for line in fh if line.startswith(">")]


def field_key(name):
    # hla_a_02_01_01_01 -> (2, 1, 1, 1); expression suffix (N/L/Q...) sorts last
    parts = name.split("_")[2:]
    return (tuple(int(re.sub(r"\D", "", x) or 0) for x in parts), bool(re.search(r"[a-z]$", parts[-1])))


def pick_allele(call, names):
    """'A*02:01' -> the lowest-numbered fasta allele in that group
    (hla_a_02_01_01_01), preferring expressed alleles over N/L/Q ones."""
    gene, fields = call.split("*", 1)
    prefix = "hla_" + gene.lower() + "_" + fields.replace(":", "_").lower()
    hits = [n for n in names if n == prefix or n.startswith(prefix + "_")]
    if not hits:
        return None
    return sorted(hits, key=lambda n: (field_key(n)[1], field_key(n)[0]))[0]


def class_i_calls(path, source, resolution):
    if source == "optitype":
        return parse_optitype_file(path, resolution)
    if source == "hlahd":
        return parse_hlahd_file(path, resolution)
    return parse_t1k_file(path, resolution, quality_threshold=0)[0]


def read_purity_ploidy(path):
    """sarek's ASCAT <pair>.purityploidy.txt: a header row naming purity
    (AberrantCellFraction) and ploidy, one value row; also accepts key/value rows."""
    rows = [re.split(r"[\t ,]+", l.strip().replace('"', "")) for l in open(path) if l.strip()]
    vals = {}
    if len(rows) >= 2 and any(re.search("ploidy", h, re.I) for h in rows[0]):
        head, data = rows[0], rows[1]
        if len(data) == len(head) + 1:   # R row name in front
            data = data[1:]
        vals = dict(zip((h.lower() for h in head), data))
    else:
        vals = {r[0].lower(): r[-1] for r in rows if len(r) >= 2}
    purity = next((v for k, v in vals.items() if re.search("aberrant|purity|rho", k)), None)
    ploidy = next((v for k, v in vals.items() if re.search("ploidy|psi", k)), None)
    try:
        return float(purity), float(ploidy)
    except (TypeError, ValueError):
        return None


def link(src, dst):
    if os.path.lexists(dst):
        os.remove(dst)
    os.symlink(os.path.abspath(src), dst)


def bam_index(bam):
    for c in (bam + ".bai", re.sub(r"\.bam$", ".bai", bam)):
        if os.path.exists(c):
            return c
    return None


def main():
    args = parse_args()
    names = fasta_names(args.hla_fasta)
    if not names:
        sys.exit(f"no sequences in {args.hla_fasta}")

    patients = OrderedDict()
    with open(args.manifest) as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t") + [""] * 5
            patients.setdefault(f[0], []).append({"sample": f[1], "role": f[2], "bam": f[3], "pp": f[4]})

    ready = 0
    for p, rows in patients.items():
        normals = [r for r in rows if r["role"] == "normal"]
        tumours = [r for r in rows if r["role"] == "tumour"]
        if len(normals) != 1 or not tumours:
            print(f"  [skip] {p}: needs one normal and >=1 tumour", file=sys.stderr)
            continue
        hits = sorted(glob.glob(os.path.join(args.calls_dir, args.source, FILE_GLOB[args.source].format(p=p))))
        if not hits:
            print(f"  [skip] {p}: no {args.source} calls in {args.calls_dir}/{args.source}", file=sys.stderr)
            continue

        indir = os.path.join(args.outdir, p, "input")
        bamdir = os.path.join(indir, "bams")
        os.makedirs(bamdir, exist_ok=True)

        missing = False
        for r in rows:
            idx = bam_index(r["bam"])
            if not os.path.exists(r["bam"]) or not idx:
                print(f"  [skip] {p}: bam or index missing for {r['sample']}: {r['bam']}", file=sys.stderr)
                missing = True
                break
            link(r["bam"], os.path.join(bamdir, r["sample"] + ".bam"))
            link(idx, os.path.join(bamdir, r["sample"] + ".bam.bai"))
        if missing:
            continue

        # HLA alleles
        calls = class_i_calls(hits[0], args.source, args.resolution)
        alleles, notes = [], []
        for gene in CLASS_I:
            pair = [a for a in calls.get(gene, (None, None)) if a]
            picked = []
            for a in pair:
                n = pick_allele(a, names)
                if n is None:
                    notes.append(f"{a} not in fasta")
                elif n not in picked:
                    picked.append(n)
            if len(picked) < 2:
                notes.append(f"HLA-{gene} homozygous/incomplete ({', '.join(pair) or 'no call'}): LOHHLA skips it")
            alleles += picked
        with open(os.path.join(indir, "hla_alleles.txt"), "w", newline="\n") as fh:
            fh.write("\n".join(alleles) + "\n")

        # purity / ploidy
        sol = []
        for t in tumours:
            pp = read_purity_ploidy(t["pp"]) if t["pp"] and os.path.exists(t["pp"]) else None
            if pp:
                sol.append(f"{t['sample']}\t{pp[0]:g}\t{pp[1]:g}")
            else:
                notes.append(f"{t['sample']}: no ASCAT purity/ploidy, allelic imbalance only")
        with open(os.path.join(indir, "solutions.txt"), "w", newline="\n") as fh:
            fh.write("tumorPurity\ttumorPloidy\n" + "".join(s + "\n" for s in sol))
        with open(os.path.join(indir, "normal.txt"), "w", newline="\n") as fh:
            fh.write(normals[0]["sample"] + "\n")

        print(f"  {p}: normal {normals[0]['sample']}, {len(tumours)} tumour(s), "
              f"{len(alleles)} alleles from {os.path.basename(hits[0])}")
        for n in notes:
            print(f"      note: {n}")
        ready += 1

    if not ready:
        sys.exit("no patient inputs written")


if __name__ == "__main__":
    main()
