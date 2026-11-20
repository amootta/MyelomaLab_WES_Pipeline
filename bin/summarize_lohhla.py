#!/usr/bin/env python3
"""
summarize_lohhla.py

Collect every LOHHLA prediction table under <lohhla-dir>/<patient>/
(<P>.<minCov>.DNA.HLAlossPrediction_CI.xls -- tab-separated despite the name),
drop the duplicate rows LOHHLA writes, and call each tumour x HLA gene:

  LOH            one allele's copy number < --loss-cn (default 0.5)
                 and PVal_unique < --pval (default 0.01)
  imbalance      PVal_unique < --pval, but neither allele below --loss-cn
                 (or no copy number: no ASCAT purity/ploidy for that tumour)
  no_change      PVal_unique >= --pval
  insufficient   no p-value (too few mismatch sites), or fewer than
                 --min-sites covered mismatch sites

The copy numbers used are HLA_type{1,2}copyNum_withBAFBin, as in
McGranahan et al. 2017 (Cell). Homozygous genes are not in LOHHLA's output, so
they are not listed here.

Writes:
  <out>.tsv       one row per tumour x gene, with the call and the numbers behind it
  <out>_wide.tsv  one row per tumour, one column per gene (HLA-A, HLA-B, HLA-C)

Only the standard library is used.
"""

import argparse
import csv
import glob
import math
import os
import re
import sys


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--lohhla-dir", required=True, help="LOHHLA output root (one folder per patient)")
    p.add_argument("--out", required=True, help="output prefix, e.g. LOHHLA/LOHHLA_summary")
    p.add_argument("--loss-cn", type=float, default=0.5)
    p.add_argument("--pval", type=float, default=0.01)
    p.add_argument("--min-sites", type=int, default=0, help="minimum numMisMatchSitesCov (default 0 = no filter)")
    return p.parse_args()


def num(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(v) else v


def col(row, *names):
    """First matching column, ignoring case (LOHHLA versions differ: PVal_unique / Pval_unique)."""
    low = {k.lower(): v for k, v in row.items() if k}
    for n in names:
        if n.lower() in low:
            return low[n.lower()]
    return None


def gene_of(allele):
    m = re.match(r"hla_([a-z]+)_", allele or "", re.I)
    return "HLA-" + m.group(1).upper() if m else "?"


def pretty(allele):
    # hla_a_02_01_01_01 -> A*02:01:01:01
    parts = (allele or "").split("_")
    if len(parts) < 3:
        return allele or ""
    return parts[1].upper() + "*" + ":".join(parts[2:])


def call(cn1, cn2, p, sites, a):
    if p is None or (a.min_sites and (sites is None or sites < a.min_sites)):
        return "insufficient"
    if p >= a.pval:
        return "no_change"
    if cn1 is not None and cn2 is not None and min(cn1, cn2) < a.loss_cn:
        return "LOH"
    return "imbalance"


def main():
    a = parse_args()
    files = sorted(glob.glob(os.path.join(a.lohhla_dir, "*", "*.DNA.HLAlossPrediction_CI.*")))
    if not files:
        sys.exit(f"no *.DNA.HLAlossPrediction_CI.* files under {a.lohhla_dir}/*/")

    rows, seen = [], set()
    for f in files:
        patient = os.path.basename(os.path.dirname(f))
        with open(f, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                region = col(r, "region") or ""
                # LOHHLA names the allele columns HLA_A_type1/2 for every gene
                t1, t2 = col(r, "HLA_A_type1"), col(r, "HLA_A_type2")
                key = (patient, region, t1, t2)
                if key in seen:
                    continue
                seen.add(key)
                cn1 = num(col(r, "HLA_type1copyNum_withBAFBin"))
                cn2 = num(col(r, "HLA_type2copyNum_withBAFBin"))
                p = num(col(r, "PVal_unique"))
                sites = num(col(r, "numMisMatchSitesCov"))
                c = call(cn1, cn2, p, sites, a)
                lost = col(r, "LossAllele") if c == "LOH" else ""
                rows.append({
                    "patient": patient,
                    "tumour": os.path.basename(region),
                    "gene": gene_of(t1),
                    "allele1": pretty(t1),
                    "allele2": pretty(t2),
                    "call": c,
                    "lost_allele": pretty(lost) if lost else "",
                    "cn_allele1": "" if cn1 is None else f"{cn1:.3g}",
                    "cn_allele2": "" if cn2 is None else f"{cn2:.3g}",
                    "pval_unique": "" if p is None else f"{p:.3g}",
                    "mismatch_sites": "" if sites is None else f"{sites:g}",
                    "source": os.path.relpath(f, a.lohhla_dir),
                })

    rows.sort(key=lambda r: (r["patient"], r["tumour"], r["gene"]))
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out + ".tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    genes = ["HLA-A", "HLA-B", "HLA-C"]
    wide = {}
    for r in rows:
        w = wide.setdefault((r["patient"], r["tumour"]), {g: "homozygous/not run" for g in genes})
        w[r["gene"]] = r["call"] + (f" ({r['lost_allele']})" if r["lost_allele"] else "")
    with open(a.out + "_wide.tsv", "w", newline="") as fh:
        fh.write("\t".join(["patient", "tumour"] + genes) + "\n")
        for (p, t), w in sorted(wide.items()):
            fh.write("\t".join([p, t] + [w[g] for g in genes]) + "\n")

    counts = {}
    for r in rows:
        counts[r["call"]] = counts.get(r["call"], 0) + 1
    print(f"{len(files)} file(s), {len(rows)} tumour x gene rows: "
          + ", ".join(f"{k} {v}" for k, v in sorted(counts.items())))
    print(f"wrote {a.out}.tsv and {a.out}_wide.tsv")


if __name__ == "__main__":
    main()
