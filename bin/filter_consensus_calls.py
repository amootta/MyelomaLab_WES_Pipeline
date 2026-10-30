#!/usr/bin/env python3
"""
filter_consensus_calls.py

Filter mgp1000's per-pair union/consensus somatic SNV TSV
(*.hq.union.consensus.somatic.snv.txt.gz, produced by
snvindel_union_consensus_polisher.R) down to variants called
("FILTER" == PASS) by at least N of the three callers whose columns
the file already carries (MUTECT_*, STRELKA_*, VARSCAN_*).

The consensus file is already a per-position union across callers:
each row carries every caller's own column block, NA-filled for
callers that didn't report that position. This script does not
re-merge anything -- it just counts, per row, how many of
MUTECT_FILTER / STRELKA_FILTER / VARSCAN_FILTER equal "PASS", adds
two columns (n_callers_pass, callers_pass e.g. "mutect,strelka"),
and keeps rows meeting --min-callers.

Usage (single file):
    python3 filter_consensus_calls.py \
        --input COS_001_018_BL_vs_COS_001_018_G.hq.union.consensus.somatic.snv.txt.gz \
        --min-callers 2 \
        --out COS_001_018_BL_vs_COS_001_018_G.callers_ge2.tsv.gz

Usage (whole cohort, one output per input, same directory naming pattern):
    python3 filter_consensus_calls.py \
        --input "/path/to/*.hq.union.consensus.somatic.snv.txt.gz" \
        --min-callers 2 \
        --outdir filtered/

Notes:
  * "Called" = that caller's own FILTER column equals PASS
    (case-insensitive). This is the strict definition -- a caller
    that reported the position but flagged/filtered it does NOT
    count as calling it. Missing (NA) never counts.
  * Add/remove callers by editing CALLER_FILTER_COLS below if a
    cohort's consensus file gains/drops a caller block.
  * Every other column is passed through untouched (read as string,
    NA tokens preserved as literal "NA" on output) so numeric
    formatting in the original file is never altered.
"""

import argparse
import glob
import os
import sys

import pandas as pd

CALLER_FILTER_COLS = {
    "mutect": "MUTECT_FILTER",
    "strelka": "STRELKA_FILTER",
    "varscan": "VARSCAN_FILTER",
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True,
                   help="Path to one consensus TSV(.gz), or a glob pattern matching several "
                        "(quote the pattern so the shell doesn't expand it, e.g. \"*.consensus.somatic.snv.txt.gz\")")
    p.add_argument("--min-callers", type=int, default=2,
                   help="Minimum number of callers with FILTER==PASS to keep a variant (default: 2, i.e. 2/3 or 3/3)")
    p.add_argument("--out", default=None,
                   help="Output path (single-file mode only). If omitted, derived from --input by inserting "
                        "'.callers_ge<N>' before the extension.")
    p.add_argument("--outdir", default=None,
                   help="Output directory (multi-file / glob mode). Each input gets the same derived filename "
                        "written into this directory. Ignored in single-file mode if --out is given.")
    p.add_argument("--keep-below-threshold-summary", action="store_true",
                   help="Also print a breakdown of variants that did NOT meet --min-callers, by caller count")
    return p.parse_args()


def derive_out_path(in_path, min_callers, outdir=None):
    base = os.path.basename(in_path)
    # strip a trailing .gz first, then the extension, so
    # foo.hq.union.consensus.somatic.snv.txt.gz -> foo.hq.union.consensus.somatic.snv.callers_ge2.txt.gz
    gz = base.endswith(".gz")
    stem = base[:-3] if gz else base
    root, ext = os.path.splitext(stem)
    new_name = f"{root}.callers_ge{min_callers}{ext}{'.gz' if gz else ''}"
    if outdir:
        os.makedirs(outdir, exist_ok=True)
        return os.path.join(outdir, new_name)
    return os.path.join(os.path.dirname(in_path) or ".", new_name)


def process_file(in_path, min_callers, out_path, show_below_summary=False):
    missing = []
    # Read everything as string so we don't touch numeric formatting;
    # pandas still converts recognized NA tokens (incl. literal "NA") to
    # real NaN even under dtype=str, which is what we need for the PASS check.
    df = pd.read_csv(in_path, sep="\t", dtype=str, compression="infer", low_memory=False)

    present_cols = {}
    for caller, col in CALLER_FILTER_COLS.items():
        if col in df.columns:
            present_cols[caller] = col
        else:
            missing.append(col)
    if missing:
        print(f"  WARNING: {in_path} is missing expected column(s) {missing} -- "
              f"proceeding with the {len(present_cols)} caller(s) that ARE present.",
              file=sys.stderr)
    if not present_cols:
        print(f"  ERROR: none of {list(CALLER_FILTER_COLS.values())} found in {in_path}; skipping.", file=sys.stderr)
        return

    pass_mask = {}
    for caller, col in present_cols.items():
        pass_mask[caller] = df[col].str.strip().str.upper() == "PASS"

    mask_df = pd.DataFrame(pass_mask)
    df["n_callers_pass"] = mask_df.sum(axis=1)
    df["callers_pass"] = mask_df.apply(lambda row: ",".join(c for c in present_cols if row[c]), axis=1)

    # place the two new columns right after ALT for visibility, if ALT exists
    cols = list(df.columns)
    for new_col in ("callers_pass", "n_callers_pass"):
        cols.remove(new_col)
    if "ALT" in cols:
        insert_at = cols.index("ALT") + 1
    else:
        insert_at = len(cols)
    cols[insert_at:insert_at] = ["n_callers_pass", "callers_pass"]
    df = df[cols]

    kept = df[df["n_callers_pass"] >= min_callers].copy()

    print(f"{in_path}")
    print(f"  total variants:        {len(df)}")
    print(f"  callers present:       {', '.join(present_cols)}")
    print(f"  kept (>= {min_callers} callers): {len(kept)}")
    vc = df["n_callers_pass"].value_counts().sort_index(ascending=False)
    for n, count in vc.items():
        print(f"    n_callers_pass={n}: {count}")
    if show_below_summary:
        below = df[df["n_callers_pass"] < min_callers]
        if len(below):
            combo_counts = below.loc[below["n_callers_pass"] > 0, "callers_pass"].value_counts()
            print("  below-threshold caller combos (excluding 0-caller rows):")
            for combo, count in combo_counts.items():
                print(f"    {combo}: {count}")

    kept.to_csv(out_path, sep="\t", index=False, na_rep="NA",
                compression="gzip" if out_path.endswith(".gz") else None)
    print(f"  wrote: {out_path}\n")


def main():
    args = parse_args()
    matches = sorted(glob.glob(args.input))
    if not matches:
        # not a glob -- treat as a literal path so a clear file-not-found surfaces
        matches = [args.input]

    if len(matches) > 1 and args.out:
        print("ERROR: --out is only valid with a single input file; use --outdir for multiple inputs.", file=sys.stderr)
        sys.exit(1)

    for in_path in matches:
        if not os.path.exists(in_path):
            print(f"ERROR: input not found: {in_path}", file=sys.stderr)
            continue
        out_path = args.out if (args.out and len(matches) == 1) else derive_out_path(in_path, args.min_callers, args.outdir)
        process_file(in_path, args.min_callers, out_path, args.keep_below_threshold_summary)


if __name__ == "__main__":
    main()