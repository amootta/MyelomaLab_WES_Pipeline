#!/usr/bin/env python3
r"""
kir_concordance.py

Compare T1K KIR genotype TSVs (one file per sample) against a "true" (e.g.
SSO/array-based) KIR typing table, across as many samples as you have T1K
output for. Produces one Excel report with a long-format comparison table,
a wide matrix that mirrors the layout of the true-typing table, and a
per-sample concordance summary.

USAGE
-----
    python3 kir_concordance.py --t1k-dir /path/to/t1k_outputs --out KIR_concordance.xlsx

    # custom filename -> sample ID pattern (default matches cos_004_003_KIR_genotype.tsv -> "004-003")
    python3 kir_concordance.py --t1k-dir . --out out.xlsx \
        --filename-pattern 'cos_(\d+)_(\d+)_KIR_genotype\.tsv'

    # use your own ground-truth table instead of the 6 embedded samples
    python3 kir_concordance.py --t1k-dir . --out out.xlsx --truth-csv truth.csv

    # loosen/tighten what counts as a confident T1K call (default: quality >= 1)
    python3 kir_concordance.py --t1k-dir . --out out.xlsx --quality-threshold 1

INPUT FORMAT ASSUMED FOR EACH T1K FILE (tab-separated, no header)
-------------------------------------------------------------------
    gene  n_alleles  allele1  abundance1  quality1  allele2  abundance2  quality2  [extra candidate columns...]

Any columns beyond the first 8 are treated as secondary/ambiguous allele
candidates (T1K sometimes appends "allele;abundance;quality" style extra
fields, e.g. for KIR3DL3) and are kept for reference only -- they don't
feed into the presence/absence call.

WHAT COUNTS AS "PRESENT"
-------------------------
A gene is called present ("+") if T1K reports at least one allele AND the
max quality across its 1-2 alleles is >= --quality-threshold (default 1).
This threshold was chosen because several loci (2DS2, 2DS3, 2DS5, 3DS1,
2DL5A/2DL5B) can get a low-abundance, quality=0 allele call that looks like
noise/off-target mapping rather than a real detection -- see the writeup
this script grew out of. If a gene clears n_alleles>=1 but not the quality
bar, it's marked "- (low-conf hit)" so you can still see what T1K found.

KIR2DS4 full-length vs deletion
--------------------------------
KIR2DS4 has a well-described 22bp exon-5 deletion that abolishes surface
expression, and typing panels often report DEL/FULL rather than +/- for
this gene. The allele -> group mapping below is taken from the GenDx KIR
Genotyping SSO kit interpretation table (2025, 2V2):
  https://labproducts.gendx.com/productfiles/interpretation_tables/KIR%202025%20Genotyping%202V2%20Interpretation%20Table.pdf
Full-length group : *001, *011, *014-*017, *019-*022
Deletion group     : *003, *004, *006-*010, *012, *013, *018
If T1K's called alleles span both groups (rare/heterozygous), the sample is
reported as "FULL/DEL".

KIR3DP1
-------
The embedded true-typing table calls 3DP1=DEL for every sample it has (it
doesn't vary in that cohort), and no public allele->group mapping for 3DP1
was found the way one exists for 2DS4. So 3DP1 is only compared on
presence/absence; its DEL/FULL sub-typing is not attempted.
"""

import argparse
import csv
import glob
import os
import re
import sys

import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

# ---------------------------------------------------------------------------
# Ground truth: the 6 samples pasted in chat. Edit this table (or pass
# --truth-csv with columns: sample,2DL1,2DL2,2DL3,3DL1,3DL2,2DS1,2DS2,2DS3,
# 2DS4,2DS5,3DS1,2DP1,3DL3,2DL4,2DL5,3DP1) to add more samples.
# ---------------------------------------------------------------------------
GENES = ["2DL1", "2DL2", "2DL3", "3DL1", "3DL2", "2DS1", "2DS2", "2DS3", "2DS4",
         "2DS5", "3DS1", "2DP1", "3DL3", "2DL4", "2DL5", "3DP1"]

DEFAULT_TRUTH = {
    "001-063": ["+", "-", "+", "+", "+", "-", "-", "-", "DEL", "-", "-", "+", "+", "+", "-", "DEL"],
    "003-009": ["+", "+", "-", "+", "+", "+", "+", "-", "DEL", "+", "-", "+", "+", "+", "+", "DEL"],
    "001-037": ["+", "-", "+", "+", "+", "-", "-", "-", "FULL/DEL", "-", "-", "+", "+", "+", "-", "DEL"],
    "001-003": ["+", "-", "+", "+", "+", "+", "-", "-", "DEL", "+", "+", "+", "+", "+", "+", "DEL"],
    "001-018": ["+", "+", "+", "+", "+", "+", "+", "+", "DEL", "+", "+", "+", "+", "+", "+", "DEL"],
    "004-003": ["+", "-", "+", "+", "+", "-", "-", "-", "DEL", "-", "-", "+", "+", "+", "-", "DEL"],
}

# gene -> which raw T1K gene row(s) feed it (2DL5 combines the A and B loci)
RAW_ROW_MAP = {
    "2DL1": ["KIR2DL1"], "2DL2": ["KIR2DL2"], "2DL3": ["KIR2DL3"],
    "3DL1": ["KIR3DL1"], "3DL2": ["KIR3DL2"], "2DS1": ["KIR2DS1"],
    "2DS2": ["KIR2DS2"], "2DS3": ["KIR2DS3"], "2DS4": ["KIR2DS4"],
    "2DS5": ["KIR2DS5"], "3DS1": ["KIR3DS1"], "2DP1": ["KIR2DP1"],
    "3DL3": ["KIR3DL3"], "2DL4": ["KIR2DL4"], "2DL5": ["KIR2DL5A", "KIR2DL5B"],
    "3DP1": ["KIR3DP1"],
}

DEL_2DS4 = {"003", "004", "006", "007", "008", "009", "010", "012", "013", "018"}
FL_2DS4 = {"001", "011", "014", "015", "016", "017", "019", "020", "021", "022"}


def allele_number(allele_str):
    """'KIR2DS4*00601' or 'KIR2DS4*006' -> '006' (first 3 digits after *)."""
    m = re.search(r"\*(\d{3})", allele_str or "")
    return m.group(1) if m else None


def ds4_group(alleles):
    groups = set()
    for a in alleles:
        num = allele_number(a)
        if num in DEL_2DS4:
            groups.add("DEL")
        elif num in FL_2DS4:
            groups.add("FULL")
    if groups == {"DEL"}:
        return "DEL"
    if groups == {"FULL"}:
        return "FULL"
    if groups == {"DEL", "FULL"}:
        return "FULL/DEL"
    return "UNKNOWN"


def parse_t1k_file(path):
    """Return {gene_name: {'n':int, 'alleles':[...], 'abundances':[...], 'qualities':[...]}}."""
    rows = {}
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for line in reader:
            if not line or not line[0].strip():
                continue
            gene = line[0].strip()
            if not gene.upper().startswith("KIR"):
                continue  # skip stray header/blank lines if present
            n_alleles = int(float(line[1])) if len(line) > 1 and line[1] not in ("", ".") else 0
            alleles, abundances, qualities = [], [], []
            # columns are (allele, abundance, quality) triplets from index 2 onward
            for i in range(2, len(line), 3):
                triplet = line[i:i + 3]
                if len(triplet) < 3:
                    break
                a, ab, q = triplet
                if a and a != ".":
                    alleles.append(a)
                    try:
                        abundances.append(float(ab))
                    except ValueError:
                        abundances.append(0.0)
                    try:
                        qualities.append(float(q))
                    except ValueError:
                        qualities.append(-1.0)
            rows[gene] = {"n": n_alleles, "alleles": alleles, "abundances": abundances, "qualities": qualities}
    return rows


def call_gene(raw, gene_names, quality_threshold):
    """Combine 1-2 raw T1K rows (for 2DL5's A/B split) into one presence call."""
    all_alleles, all_quals = [], []
    total_n = 0
    for g in gene_names:
        r = raw.get(g, {"n": 0, "alleles": [], "abundances": [], "qualities": []})
        total_n += r["n"]
        all_alleles.extend(r["alleles"])
        all_quals.extend(r["qualities"])

    if total_n == 0 or not all_alleles:
        return "-", "", ""

    max_q = max(all_quals) if all_quals else -1
    alleles_str = " / ".join(all_alleles)
    qual_str = ", ".join(f"q={q:g}" for q in all_quals)

    if max_q >= quality_threshold:
        call = "+"
    else:
        call = "- (low-conf hit)"
    return call, alleles_str, qual_str


def compare_to_truth(gene, true_val, t1k_call, t1k_alleles):
    if gene == "2DS4":
        if t1k_call.startswith("-"):
            return "MATCH" if true_val == "-" else f"MISMATCH (T1K=absent, true={true_val})"
        group = ds4_group([a for a in t1k_alleles.split(" / ") if a])
        return "MATCH" if group == true_val else f"MISMATCH (T1K={group})"
    if gene == "3DP1":
        true_present = true_val in ("+", "DEL", "FULL", "FULL/DEL")
        t1k_present = t1k_call.startswith("+")
        if true_present and t1k_present:
            return "MATCH (presence only; DEL/FULL sub-type not checked)"
        return "MATCH" if true_present == t1k_present else "MISMATCH"
    # plain +/- genes
    true_symbol = true_val[0]
    t1k_symbol = t1k_call[0]
    return "MATCH" if true_symbol == t1k_symbol else "MISMATCH"


def extract_sample_id(filename, pattern):
    m = re.search(pattern, os.path.basename(filename), re.IGNORECASE)
    if not m:
        return None
    return "-".join(m.groups())


def load_truth(truth_csv):
    if truth_csv:
        df = pd.read_csv(truth_csv, dtype=str).set_index("sample")
        return {idx: [row[g] for g in GENES] for idx, row in df.iterrows()}
    return DEFAULT_TRUTH


# ---------------------------------------------------------------------------
# Excel report
# ---------------------------------------------------------------------------
FONT = "Arial"
HDR_FILL = PatternFill("solid", fgColor="1F4E78")
HDR_FONT = Font(name=FONT, bold=True, color="FFFFFF", size=10)
FILL_MATCH = PatternFill("solid", fgColor="C6EFCE")
FILL_MISMATCH = PatternFill("solid", fgColor="FFC7CE")
FILL_CAVEAT = PatternFill("solid", fgColor="FFEB9C")
THIN = Side(style="thin", color="B7B7B7")
BORDER = Border(left=THIN, right=THIN, top=THIN, bottom=THIN)


def style_header(ws, row, ncols):
    for c in range(1, ncols + 1):
        cell = ws.cell(row=row, column=c)
        cell.fill = HDR_FILL
        cell.font = HDR_FONT
        cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        cell.border = BORDER


def autosize(ws, widths):
    for i, w in enumerate(widths, start=1):
        ws.column_dimensions[get_column_letter(i)].width = w


def build_report(records, out_path, quality_threshold, no_truth_samples):
    wb = Workbook()

    # --- long format ---------------------------------------------------
    ws = wb.active
    ws.title = "Comparison_long"
    headers = ["Sample", "Gene", "True", "T1K_call", "T1K_confidence", "T1K_alleles", "Match", "Notes"]
    ws.append(headers)
    style_header(ws, 1, len(headers))
    for rec in records:
        ws.append([rec["sample"], "KIR" + rec["gene"], rec["true"], rec["t1k_call"],
                   rec["t1k_qual"], rec["t1k_alleles"], rec["match"], rec["notes"]])
    for row in ws.iter_rows(min_row=2, max_row=ws.max_row, min_col=1, max_col=len(headers)):
        for cell in row:
            cell.font = Font(name=FONT, size=10)
            cell.border = BORDER
            cell.alignment = Alignment(vertical="center", wrap_text=True)
        match_val = str(row[6].value)
        fill = FILL_MATCH if match_val.startswith("MATCH") and "presence only" not in match_val else \
            (FILL_CAVEAT if "presence only" in match_val else FILL_MISMATCH)
        row[6].fill = fill
    autosize(ws, [12, 10, 10, 20, 24, 30, 30, 55])
    ws.freeze_panes = "A2"

    # --- wide matrix (mirrors the true-typing table layout) ------------
    ws2 = wb.create_sheet("Comparison_wide")
    header_row = ["Sample"]
    for g in GENES:
        header_row += [f"{g}_True", f"{g}_T1K", f"{g}_Match"]
    ws2.append(header_row)
    style_header(ws2, 1, len(header_row))

    by_sample = {}
    for rec in records:
        by_sample.setdefault(rec["sample"], {})[rec["gene"]] = rec

    for sample in sorted(by_sample):
        row_vals = [sample]
        for g in GENES:
            rec = by_sample[sample].get(g)
            if rec is None:
                row_vals += ["", "", ""]
            else:
                row_vals += [rec["true"], rec["t1k_call"], rec["match"]]
        ws2.append(row_vals)

    for r in range(2, ws2.max_row + 1):
        for gi, g in enumerate(GENES):
            match_col = 2 + gi * 3 + 2  # 1-indexed column of the *_Match cell
            cell = ws2.cell(row=r, column=match_col)
            v = str(cell.value)
            if v.startswith("MATCH") and "presence only" not in v:
                fill = FILL_MATCH
            elif "presence only" in v:
                fill = FILL_CAVEAT
            elif v.startswith("MISMATCH"):
                fill = FILL_MISMATCH
            else:
                fill = None
            if fill:
                for c in (match_col - 2, match_col - 1, match_col):
                    ws2.cell(row=r, column=c).fill = fill
    for row in ws2.iter_rows(min_row=1, max_row=ws2.max_row, min_col=1, max_col=len(header_row)):
        for cell in row:
            cell.font = Font(name=FONT, size=9)
            cell.border = BORDER
            cell.alignment = Alignment(horizontal="center")
    autosize(ws2, [10] + [9] * (len(GENES) * 3))
    ws2.freeze_panes = "B2"

    # --- summary ---------------------------------------------------------
    ws3 = wb.create_sheet("Summary")
    ws3.append(["Sample", "N_genes_compared", "N_match", "N_mismatch", "N_caveat", "Concordance_%", "Mismatched_genes"])
    style_header(ws3, 1, 7)
    for sample in sorted(by_sample):
        genes_here = by_sample[sample]
        n_match = sum(1 for r in genes_here.values() if r["match"].startswith("MATCH") and "presence only" not in r["match"])
        n_caveat = sum(1 for r in genes_here.values() if "presence only" in r["match"])
        n_mismatch = sum(1 for r in genes_here.values() if r["match"].startswith("MISMATCH"))
        n_no_truth = sum(1 for r in genes_here.values() if r["match"] == "no ground truth")
        n_scored = n_match + n_caveat + n_mismatch
        mismatched = ", ".join("KIR" + g for g, r in genes_here.items() if r["match"].startswith("MISMATCH"))
        pct = round(100 * (n_match + n_caveat) / n_scored, 1) if n_scored else "n/a (no ground truth)"
        ws3.append([sample, n_scored, n_match, n_mismatch, n_caveat, pct, mismatched])
    if no_truth_samples:
        r = ws3.max_row + 2
        ws3.cell(row=r, column=1, value="Samples with T1K output but no ground-truth row (T1K calls still parsed, not scored):").font = Font(name=FONT, italic=True, size=9)
        ws3.cell(row=r + 1, column=1, value=", ".join(sorted(no_truth_samples))).font = Font(name=FONT, size=9)
    for row in ws3.iter_rows(min_row=1, max_row=ws3.max_row, min_col=1, max_col=7):
        for cell in row:
            cell.font = Font(name=FONT, size=10)
            cell.border = BORDER
    autosize(ws3, [12, 16, 10, 12, 10, 14, 50])
    ws3.freeze_panes = "A2"

    # --- notes -------------------------------------------------------------
    ws4 = wb.create_sheet("Notes_sources")
    lines = [
        ("Presence/absence threshold", f"quality >= {quality_threshold} across a gene's 1-2 alleles counts as present; "
                                        "below that it's flagged '- (low-conf hit)' and still treated as absent for scoring."),
        ("KIR2DS4 full-length vs deletion", "Per the GenDx KIR Genotyping SSO kit interpretation table (2025, 2V2): "
                                             "full-length group = *001, *011, *014-*017, *019-*022; deletion group = "
                                             "*003, *004, *006-*010, *012, *013, *018."),
        ("KIR3DP1", "Only compared on presence/absence -- no public allele->group (fl/del) mapping was found for 3DP1, "
                     "so any DEL/FULL label in your true-typing table isn't independently checked against T1K's alleles."),
        ("Regenerating this report", "Edit the DEFAULT_TRUTH dict in this script (or pass --truth-csv) to add samples; "
                                      "point --t1k-dir at the folder of *_KIR_genotype.tsv files."),
    ]
    ws4.append(["Topic", "Detail"])
    style_header(ws4, 1, 2)
    for t, d in lines:
        ws4.append([t, d])
    for row in ws4.iter_rows(min_row=1, max_row=ws4.max_row, min_col=1, max_col=2):
        for cell in row:
            cell.font = Font(name=FONT, size=10)
            cell.alignment = Alignment(vertical="top", wrap_text=True)
            cell.border = BORDER
    autosize(ws4, [28, 100])

    wb.save(out_path)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--t1k-dir", required=True, help="Folder containing *_KIR_genotype.tsv files (searched recursively for that glob).")
    ap.add_argument("--glob", default="*KIR_genotype*.tsv", help="Glob pattern for T1K files within --t1k-dir (default: *KIR_genotype*.tsv).")
    ap.add_argument("--out", required=True, help="Output .xlsx path.")
    ap.add_argument("--quality-threshold", type=float, default=1.0, help="Min T1K allele quality to call a gene present (default 1).")
    ap.add_argument("--filename-pattern", default=r"cos_(\d+)_(\d+)_KIR_genotype",
                     help=r"Regex with 2 capture groups pulled from the filename and joined with '-' to make the sample ID (default matches cos_004_003_KIR_genotype.tsv -> 004-003).")
    ap.add_argument("--truth-csv", default=None, help="Optional CSV (columns: sample,2DL1,2DL2,...,3DP1) to use instead of the 6 embedded samples.")
    args = ap.parse_args()

    truth = load_truth(args.truth_csv)

    files = sorted(glob.glob(os.path.join(args.t1k_dir, "**", args.glob), recursive=True))
    if not files:
        sys.exit(f"No files matching '{args.glob}' found under {args.t1k_dir}")

    records = []
    no_truth_samples = set()
    for f in files:
        sample = extract_sample_id(f, args.filename_pattern)
        if sample is None:
            print(f"WARNING: could not extract a sample ID from '{f}' using pattern "
                  f"'{args.filename_pattern}' -- skipping. Pass --filename-pattern to fix.", file=sys.stderr)
            continue
        raw = parse_t1k_file(f)
        true_row = truth.get(sample)
        if true_row is None:
            no_truth_samples.add(sample)
        for gi, gene in enumerate(GENES):
            gene_names = RAW_ROW_MAP[gene]
            t1k_call, t1k_alleles, t1k_qual = call_gene(raw, gene_names, args.quality_threshold)
            true_val = true_row[gi] if true_row else "?"
            match = compare_to_truth(gene, true_val, t1k_call, t1k_alleles) if true_row else "no ground truth"
            notes = ""
            if gene == "2DL5" and "low-conf hit" in t1k_call:
                notes = "Check abundance -- a 2DL5A hit >~5 with quality 0 is borderline, not clearly noise."
            records.append({
                "sample": sample, "gene": gene, "true": true_val,
                "t1k_call": t1k_call, "t1k_alleles": t1k_alleles, "t1k_qual": t1k_qual,
                "match": match, "notes": notes,
            })

    build_report(records, args.out, args.quality_threshold, no_truth_samples)
    print(f"Wrote {args.out} ({len(files)} T1K file(s) processed, "
          f"{len(files) - len(no_truth_samples)} matched to ground truth).")


if __name__ == "__main__":
    main()