#!/usr/bin/env python3
"""
hla_concordance.py

Compare HLA class I / non-classical typing calls from three tools
(OptiType, HLA-HD, T1K) against a ground-truth typing table (SSO/array/
whatever reference Amatta is using), and produce a single color-coded
Excel concordance report. Same overall design as kir_concordance.py:
one script, run it over a folder of each tool's output, point it at the
ground-truth table, get a workbook back.

------------------------------------------------------------------------
INPUT FORMATS EXPECTED (one file per sample per tool)
------------------------------------------------------------------------

OptiType (--optitype-dir), e.g. COS_001_037_result.tsv
    A tab-separated file that looks like a pandas DataFrame dump with an
    index column, one data row:
        \tA1\tA2\tB1\tB2\tC1\tC2\tReads\tObjective
        0\tA*02:01\tA*02:01\tB*18:01\tB*35:01\tC*04:01\tC*05:01\t3881.0\t3811.14...
    Only genes A, B, C are typed by OptiType.

HLA-HD (--hlahd-dir), e.g. COS_001_037_hlahd_final.txt (name is whatever
    you use; pattern is configurable). One gene per line, tab-separated,
    no header, gene code first then 2 (occasionally more) alleles:
        A\tHLA-A*02:01:01\t-
        B\tHLA-B*35:01:01\tHLA-B*18:01:01
        ...
        E\tHLA-E*01:01:01\tHLA-E*01:03:02\tHLA-E*01:01:02\tHLA-E*01:03:01
    "-" and "Not typed" both mean no call. Where more than 2 alleles are
    present (as HLA-HD sometimes emits for E/G/DRB-region genes) only the
    first two are used as allele 1 / allele 2 -- see Notes_sources sheet.

T1K (--t1k-dir), e.g. COS_001_037_HLA_genotype.tsv -- same shape as the
    KIR T1K genotype files: gene, n_alleles, then (allele, abundance,
    quality) triplets for up to 2 alleles, with an optional 3rd
    semicolon-packed "extra candidate" field tacked on the end that is
    ignored for the primary calls:
        HLA-A\t2\tHLA-A*24:02:01\t320.948316\t60\tHLA-A*31:01:02\t218.322963\t60\t
        HLA-L\t2\tHLA-L*01:01:02\t201.837294\t60\tHLA-L*01:03\t201.038564\t60\tHLA-L*01:01:01;68.921656;30

Ground truth (--truth-xlsx or --truth-csv): a table shaped exactly like
    the one Amatta uses for the true typing, with SAMPLE NAME as the
    first column and two header rows above the data:
        row 1 (gene group, merged cells):  <blank>  HLA-F  <blank>  HLA-G  <blank>  HLA-A  ...  HLA-C  ...  HLA-B  ...
        row 2 (sub-columns):   Sample Name  Allele 1  Allele 2 | Allele 1  Allele 2 | Allele 1  Group  Bw4  Allele 2  Group  Bw4  Genotype  Bw4 Genotype | Allele 1  Group  Allele 2  Group  Genotype | Allele 1  Group  Allele 2  Group  Genotype
    Run with --truth-template to generate a blank xlsx with exactly this
    header shape to fill in / to sanity-check your existing sheet against.

------------------------------------------------------------------------
USAGE
------------------------------------------------------------------------
    python3 hla_concordance.py \
        --optitype-dir /path/to/optitype_results \
        --hlahd-dir /path/to/hlahd_results \
        --t1k-dir /path/to/t1k_results \
        --truth-xlsx /path/to/ground_truth.xlsx \
        --out HLA_typing_comparison.xlsx

Any of --optitype-dir / --hlahd-dir / --t1k-dir can be omitted if you
don't have that tool's output yet -- the report just shows N/A for that
method/gene combination. Ground truth is required (or use --truth-template
to generate a blank one first).
"""

import argparse
import glob
import os
import re
import sys

import pandas as pd
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

# ---------------------------------------------------------------------------
# Constants / gene schema
# ---------------------------------------------------------------------------

# Order + shape of the ground-truth table's column blocks, exactly as pasted.
GENE_ORDER = ["F", "G", "A", "C", "B"]
GENE_SCHEMAS = {
    "F": ["Allele 1", "Allele 2"],
    "G": ["Allele 1", "Allele 2"],
    "A": ["Allele 1", "Group (A1)", "Bw4 (A1)", "Allele 2", "Group (A2)", "Bw4 (A2)", "Genotype", "Bw4 Genotype"],
    "C": ["Allele 1", "Group (A1)", "Allele 2", "Group (A2)", "Genotype"],
    "B": ["Allele 1", "Group (A1)", "Allele 2", "Group (A2)", "Genotype"],
}

# Genes actually scored for concordance (must have truth data for these).
CORE_GENES = ["A", "B", "C", "F", "G"]

MISSING_TOKENS = {"", "-", ".", "not typed", "nottyped", "nan", "none", "na", "n/a"}

FILL_MATCH = PatternFill(start_color="C6EFCE", end_color="C6EFCE", fill_type="solid")
FILL_PARTIAL = PatternFill(start_color="FFEB9C", end_color="FFEB9C", fill_type="solid")
FILL_MISMATCH = PatternFill(start_color="FFC7CE", end_color="FFC7CE", fill_type="solid")
FILL_NA = PatternFill(start_color="D9D9D9", end_color="D9D9D9", fill_type="solid")
FILL_HEADER = PatternFill(start_color="4472C4", end_color="4472C4", fill_type="solid")
FONT_HEADER = Font(color="FFFFFF", bold=True)


# ---------------------------------------------------------------------------
# Allele normalisation / matching
# ---------------------------------------------------------------------------

def normalize_allele(raw, resolution=2, gene=None):
    """Normalise an HLA allele string to 'LOCUS*ff:ff[:ff...]' at the given
    field resolution. Returns None for missing/absent calls.

    Accepts input with or without an 'HLA-' prefix, and with or without the
    locus letter already attached (e.g. 'HLA-A*02:01:01', 'A*02:01',
    '02:01' -- the last needs `gene` supplied to reattach the locus).
    """
    if raw is None:
        return None
    s = str(raw).strip()
    if s.lower() in MISSING_TOKENS:
        return None
    s = s.replace("HLA-", "").replace("hla-", "")
    if "*" in s:
        locus, fields = s.split("*", 1)
    elif gene:
        locus, fields = gene, s
    else:
        return None
    locus = locus.strip()
    parts = [p for p in fields.split(":") if p != ""]
    if not parts:
        return None
    parts = parts[:resolution]
    return f"{locus}*{':'.join(parts)}"


def fill_homozygous(pair):
    """If exactly one side of a (a1, a2) pair is called and the other is
    missing, treat the locus as homozygous (duplicate the called allele)
    rather than as a partial/missing call. This is how OptiType/HLA-HD/T1K
    all represent a homozygous locus when they don't explicitly repeat the
    allele twice -- HLA-HD and T1K commonly report just one allele plus a
    '-' / '.' placeholder for the second slot in that case.
    """
    a1, a2 = pair
    if a1 and not a2:
        return (a1, a1)
    if a2 and not a1:
        return (a2, a2)
    return (a1, a2)


def classify_match(truth_pair, call_pair):
    """Compare (truth_a1, truth_a2) vs (call_a1, call_a2) as unordered sets.
    Returns one of 'Match', 'Partial', 'Mismatch', 'N/A'.
    A single-allele call (other slot missing) is treated as an implicit
    homozygous call -- see fill_homozygous().
    """
    truth_pair = fill_homozygous(truth_pair)
    call_pair = fill_homozygous(call_pair)
    t = sorted(a for a in truth_pair if a)
    c = sorted(a for a in call_pair if a)
    if not t or not c:
        return "N/A"
    # Multiset-aware overlap (handles homozygous truth/call correctly)
    t_remaining = list(t)
    hits = 0
    for allele in c:
        if allele in t_remaining:
            t_remaining.remove(allele)
            hits += 1
    if hits == 0:
        return "Mismatch"
    if hits >= max(len(t), len(c)):
        return "Match"
    return "Partial"


# ---------------------------------------------------------------------------
# Sample ID extraction
# ---------------------------------------------------------------------------

def extract_sample_id(filename, pattern):
    base = os.path.basename(filename)
    m = re.search(pattern, base, flags=re.IGNORECASE)
    if m:
        groups = [g for g in m.groups() if g]
        if groups:
            return "-".join(groups)
        return m.group(0)
    stem = re.sub(r"\.(tsv|txt|csv)$", "", base, flags=re.IGNORECASE)
    print(f"  [warn] could not match sample-ID pattern against '{base}', using stem '{stem}'", file=sys.stderr)
    return stem


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

def parse_optitype_file(path, resolution):
    df = pd.read_csv(path, sep="\t", index_col=0)
    if df.empty:
        return {}
    if "Objective" in df.columns:
        df = df.sort_values("Objective", ascending=False)
    row = df.iloc[0]
    calls = {}
    for gene in ["A", "B", "C"]:
        a1 = normalize_allele(row.get(f"{gene}1"), resolution, gene)
        a2 = normalize_allele(row.get(f"{gene}2"), resolution, gene)
        calls[gene] = (a1, a2)
    return calls


def parse_hlahd_file(path, resolution):
    calls = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            fields = line.split("\t")
            gene = fields[0].strip()
            normed = [normalize_allele(a, resolution, gene) for a in fields[1:]]
            normed = [a for a in normed if a]
            a1 = normed[0] if len(normed) > 0 else None
            a2 = normed[1] if len(normed) > 1 else None
            calls[gene] = (a1, a2)
            if len(normed) > 2:
                calls.setdefault("_extra_alleles", {})[gene] = normed[2:]
    return calls


def parse_t1k_file(path, resolution, quality_threshold):
    calls = {}
    low_conf = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            fields = line.split("\t")
            gene_raw = fields[0].strip()
            gene = gene_raw.replace("HLA-", "").replace("hla-", "")
            if len(fields) < 8:
                continue

            def resolve(raw, q_raw):
                allele = normalize_allele(raw, resolution, gene)
                if allele is None:
                    return None, False
                try:
                    qv = float(q_raw)
                except (TypeError, ValueError):
                    qv = None
                lc = qv is not None and qv < quality_threshold
                return allele, lc

            a1, lc1 = resolve(fields[2], fields[4])
            a2, lc2 = resolve(fields[5], fields[7])
            calls[gene] = (a1, a2)
            low_conf[gene] = (lc1, lc2)
    return calls, low_conf


def load_calls_dir(directory, pattern, parser, *parser_extra_args):
    """Run `parser(path, *parser_extra_args)` over every file in `directory`
    matching the sample-ID `pattern`, keyed by extracted sample ID.
    """
    result = {}
    if not directory:
        return result
    files = sorted(glob.glob(os.path.join(directory, "*")))
    if not files:
        print(f"  [warn] no files found in {directory}", file=sys.stderr)
    for f in files:
        if os.path.isdir(f):
            continue
        sample_id = extract_sample_id(f, pattern)
        try:
            result[sample_id] = parser(f, *parser_extra_args)
        except Exception as e:
            print(f"  [warn] failed to parse {f}: {e}", file=sys.stderr)
    return result


# ---------------------------------------------------------------------------
# Ground truth
# ---------------------------------------------------------------------------

def load_truth(path, resolution):
    if str(path).lower().endswith((".xlsx", ".xls")):
        raw = pd.read_excel(path, header=None, dtype=str)
    else:
        raw = pd.read_csv(path, header=None, dtype=str, sep=None, engine="python")

    gene_row_raw = raw.iloc[0].tolist()
    sub_row_raw = raw.iloc[1].tolist()

    filled_genes = []
    last = ""
    for g in gene_row_raw:
        g = "" if pd.isna(g) else str(g).strip()
        if g:
            last = g.replace("HLA-", "")
        filled_genes.append(last if g or last else "")
    # Reset "last" propagation before the first real gene label appears
    # (keeps the Sample Name column un-labelled instead of inheriting
    # nothing spurious -- filled_genes[0] should stay "").
    filled_genes[0] = ""

    data = raw.iloc[2:].reset_index(drop=True)

    truth = {}
    sample_col_idx = 0
    for _, row in data.iterrows():
        sid_raw = row[sample_col_idx]
        if pd.isna(sid_raw) or not str(sid_raw).strip():
            continue
        sid = str(sid_raw).strip()
        gene_data = {}
        for gene in GENE_ORDER:
            idxs = [i for i, g in enumerate(filled_genes) if g == gene]
            schema = GENE_SCHEMAS[gene]
            if not idxs:
                gene_data[gene] = {k: None for k in schema}
                continue
            if len(idxs) != len(schema):
                print(
                    f"  [warn] sample {sid}: expected {len(schema)} columns for HLA-{gene}, "
                    f"found {len(idxs)} -- truncating/padding positionally",
                    file=sys.stderr,
                )
            entry = {}
            for pos, key in enumerate(schema):
                if pos < len(idxs):
                    entry[key] = row[idxs[pos]]
                else:
                    entry[key] = None
            gene_data[gene] = entry
        truth[sid] = gene_data

    # Normalise allele fields now that we have gene context
    for sid, genes in truth.items():
        for gene, entry in genes.items():
            entry["Allele 1"] = normalize_allele(entry.get("Allele 1"), resolution, gene)
            entry["Allele 2"] = normalize_allele(entry.get("Allele 2"), resolution, gene)

    return truth


def write_truth_template(path):
    rows = []
    row1 = [""]
    row2 = ["Sample Name"]
    for gene in GENE_ORDER:
        schema = GENE_SCHEMAS[gene]
        row1.append(f"HLA-{gene}")
        row1.extend([""] * (len(schema) - 1))
        row2.extend(schema)
    rows.append(row1)
    rows.append(row2)
    df = pd.DataFrame(rows)
    df.to_excel(path, header=False, index=False)
    print(f"Wrote blank ground-truth template to {path}")
    print("Fill in one row per sample under the two header rows (Sample Name in column A), "
          "then pass this file back in with --truth-xlsx.")


# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------

def build_comparison(truth, optitype_calls, hlahd_calls, t1k_calls, t1k_lowconf, compare_genes):
    sample_ids = sorted(set(truth) | set(optitype_calls) | set(hlahd_calls) | set(t1k_calls))
    long_rows = []
    for sid in sample_ids:
        t_genes = truth.get(sid, {})
        o_genes = optitype_calls.get(sid, {})
        h_genes = hlahd_calls.get(sid, {})
        k_genes = t1k_calls.get(sid, {})
        k_lc = t1k_lowconf.get(sid, {})
        for gene in compare_genes:
            t_entry = t_genes.get(gene, {})
            t_pair = (t_entry.get("Allele 1"), t_entry.get("Allele 2"))

            o_pair = o_genes.get(gene, (None, None))
            h_pair = h_genes.get(gene, (None, None))
            k_pair = k_genes.get(gene, (None, None))
            k_lc_pair = k_lc.get(gene, (False, False))

            row = {
                "Sample": sid,
                "Gene": f"HLA-{gene}",
                "Truth_Allele1": t_pair[0],
                "Truth_Allele2": t_pair[1],
                "Truth_Group_A1": t_entry.get("Group (A1)"),
                "Truth_Group_A2": t_entry.get("Group (A2)"),
                "Truth_Bw4_A1": t_entry.get("Bw4 (A1)"),
                "Truth_Bw4_A2": t_entry.get("Bw4 (A2)"),
                "Truth_Genotype": t_entry.get("Genotype"),
                "Truth_Bw4_Genotype": t_entry.get("Bw4 Genotype"),
                "OptiType_Allele1": o_pair[0],
                "OptiType_Allele2": o_pair[1],
                "OptiType_Match": classify_match(t_pair, o_pair) if gene in ("A", "B", "C") else "N/A (not typed by tool)",
                "HLAHD_Allele1": h_pair[0],
                "HLAHD_Allele2": h_pair[1],
                "HLAHD_Match": classify_match(t_pair, h_pair),
                "T1K_Allele1": k_pair[0],
                "T1K_Allele2": k_pair[1],
                "T1K_Match": classify_match(t_pair, k_pair),
                "T1K_LowConfidence": any(k_lc_pair),
            }
            long_rows.append(row)
    return pd.DataFrame(long_rows)


def build_summary(long_df):
    rows = []
    for method in ["OptiType", "HLAHD", "T1K"]:
        col = f"{method}_Match"
        for sid, sub in long_df.groupby("Sample"):
            counts = sub[col].value_counts()
            assessable = counts.get("Match", 0) + counts.get("Partial", 0) + counts.get("Mismatch", 0)
            score = counts.get("Match", 0) + 0.5 * counts.get("Partial", 0)
            pct = round(100 * score / assessable, 1) if assessable else None
            mismatched_genes = sub.loc[sub[col].isin(["Partial", "Mismatch"]), "Gene"].tolist()
            rows.append({
                "Sample": sid,
                "Method": method,
                "Genes_assessed": assessable,
                "Full_matches": int(counts.get("Match", 0)),
                "Partial_matches": int(counts.get("Partial", 0)),
                "Mismatches": int(counts.get("Mismatch", 0)),
                "Concordance_%": pct,
                "Discordant_genes": ", ".join(mismatched_genes) if mismatched_genes else "",
            })
    return pd.DataFrame(rows).sort_values(["Sample", "Method"]).reset_index(drop=True)


def build_extended(optitype_calls, hlahd_calls, t1k_calls, compare_genes):
    """All loci typed by HLA-HD / T1K outside the core compared set --
    informational only, no ground truth available for these."""
    sample_ids = sorted(set(hlahd_calls) | set(t1k_calls))
    rows = []
    for sid in sample_ids:
        h_genes = {g: v for g, v in hlahd_calls.get(sid, {}).items() if g != "_extra_alleles"}
        k_genes = t1k_calls.get(sid, {})
        all_genes = sorted(set(h_genes) | set(k_genes))
        for gene in all_genes:
            if gene in compare_genes:
                continue
            h_pair = h_genes.get(gene, (None, None))
            k_pair = k_genes.get(gene, (None, None))
            if h_pair == (None, None) and k_pair == (None, None):
                continue
            rows.append({
                "Sample": sid,
                "Gene": gene,
                "HLAHD_Allele1": h_pair[0],
                "HLAHD_Allele2": h_pair[1],
                "T1K_Allele1": k_pair[0],
                "T1K_Allele2": k_pair[1],
            })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Excel output
# ---------------------------------------------------------------------------

def autosize(ws, df):
    for i, col in enumerate(df.columns, start=1):
        col_max = df[col].map(lambda x: len(str(x))).max() if len(df) else 0
        width = max(len(str(col)), col_max) + 2
        ws.column_dimensions[get_column_letter(i)].width = min(width, 40)


def style_header(ws, ncols):
    for c in range(1, ncols + 1):
        cell = ws.cell(row=1, column=c)
        cell.fill = FILL_HEADER
        cell.font = FONT_HEADER
        cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
    ws.freeze_panes = "A2"


def color_match_columns(ws, df, match_cols):
    fills = {"Match": FILL_MATCH, "Partial": FILL_PARTIAL, "Mismatch": FILL_MISMATCH}
    for col_name in match_cols:
        if col_name not in df.columns:
            continue
        col_idx = df.columns.get_loc(col_name) + 1
        for r, val in enumerate(df[col_name], start=2):
            fill = fills.get(val, FILL_NA)
            ws.cell(row=r, column=col_idx).fill = fill


def write_workbook(out_path, long_df, summary_df, extended_df, notes_lines):
    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        long_df.to_excel(writer, sheet_name="Comparison_long", index=False)
        summary_df.to_excel(writer, sheet_name="Summary", index=False)
        if not extended_df.empty:
            extended_df.to_excel(writer, sheet_name="Extended_other_loci", index=False)
        pd.DataFrame({"Notes": notes_lines}).to_excel(writer, sheet_name="Notes_sources", index=False, header=False)

        wb = writer.book

        ws = wb["Comparison_long"]
        style_header(ws, len(long_df.columns))
        color_match_columns(ws, long_df, ["OptiType_Match", "HLAHD_Match", "T1K_Match"])
        autosize(ws, long_df)

        ws = wb["Summary"]
        style_header(ws, len(summary_df.columns))
        autosize(ws, summary_df)
        pct_col = summary_df.columns.get_loc("Concordance_%") + 1
        for r, val in enumerate(summary_df["Concordance_%"], start=2):
            cell = ws.cell(row=r, column=pct_col)
            if val is None:
                cell.fill = FILL_NA
            elif val >= 99.9:
                cell.fill = FILL_MATCH
            elif val >= 50:
                cell.fill = FILL_PARTIAL
            else:
                cell.fill = FILL_MISMATCH

        if not extended_df.empty:
            ws = wb["Extended_other_loci"]
            style_header(ws, len(extended_df.columns))
            autosize(ws, extended_df)

        ws = wb["Notes_sources"]
        ws.column_dimensions["A"].width = 120
        for r in range(1, len(notes_lines) + 1):
            ws.cell(row=r, column=1).alignment = Alignment(wrap_text=True, vertical="top")

    print(f"Wrote {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

DEFAULT_NOTES = [
    "HLA typing concordance: OptiType vs HLA-HD vs T1K vs ground truth",
    "",
    "Methodology:",
    f"- All alleles truncated to {{resolution}}-field resolution before comparison (e.g. A*02:01:01 -> A*02:01), "
    "since OptiType only reports 2-field calls and this is the level the truth table appears to use.",
    "- Match/Partial/Mismatch/N/A is computed per gene per sample per method as an UNORDERED comparison of the "
    "allele pair against truth (a swapped Allele1/Allele2 order still counts as a Match). 'Partial' means exactly "
    "one of the two alleles agrees; 'Mismatch' means neither does; 'N/A' means truth or the tool call was missing.",
    "- If a locus has only ONE allele called (the other slot is blank/'-'/'.') on either the truth or the tool "
    "side, that's treated as an implicit homozygous call (the single allele is compared as if duplicated) rather "
    "than penalised as a partial call -- this is how OptiType/HLA-HD/T1K usually represent homozygosity when they "
    "don't explicitly repeat the allele twice. A genuinely different single call still shows up as Partial or "
    "Mismatch against a heterozygous truth, since only one of the two truth alleles can match.",
    "- Genes scored for concordance: HLA-A, -B, -C, -F, -G (the loci present in the ground-truth table). OptiType "
    "only types A/B/C, so its F and G cells are marked 'N/A (not typed by tool)' rather than scored as mismatches.",
    "- Genes HLA-HD / T1K type beyond A/B/C/F/G (E, H, J, K, L, DRB1, DQA1/B1, DPA1/B1, DMA/B, DOA/B, MICA/B, "
    "TAP1/2, etc.) are listed in the 'Extended_other_loci' sheet for reference only -- there is no ground truth to "
    "score them against here.",
    "- HLA-HD occasionally reports more than 2 alleles for a locus (e.g. HLA-E in the example data had 4 entries). "
    "Only the first two are used as Allele 1 / Allele 2; anything beyond that is dropped (not currently surfaced "
    "anywhere -- flag if you want the extras captured too).",
    "- T1K quality/abundance columns are parsed and a call is flagged T1K_LowConfidence=True if either allele's "
    "reported quality is below --quality-threshold (default 1). Unlike the KIR script, this does NOT force the "
    "call to 'absent' for these loci, since HLA-A/B/C/F/G are always diploid (not presence/absence like KIR) -- "
    "it's only a confidence flag for you to sanity-check, matching how DRB3 quality=0 showed up in the example data.",
    "- T1K's trailing semicolon-packed 'extra candidate' field (3rd+ allele guess, e.g. "
    "'HLA-L*01:01:01;68.921656;30') is parsed out but not used for the primary Allele1/Allele2 call, same as the "
    "KIR script.",
    "- Ground-truth Group and Bw4 columns are carried through from the truth table as-is for reference "
    "(Truth_Group_A1/A2, Truth_Bw4_A1/A2, Truth_Genotype, Truth_Bw4_Genotype columns in Comparison_long) but are "
    "NOT independently recomputed from the tool calls -- that needs an allele-to-G-group / Bw4-epitope lookup "
    "table which wasn't provided. Add one later if you want OptiType/HLA-HD/T1K's own Bw4 status verified rather "
    "than just copied from truth.",
    "",
    "Sheets:",
    "- Comparison_long: one row per sample x gene, truth + all 3 tool calls + per-tool match verdict, color-coded.",
    "- Summary: per sample x method, concordance % (Match=1, Partial=0.5, Mismatch=0, out of genes assessable) "
    "and list of discordant genes.",
    "- Extended_other_loci: everything HLA-HD/T1K typed outside A/B/C/F/G, informational only.",
]


def main():
    p = argparse.ArgumentParser(description="Compare OptiType/HLA-HD/T1K HLA typing against ground truth.")
    p.add_argument("--optitype-dir", help="Folder of OptiType *_result.tsv files")
    p.add_argument("--optitype-pattern", default=r"(\d+)_(\d+)_result", help="Regex to pull sample ID from OptiType filenames")
    p.add_argument("--hlahd-dir", help="Folder of HLA-HD final-result text files")
    p.add_argument("--hlahd-pattern", default=r"(\d+)_(\d+)", help="Regex to pull sample ID from HLA-HD filenames")
    p.add_argument("--t1k-dir", help="Folder of T1K *_HLA_genotype.tsv files")
    p.add_argument("--t1k-pattern", default=r"(\d+)_(\d+)_HLA_genotype", help="Regex to pull sample ID from T1K filenames")
    p.add_argument("--truth-xlsx", help="Ground truth workbook (2 header rows, see module docstring)")
    p.add_argument("--truth-csv", help="Ground truth CSV/TSV, same 2-header-row shape as --truth-xlsx")
    p.add_argument("--truth-template", metavar="PATH", help="Write a blank ground-truth template to PATH and exit")
    p.add_argument("--resolution", type=int, default=2, help="Allele field resolution to compare at (default 2, e.g. A*02:01)")
    p.add_argument("--quality-threshold", type=float, default=1.0, help="T1K quality below this -> flagged T1K_LowConfidence (default 1)")
    p.add_argument("--out", default="HLA_typing_comparison.xlsx", help="Output workbook path")
    args = p.parse_args()

    if args.truth_template:
        write_truth_template(args.truth_template)
        return

    if not args.truth_xlsx and not args.truth_csv:
        p.error("Need --truth-xlsx or --truth-csv (or run with --truth-template to generate one first)")

    truth_path = args.truth_xlsx or args.truth_csv
    print(f"Loading ground truth from {truth_path} ...")
    truth = load_truth(truth_path, args.resolution)
    print(f"  {len(truth)} sample(s) in ground truth")

    print("Loading OptiType calls ...")
    optitype_calls = load_calls_dir(args.optitype_dir, args.optitype_pattern, parse_optitype_file, args.resolution)
    print(f"  {len(optitype_calls)} sample(s)")

    print("Loading HLA-HD calls ...")
    hlahd_calls = load_calls_dir(args.hlahd_dir, args.hlahd_pattern, parse_hlahd_file, args.resolution)
    print(f"  {len(hlahd_calls)} sample(s)")

    print("Loading T1K calls ...")
    t1k_calls = {}
    t1k_lowconf = {}
    if args.t1k_dir:
        for f in sorted(glob.glob(os.path.join(args.t1k_dir, "*"))):
            if os.path.isdir(f):
                continue
            sid = extract_sample_id(f, args.t1k_pattern)
            try:
                calls, lc = parse_t1k_file(f, args.resolution, args.quality_threshold)
                t1k_calls[sid] = calls
                t1k_lowconf[sid] = lc
            except Exception as e:
                print(f"  [warn] failed to parse {f}: {e}", file=sys.stderr)
    print(f"  {len(t1k_calls)} sample(s)")

    long_df = build_comparison(truth, optitype_calls, hlahd_calls, t1k_calls, t1k_lowconf, CORE_GENES)
    summary_df = build_summary(long_df)
    extended_df = build_extended(optitype_calls, hlahd_calls, t1k_calls, CORE_GENES)

    notes = [DEFAULT_NOTES[0], DEFAULT_NOTES[1], DEFAULT_NOTES[2]]
    notes += [line.format(resolution=args.resolution) if "{resolution}" in line else line for line in DEFAULT_NOTES[3:]]

    write_workbook(args.out, long_df, summary_df, extended_df, notes)

    print("\nOverall concordance by method:")
    overall = summary_df.groupby("Method").apply(
        lambda d: round(100 * (d["Full_matches"].sum() + 0.5 * d["Partial_matches"].sum()) / d["Genes_assessed"].sum(), 1)
        if d["Genes_assessed"].sum() else None
    )
    for method, pct in overall.items():
        print(f"  {method}: {pct}%")


if __name__ == "__main__":
    main()