#!/usr/bin/env python3
"""Per-genome KEGG module completeness from the RAW DRAM annotations.

Why not the distillate workbook. Its curated module sheets do not carry
K01034/K01035, the butyryl-CoA:acetate CoA-transferase route, which is the
commoner of the two butyrate routes in the gut. A module layer built from that
workbook is structurally blind to it. The raw per-gene table does carry it
(25 and 24 genomes here), so module calls are made from `genome_ko_counts.tsv.gz`,
which `build_full_ko_table.py` derived from the 1.37 M-row annotation file.

Why not fraction-of-KOs-present. A module that offers three alternative enzymes
for one step contributes three KOs but needs only one of them. Scoring the
fraction of a module's KOs that a genome carries therefore penalises genomes
that simply picked one alternative, and the error is largest in exactly the
modules with the most alternatives. This evaluates the real KEGG DEFINITION
expression instead.

DEFINITION grammar, as KEGG writes it:
    space   sequential steps, all required       "K00844 K01810"
    comma   alternatives, any one suffices       "K00844,K12407"
    plus    subunits of one complex, all needed  "K00844+K12407"
    minus   optional component, ignored          "K00844-K12407"
    ()      grouping
    --      an unresolved step, counted as unsatisfiable

Completeness is the fraction of TOP-LEVEL steps satisfied, which is how KEGG
itself scores a module, not the fraction of KOs matched.

Outputs -> <outdir>/genome_module_completeness.tsv.gz  (long, one row per
           genome x module with completeness > 0)
           module_definitions_parsed.tsv               (audit of the parse)

Usage: build_kegg_module_completeness.py <outdir>
"""
import csv
import gzip
import os
import re
import sys
from collections import defaultdict

ROOT = os.path.dirname(os.path.abspath(__file__))
KEGG = os.path.join(ROOT, "result2/n76/fig4/gene_catalogue/kegg_reference_20260919")
KO_TABLE = os.path.join(ROOT, "result2/n76/metabolic_support/annotations_full/genome_ko_counts.tsv.gz")
DEFS = os.path.join(KEGG, "kegg_module_definitions.tsv")
HIER = os.path.join(KEGG, "ko00002_module_hierarchy.json")

KO_RE = re.compile(r"K\d{5}")


def split_top(expr, sep):
    """Split on `sep` only at bracket depth zero."""
    out, depth, cur = [], 0, []
    for ch in expr:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        if ch == sep and depth == 0:
            out.append("".join(cur))
            cur = []
        else:
            cur.append(ch)
    out.append("".join(cur))
    return [p for p in out if p.strip()]


def evaluate(expr, kos):
    """True if the genome's KO set satisfies this definition fragment."""
    expr = expr.strip()
    if not expr:
        return True
    # strip one layer of fully-enclosing brackets
    while expr.startswith("(") and expr.endswith(")"):
        depth = 0
        enclosing = True
        for i, ch in enumerate(expr):
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
                if depth == 0 and i < len(expr) - 1:
                    enclosing = False
                    break
        if not enclosing:
            break
        expr = expr[1:-1].strip()

    parts = split_top(expr, " ")
    if len(parts) > 1:                               # AND over sequential steps
        return all(evaluate(p, kos) for p in parts)
    parts = split_top(expr, ",")
    if len(parts) > 1:                               # OR over alternatives
        return any(evaluate(p, kos) for p in parts)
    parts = split_top(expr, "+")
    if len(parts) > 1:                               # AND over complex subunits
        return all(evaluate(p, kos) for p in parts)
    parts = split_top(expr, "-")
    if len(parts) > 1:                               # trailing parts are optional
        return evaluate(parts[0], kos)

    if expr.startswith("--") or expr == "--":
        return False                                 # unresolved step
    if expr.startswith("M"):
        return False                                 # nested module, not expanded
    return expr in kos


def steps(definition):
    """Top-level steps of a module definition."""
    return split_top(definition.strip(), " ")


def main():
    outdir = sys.argv[1]
    os.makedirs(outdir, exist_ok=True)

    defs = {}
    with open(DEFS) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            defs[r["module_id"]] = (r["module_name"], r["definition"])
    print(f"module definitions: {len(defs)}")

    # classification from the BRITE hierarchy, for grouping in the figures
    klass = {}
    try:
        import json
        def walk(node, path=()):
            ch = node.get("children")
            nm = node.get("name", "")
            if not ch:
                mid = nm.split()[0] if nm else ""
                if mid.startswith("M"):
                    klass[mid] = " | ".join(path[1:4])
                return
            for c in ch:
                walk(c, path + (nm,))
        walk(json.load(open(HIER)))
    except Exception as exc:
        print("hierarchy unavailable:", exc)

    per_genome = defaultdict(set)
    with gzip.open(KO_TABLE, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            per_genome[r["catalog_id"]].add(r["ko"])
    print(f"genomes: {len(per_genome)}")

    with open(os.path.join(outdir, "module_definitions_parsed.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["module_id", "module_name", "n_steps", "n_kos", "class", "definition"])
        for mid, (name, d) in sorted(defs.items()):
            w.writerow([mid, name, len(steps(d)), len(set(KO_RE.findall(d))),
                        klass.get(mid, ""), d])

    dest = os.path.join(outdir, "genome_module_completeness.tsv.gz")
    n = 0
    with gzip.open(dest, "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["catalog_id", "module_id", "module_name", "class",
                    "steps_total", "steps_present", "completeness", "n_ko_hit"])
        for gid, kos in sorted(per_genome.items()):
            for mid, (name, d) in defs.items():
                st = steps(d)
                if not st:
                    continue
                hit = sum(1 for s in st if evaluate(s, kos))
                if hit == 0:
                    continue
                nko = len(set(KO_RE.findall(d)) & kos)
                w.writerow([gid, mid, name, klass.get(mid, ""), len(st), hit,
                            f"{hit/len(st):.4f}", nko])
                n += 1
    print(f"WROTE {dest}  ({n} genome x module rows)")


if __name__ == "__main__":
    main()
