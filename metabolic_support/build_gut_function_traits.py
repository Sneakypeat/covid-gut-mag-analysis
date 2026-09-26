#!/usr/bin/env python3
"""Nonexclusive annotation evidence; no ecological/clinical guild inference."""
import argparse
import collections
import csv
import datetime
import hashlib
import json
import math
from pathlib import Path
import re
import xml.etree.ElementTree as ET
import zipfile

NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
PLANT_TERMS = ("Amorphous Cellulose", "Crystalline Cellulose", "Arabinan",
               "Beta-galactan (pectic galactan)", "Beta-mannan", "Starch",
               "Pectin", "Xylan", "Xyloglucan", "Mixed-Linkage glucans")
ACIDS = {"lac__D", "lac__L", "succ", "ac"}
CORE_MARKERS = {"K00634", "K00929", "K00625", "K00925", "GH33", "GH29",
                "GH95", "GH20", "GH35"}
TRAITS = (
    ("plant_backbone_cazyme", "Plant-glycan CAZyme repertoire",
     "any positive marker with frozen plant Backbone Cleavage annotation",
     "Family-level association; not validated extracellular substrate degradation"),
    ("butyrate_terminal_pair", "Butyrate kinase-route terminal pair",
     "K00634 AND K00929", "Not a complete pathway or census of butyrate producers"),
    ("acetate_terminal_pair", "Acetate terminal pair", "K00625 AND K00925",
     "Not net acetate production, flux or complete upstream pathway"),
    ("strict_acid_uptake", "Strict organic-acid uptake potential",
     "strict requirements intersect {lac__D,lac__L,succ,ac}",
     "Source-seed plus inward gene-grounded transport; not consumption or cross-feeding"),
    ("mucin_gh_repertoire", "Mucin-associated GH repertoire",
     "GH33 AND (GH29 OR GH95) AND (GH20 OR GH35)",
     "Putative non-specific repertoire; not necessary/sufficient for mucin degradation"),
)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_table(path, delimiter="\t"):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def keyed(rows, key):
    out = {row[key]: row for row in rows}
    if len(out) != len(rows):
        raise ValueError(f"Duplicate {key}")
    return out


def write_table(path, rows, fields=None):
    if not rows and not fields:
        raise ValueError(f"No schema for {path}")
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def column_name(cell):
    return re.sub(r"\d+$", "", cell.attrib["r"])


def cell_value(cell, strings):
    if cell.find("m:f", NS) is not None:
        raise ValueError("Formula cell in frozen count workbook")
    value = cell.find("m:v", NS)
    if cell.attrib.get("t") == "s":
        return strings[int(value.text)] if value is not None else ""
    inline = cell.find("m:is", NS)
    if inline is not None:
        return "".join(inline.itertext())
    return value.text if value is not None else ""


def numeric_count(value):
    if value == "":
        return None
    count = float(value)
    if not math.isfinite(count) or count < 0 or count != int(count):
        raise ValueError(f"Invalid gene count: {value}")
    return int(count)


def read_workbook(path, genomes):
    """Read only the MISC/carbon-utilization count sheets; exact genome join."""
    counts, definitions, plant = {}, [], set()
    with zipfile.ZipFile(path) as archive:
        strings = []
        if "xl/sharedStrings.xml" in archive.namelist():
            strings = ["".join(x.itertext()) for x in ET.fromstring(archive.read("xl/sharedStrings.xml"))]
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        sheets = list(workbook.find("m:sheets", NS))
        if [s.attrib["name"] for s in sheets[:2]] != ["MISC", "carbon utilization"]:
            raise ValueError("Unexpected DRAM workbook sheet mapping")
        for number, sheet in [(1, "MISC"), (2, "carbon utilization")]:
            rows = ET.fromstring(archive.read(f"xl/worksheets/sheet{number}.xml")).find("m:sheetData", NS)
            header = {column_name(c): cell_value(c, strings) for c in rows[0]}
            if len(set(header.values())) != len(header):
                raise ValueError("Duplicate workbook header")
            if not genomes <= set(header.values()):
                raise ValueError("Missing catalogue genomes in workbook")
            known_metadata = {"gene_id", "gene_description", "module", "header", "subheader",
                              "EC", "Notes", "oxygen", "specific_reaction"}
            if set(header.values()) - genomes - known_metadata:
                raise ValueError("Unexpected genome or metadata column")
            gene_columns = {col: name for col, name in header.items() if name in genomes}
            for xmlrow in list(rows)[1:]:
                values = {header[column_name(c)]: cell_value(c, strings) for c in xmlrow}
                gene = values.get("gene_id", "")
                if not gene:
                    continue
                if re.fullmatch(r"(?:GH|AA|PL|CE|CBM|GT)\d+|K\d{5}", gene) is None:
                    continue
                vector = {name: numeric_count(values.get(name, "")) for name in genomes}
                if gene in counts and counts[gene] != vector:
                    raise ValueError(f"Conflicting repeated gene counts: {gene}")
                counts[gene] = vector
                subheader = values.get("subheader", "")
                is_plant = any(f"{term} Backbone Cleavage" in subheader for term in PLANT_TERMS)
                if is_plant:
                    plant.add(gene)
                definitions.append(dict(gene_id=gene, sheet=sheet,
                    worksheet_row=xmlrow.attrib["r"], gene_description=values.get("gene_description", ""),
                    module=values.get("module", ""), substrate_description=subheader,
                    plant_backbone_selected=is_plant))
    if not CORE_MARKERS <= counts.keys() or not plant:
        raise ValueError("Required marker definitions are absent from the workbook")
    return counts, definitions, plant


def has(count):
    return None if count is None else count > 0


def any_evidence(values):
    values = list(values)
    return True if True in values else (None if None in values else False)


def all_evidence(values):
    values = list(values)
    return False if False in values else (None if None in values else True)


def trait_values(markers, plant, requirements):
    get = lambda gene: has(markers.get(gene))
    return {
        "plant_backbone_cazyme": any_evidence(get(gene) for gene in plant),
        "butyrate_terminal_pair": all_evidence([get("K00634"), get("K00929")]),
        "acetate_terminal_pair": all_evidence([get("K00625"), get("K00925")]),
        "strict_acid_uptake": None if requirements is None else bool(ACIDS & set(requirements)),
        "mucin_gh_repertoire": all_evidence([get("GH33"), any_evidence([get("GH29"), get("GH95")]),
                                            any_evidence([get("GH20"), get("GH35")])]),
    }


def state(value):
    return "unknown" if value is None else ("positive" if value else "evidence_not_detected")


def direction(row):
    if row is None or row.get("diff_GroupCase", "").upper() not in {"TRUE", "FALSE"}:
        return "Not tested"
    if row["diff_GroupCase"].upper() == "FALSE":
        return "Not significant"
    lfc = float(row["lfc_GroupCase"])
    if not math.isfinite(lfc) or lfc == 0:
        raise ValueError("Invalid significant ANCOM effect")
    return "Enriched" if lfc > 0 else "Depleted"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    base = root / "result2/n76/metabolic_support"
    out = args.out or base / "gut_traits_v1"
    paths = dict(workbook=root / "result2/inputfiles/new/DRAM_metabolism_summary.xlsx",
                 product=root / "result2/inputfiles/new/DRAM_metabolism_product.tsv",
                 taxonomy=root / "result2/inputfiles/new/MAG_quality_taxonomy.tsv",
                 primary=base / "MAG_manifest_primary_bacteria294.tsv",
                 sets=base / "graph/model_support_sets.json", graph_gate=base / "graph/extraction_audit.json",
                 ancom=root / "result2/n76/ANCOMBC2_MAG_results_n76.csv",
                 methods=base / "gut_traits_v1/TRAIT_METHODS_v1.md", script=Path(__file__))
    if not paths["methods"].is_file():
        raise ValueError("Trait method lock must exist before aggregation")
    if out.exists() and any(p.name != "TRAIT_METHODS_v1.md" for p in out.iterdir()):
        raise ValueError("Refusing to overwrite existing trait outputs; choose a fresh --out")
    hashes = {key: sha(path) for key, path in paths.items()}
    genomes = set(keyed(read_table(paths["product"]), "genome"))
    primary = set(keyed(read_table(paths["primary"]), "catalog_id"))
    taxonomy = keyed(read_table(paths["taxonomy"]), "catalog_id")
    if len(genomes) != 584 or len(primary) != 294 or set(taxonomy) != genomes or not primary <= genomes:
        raise ValueError("Catalogue or primary manifest identity mismatch")
    counts, definitions, plant = read_workbook(paths["workbook"], genomes)
    support = json.loads(paths["sets"].read_text())
    if len(support) != 339 or not primary <= support.keys():
        raise ValueError("Expected complete model-support contract")
    for genome, record in support.items():
        if genome != record["catalog_id"] or not isinstance(record.get("requirements"), list):
            raise ValueError("Malformed strict requirement record")
    # Assign every trait before loading disease-differential labels.
    calls = {g: trait_values({gene: vector[g] for gene, vector in counts.items()}, plant,
                            support[g]["requirements"] if g in support else None) for g in genomes}
    ancom = keyed(read_table(paths["ancom"], delimiter=","), "taxon")
    if set(ancom) - genomes:
        raise ValueError("Unexpected ANCOM genome identifier")
    long_rows, wide_rows, evidence_rows = [], [], []
    relevant = sorted(plant | CORE_MARKERS)
    for genome in sorted(genomes):
        meta = taxonomy[genome]
        overlay = direction(ancom.get(genome))
        common = dict(catalog_id=genome, primary=genome in primary,
                      phylum=meta["phylum"], family=meta["family"],
                      enterobacteriaceae_overlay=meta["family"] == "Enterobacteriaceae",
                      ancom_direction=overlay)
        wide_rows.append(dict(common, **{trait: state(value) for trait, value in calls[genome].items()}))
        for trait, label, rule, caveat in TRAITS:
            long_rows.append(dict(common, trait_id=trait, trait_label=label, state=state(calls[genome][trait])))
        for gene in relevant:
            count = counts[gene][genome]
            evidence_rows.append(dict(catalog_id=genome, gene_id=gene,
                                     count="" if count is None else count,
                                     evidence_state=state(has(count))))
        if genome in support:
            for acid in sorted(ACIDS):
                evidence_rows.append(dict(catalog_id=genome, gene_id="strict_requirement:" + acid,
                    count=int(acid in support[genome]["requirements"]),
                    evidence_state=state(acid in support[genome]["requirements"])))
    groups = collections.defaultdict(list)
    for row in long_rows:
        if row["primary"]:
            groups[(row["trait_id"], row["trait_label"], row["ancom_direction"])].append(row["state"])
    summary = []
    for (trait, label, category), states in sorted(groups.items()):
        n = len(states); positive = states.count("positive"); unknown = states.count("unknown")
        summary.append(dict(trait_id=trait, trait_label=label, ancom_direction=category,
                            n_group=n, n_positive=positive, n_not_detected=states.count("evidence_not_detected"),
                            n_unknown=unknown, fraction_positive_all=positive/n,
                            fraction_positive_known=positive/(n-unknown) if n > unknown else "",
                            denominator="all primary MAGs in the ANCOM category; overlapping traits"))
    out.mkdir(parents=True, exist_ok=True)
    tables = {"genome_traits_long.tsv": long_rows, "genome_traits_wide.tsv": wide_rows,
              "primary294_traits_long.tsv": [r for r in long_rows if r["primary"]],
              "primary294_traits_wide.tsv": [r for r in wide_rows if r["primary"]],
              "trait_group_summary.tsv": summary, "gene_evidence_counts.tsv": evidence_rows,
              "marker_definitions.tsv": [d for d in definitions if d["gene_id"] in relevant],
              "trait_definitions.tsv": [dict(trait_id=t, trait_label=l, rule=r, caveat=c,
                  marker_ids=";".join(sorted(plant)) if t == "plant_backbone_cazyme" else "",
                  sources_file="TRAIT_METHODS_v1.md") for t,l,r,c in TRAITS]}
    for name, rows in tables.items():
        write_table(out / name, rows)
    if {key: sha(path) for key, path in paths.items()} != hashes:
        raise ValueError("Input changed during trait build")
    gate = dict(complete=True, stage="complete", created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                catalogue_n=len(genomes), primary_n=len(primary), primary_model_n=len(primary), model_n=len(support), trait_n=len(TRAITS),
                primary_rows=len(tables["primary294_traits_long.tsv"]),
                primary_ids_sha256=hashlib.sha256("\n".join(sorted(primary)).encode()).hexdigest(),
                primary_unknown_n=sum(r["state"] == "unknown" for r in tables["primary294_traits_long.tsv"]),
                overlapping_traits=True, significance_selection=False, statistical_tests_run=False,
                input_paths={k:str(v) for k,v in paths.items()}, input_sha256=hashes,
                output_sha256={name:sha(out/name) for name in tables}, plant_marker_ids=sorted(plant),
                interpretation="annotation/model evidence; not exclusive guilds, pathogens, flux or observed cross-feeding")
    (out / "analysis_complete.json").write_text(json.dumps(gate, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k:gate[k] for k in ["complete","catalogue_n","primary_n","model_n","trait_n","primary_unknown_n"]}))


if __name__ == "__main__":
    main()
