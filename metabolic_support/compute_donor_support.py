#!/usr/bin/env python3
"""Frozen, donor-replicated topological support analysis for COVID Figure 5.

Only consumes validated graph outputs. Never reconstructs, gap-fills, changes
the graph whitelist, selects models by ANCOM, or edits manuscript files.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import scipy
from scipy.stats import mannwhitneyu, wilcoxon

PROJECT = Path(__file__).resolve().parents[2]
ROOT = PROJECT / "result2/n76/metabolic_support"
SEED = 42
N_RANDOM = 1000
N_BOOT = 2000


def require(value, message):
    if not value:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def truth(values):
    return values.astype(str).str.lower().eq("true")


def write_table(path, rows, columns=None):
    frame = pd.DataFrame(rows, columns=columns)
    frame.to_csv(path, sep="\t", index=False, na_rep="NA", float_format="%.12g")
    return frame


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def bh(p):
    p = np.asarray(p, dtype=float)
    out = np.full(len(p), np.nan)
    ok = np.where(np.isfinite(p))[0]
    order = ok[np.argsort(p[ok])]
    if len(order):
        q = p[order] * len(order) / np.arange(1, len(order) + 1)
        out[order] = np.minimum(1, np.minimum.accumulate(q[::-1])[::-1])
    return out


class SupportKernel:
    """Fixed recipient-edge denominator, vectorised over random supplier sets."""

    def __init__(self, identifiers, support, recipients, kind="strict"):
        self.ids = list(identifiers)
        self.index = {k: i for i, k in enumerate(self.ids)}
        self.recipients = list(recipients)
        require(bool(recipients), "Cannot score an empty recipient panel")
        keymap = {"strict": ("requirements", "providers"),
                  "topology_only": ("topology_only_requirements", "topology_only_providers"),
                  "stock": ("stock_requirements", "stock_providers")}
        req_key, pro_key = keymap[kind]
        self.metabolites = sorted(set().union(*[
            set(support[k][req_key]) | set(support[k][pro_key]) for k in self.ids]))
        mi = {m: i for i, m in enumerate(self.metabolites)}
        self.P = np.zeros((len(self.ids), len(mi)), dtype=np.int16)
        for k in self.ids:
            self.P[self.index[k], [mi[m] for m in support[k][pro_key]]] = 1
        edge_r, edge_m, edge_local_r, weights = [], [], [], []
        for ri, k in enumerate(recipients):
            demands = sorted(set(support[k][req_key]))
            require(bool(demands), f"Empty recipient requirement set: {k}")
            for m in demands:
                edge_r.append(self.index[k])
                edge_m.append(mi[m])
                edge_local_r.append(ri)
                weights.append(1 / (len(recipients) * len(demands)))
        self.er = np.asarray(edge_r)
        self.em = np.asarray(edge_m)
        self.local_r = np.asarray(edge_local_r)
        self.weights = np.asarray(weights)
        self.self_provides = self.P[self.er, self.em]
        require(np.isclose(self.weights.sum(), 1), "Support weights do not sum to 1")

    def counts(self, suppliers):
        sets = np.atleast_2d(suppliers).astype(np.int16)
        counts = (sets @ self.P)[:, self.em]
        counts -= sets[:, self.er] * self.self_provides
        require(np.all(counts >= 0), "Negative supplier count")
        return counts

    def scores(self, suppliers):
        counts = self.counts(suppliers)
        return {
            "mss_reference": np.round((counts > 0) @ self.weights, 12),
            "redundancy_mean": np.round(counts @ self.weights, 12),
            "redundancy_ge2": np.round((counts >= 2) @ self.weights, 12),
        }

    def detail(self, suppliers, abundance, families):
        scores = {key: float(value[0]) for key, value in self.scores(suppliers).items()}
        counts = self.counts(suppliers)[0]
        supported = counts > 0
        recipient_support = np.bincount(self.local_r, weights=supported * self.weights,
                                       minlength=len(self.recipients)) * len(self.recipients)
        resident = np.asarray([suppliers[self.index[k]] for k in self.recipients], bool)
        scores["mss_resident"] = float(recipient_support[resident].mean()) if resident.any() else np.nan
        recipient_abundance = np.asarray([abundance[self.index[k]] for k in self.recipients])
        weights = recipient_abundance * resident
        scores["mss_resident_abundance_weighted"] = (
            float(np.average(recipient_support, weights=weights)) if weights.sum() else np.nan)
        family_counts = np.zeros(len(self.er), dtype=int)
        for family in sorted(set(families)):
            family_set = suppliers & (families == family)
            family_counts += self.counts(family_set)[0] > 0
        supported_weight = float(self.weights[supported].sum())
        scores["single_family_share"] = (float(self.weights[(family_counts == 1) & supported].sum()) /
                                         supported_weight if supported_weight else np.nan)
        return scores


def strata_for(manifest, qc, abundance, pool, include_phylum=False):
    """Cutpoints defined label-blind in the dated implementation lock."""
    genes = pd.to_numeric(qc.loc[manifest.index, "model_genes"]).to_numpy()
    primary = truth(manifest["support_eligible_primary"]).to_numpy()
    gene_cut = float(np.median(genes[primary]))
    cuts = np.quantile(abundance[pool], [1 / 3, 2 / 3]) if pool.any() else [0., 0.]
    comp = pd.to_numeric(manifest.completeness).to_numpy()
    con = pd.to_numeric(manifest.contamination).to_numpy()
    strata = []
    for i in range(len(manifest)):
        tokens = [str(int(comp[i] > 95)), str(int(con[i] >= 2)),
                  str(int(genes[i] >= gene_cut)), str(int(np.searchsorted(cuts, abundance[i], side="left")))]
        if include_phylum:
            tokens.append(str(manifest.iloc[i].phylum))
        strata.append("|".join(tokens))
    return np.asarray(strata), {"gene_median": gene_cut, "abundance_tertile_1": float(cuts[0]),
                               "abundance_tertile_2": float(cuts[1]), "exact_phylum": include_phylum}


def sample_matched(target, pool, strata, rng, n=N_RANDOM):
    """Independent uniform within-stratum sets, target identities eligible."""
    require(np.all(~target | pool), "Targets must be eligible under matched null")
    result = np.zeros((n, len(target)), dtype=bool)
    movable = 0
    audit = []
    for label in sorted(set(strata[target])):
        candidates = np.flatnonzero(pool & (strata == label))
        k = int(np.sum(target & (strata == label)))
        require(k <= len(candidates), "Unattainable exact stratum match")
        audit.append({"stratum": label, "target_n": k, "candidate_n": len(candidates)})
        if k == len(candidates):
            result[:, candidates] = True
        else:
            movable += k
            # Independent uniform subsets using iid random ranks; no donor labels.
            ranks = rng.random((n, len(candidates)))
            selected = np.argpartition(ranks, k - 1, axis=1)[:, :k]
            result[np.arange(n)[:, None], candidates[selected]] = True
    require(np.all(result.sum(axis=1) == target.sum()), "Matched cardinality drift")
    return result, movable, audit


def bootstrap_mean(values, groups, rng, n=N_BOOT):
    estimates = np.zeros(n)
    weights = 0
    for group in sorted(set(groups)):
        x = np.asarray(values)[np.asarray(groups) == group]
        x = x[np.isfinite(x)]
        if len(x):
            estimates += rng.choice(x, (n, len(x)), replace=True).sum(axis=1)
            weights += len(x)
    return estimates / weights if weights else np.full(n, np.nan)


def contrast(frame, outcome, scenario, comparison, rng):
    f = frame.loc[np.isfinite(frame[outcome])].copy()
    case = f.loc[f.group == "Case"]
    control = f.loc[f.group == "Control"]
    if comparison != "all_controls":
        control = control.loc[control.source == comparison]
    if not len(case) or not len(control):
        return None
    a, b = case[outcome].to_numpy(), control[outcome].to_numpy()
    boot = bootstrap_mean(a, case.source, rng) - bootstrap_mean(b, control.source, rng)
    return {"scenario": scenario, "outcome": outcome, "comparison": comparison,
            "case_n": len(a), "control_n": len(b), "case_mean": float(a.mean()),
            "control_mean": float(b.mean()), "difference": float(a.mean() - b.mean()),
            "ci_low": float(np.quantile(boot, .025)), "ci_high": float(np.quantile(boot, .975)),
            "p": float(mannwhitneyu(a, b, alternative="two-sided", method="asymptotic").pvalue),
            "case_median": float(np.median(a)), "control_median": float(np.median(b)),
            "interpretation": "observational_source_confounded_donor_contrast"}


def load_inputs():
    paths = {
        "manifest": ROOT / "MAG_manifest_sensitivity_bacteria339.tsv",
        "sets": ROOT / "graph/model_support_sets.json",
        "graph_qc": ROOT / "graph/model_graph_qc.tsv",
        "graph_audit": ROOT / "graph/extraction_audit.json",
        "reconstruction_status": ROOT / "model_reconstruction_status_339.tsv",
        "reconstruction_summary": ROOT / "reconstruction_solver_summary.json",
        "whitelist": ROOT / "graph/metabolite_whitelist.tsv",
        "lock": ROOT / "MATCHING_IMPLEMENTATION_LOCK_20260921.md",
        "graph_lock": ROOT / "graph/extraction_lock.json",
        "protocol": ROOT / "ANALYSIS_PROTOCOL_v1.md",
        "sensitivity_note": ROOT / "SENSITIVITY_IMPLEMENTATION_NOTE_20260921.md",
        "abundance": PROJECT / "result2/inputfiles/new/MAG_relative_abundance_percent.tsv",
        "metadata": PROJECT / "result2/inputfiles/new/sample_metadata_76.tsv",
        "ancom": PROJECT / "result2/n76/ANCOMBC2_MAG_results_n76.csv",
    }
    for label, path in paths.items():
        require(path.is_file(), f"Missing {label}: {path}")
    manifest = pd.read_csv(paths["manifest"], sep="\t", keep_default_na=False).set_index("catalog_id")
    require(len(manifest) == 339 and manifest.index.is_unique, "Expected339 unique models")
    require(truth(manifest.support_eligible_primary).sum() == 294, "Primary294 drift")
    support = json.loads(paths["sets"].read_text())
    require(set(support) == set(manifest.index), "Graph outputs incomplete: need all339 before scoring")
    graph_audit = json.loads(paths["graph_audit"].read_text())
    require(graph_audit["stage"] == "full" and graph_audit["complete"] and
            graph_audit["completed_models"] == 339 and not graph_audit["failed_models"],
            "Full extraction completion gate failed")
    status = pd.read_csv(paths["reconstruction_status"], sep="\t").set_index("catalog_id")
    require(status.index.is_unique and set(status.index) == set(manifest.index), "Model status identity mismatch")
    require(status.status.eq("complete_hash_verified").all(), "Unverified reconstructed models")
    require(json.loads(paths["reconstruction_summary"].read_text())["status"] == "complete",
            "Reconstruction summary is only a partial snapshot")
    graph_lock = json.loads(paths["graph_lock"].read_text())
    require(graph_lock["whitelist_sha256"] == sha(paths["whitelist"]), "Frozen whitelist hash drift")
    whitelist = set(pd.read_csv(paths["whitelist"], sep="\t").metabolite_id)
    for cid, item in support.items():
        require(item["model_sha256"] == status.loc[cid, "model_sha256"], f"Graph/model hash drift: {cid}")
        require(item["graph_method"] == graph_lock["method"], f"Graph method mismatch: {cid}")
        require(item["whitelist_sha256"] == graph_lock["whitelist_sha256"], f"Model whitelist drift: {cid}")
        for field in ["requirements", "providers", "topology_only_requirements", "topology_only_providers",
                      "stock_requirements", "stock_providers"]:
            require(field in item, f"Missing planned graph sensitivity {field}: {cid}")
            require(set(item[field]) <= whitelist, f"Unverified metabolite in {field}: {cid}")
            require(len(item[field]) == len(set(item[field])), f"Duplicated {field}: {cid}")
    qc = pd.read_csv(paths["graph_qc"], sep="\t", keep_default_na=False).set_index("catalog_id")
    require(qc.index.is_unique and set(qc.index) == set(manifest.index), "Graph QC identity mismatch")
    require(truth(qc.technical_qc_pass).all(), "Model QC failure: explicit review required")
    require(qc.graph_status.eq("complete").all(), "Incomplete graphs")
    metadata = pd.read_csv(paths["metadata"], sep="\t", keep_default_na=False).set_index("sample")
    require(len(metadata) == 76 and metadata.index.is_unique, "Expected76 unique donors")
    require(metadata.group.value_counts().to_dict() == {"Case": 38, "Control": 38}, "Arm counts drift")
    abundance_all = pd.read_csv(paths["abundance"], sep="\t").set_index("MAG")
    require(len(abundance_all) == 584 and abundance_all.index.is_unique, "Catalogue584 identity drift")
    require(set(metadata.index) == set(abundance_all.columns), "Donor abundance identity mismatch")
    abundance = abundance_all.loc[manifest.index, metadata.index].to_numpy(dtype=float)
    require(np.all(np.isfinite(abundance)) and np.all(abundance >= 0), "Bad abundance values")
    ancom = pd.read_csv(paths["ancom"]).set_index("taxon")
    require(ancom.index.is_unique, "Duplicated ANCOM identities")
    require(set(manifest.index) <= set(ancom.index), "Models missing from ANCOM annotation")
    return paths, manifest, support, qc, metadata, abundance_all, abundance, ancom


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "support")
    args = parser.parse_args()
    out = args.output.resolve()
    paths, manifest, support, qc, metadata, abundance_all, abundance, ancom = load_inputs()
    out.mkdir(parents=True, exist_ok=True)
    require(not (out / "analysis_complete.json").exists(), "Completed analysis exists; do not silently overwrite")
    input_hashes = {name: sha(path) for name, path in paths.items()}
    input_hashes["script"] = sha(__file__)
    write_json(out / "analysis_inputs.json", {"sha256": input_hashes, "python": sys.version,
                                              "numpy": np.__version__, "pandas": pd.__version__,
                                              "scipy": scipy.__version__,
                                              "platform": platform.platform()})
    ids = list(manifest.index)
    primary = truth(manifest.support_eligible_primary).to_numpy()
    detected = abundance >= .01
    families = np.asarray([value if value else f"unclassified__{cid}"
                           for cid, value in zip(ids, manifest.family)])
    enterobacter = families == "Enterobacteriaceae"
    recipients = [k for i, k in enumerate(ids) if primary[i] and detected[i].sum() >= 16
                  and support[k]["requirements"]]
    whitelist_n = len(pd.read_csv(paths["whitelist"], sep="\t"))
    gate = {"feasible": bool(recipients), "recipient_n": len(recipients),
            "validated_model_n": len(ids), "primary_model_n": int(primary.sum()),
            "models_with_supported_requirements": sum(bool(support[k]["requirements"]) for k in ids),
            "models_with_supported_provision": sum(bool(support[k]["providers"]) for k in ids),
            "whitelist_metabolite_n": whitelist_n, "sha256": input_hashes,
            "feasibility_reason": "nonempty_fixed_transport_supported_recipient_panel" if recipients else
                                  "no_primary_model_with_nonempty_requirements_and_prevalence_at_least16",
            "scope": "transport_grounded_topological_support_not_auxotrophy_flux_or_measured_crossfeeding"}
    recipient_audit = []
    for i, k in enumerate(ids):
        reasons = []
        if not primary[i]: reasons.append("not_primary_GUNC_pass")
        if detected[i].sum() < 16: reasons.append("pooled_detection_below16")
        if not support[k]["requirements"]: reasons.append("empty_strict_requirements")
        recipient_audit.append({"catalog_id": k, "primary": bool(primary[i]),
                                "pooled_prevalence": int(detected[i].sum()),
                                "n_requirements": len(support[k]["requirements"]),
                                "included": not reasons, "exclusion_reason": ";".join(reasons)})
    write_table(out / "recipient_audit.tsv", recipient_audit)
    if not recipients:
        write_json(out / "analysis_complete.json", {**gate, "stage": "scientific_feasibility_failure"})
        print(json.dumps(gate, indent=2))
        return
    kernel = SupportKernel(ids, support, recipients)
    primary_pool = primary.copy()
    ceiling_counts = kernel.counts(primary_pool)[0]
    ceiling_scores = {key: float(value[0]) for key, value in kernel.scores(primary_pool).items()}
    write_table(out / "reference_requirement_ceiling.tsv", [
        {"recipient_id": ids[ri], "metabolite_id": kernel.metabolites[mi],
         "primary_catalogue_supplier_n": int(count), "recipient_weighted_score_weight": float(weight)}
        for ri, mi, count, weight in zip(kernel.er, kernel.em, ceiling_counts, kernel.weights)])
    gate["primary_catalogue_support_ceiling"] = ceiling_scores
    gate["reference_requirement_pairs"] = len(kernel.er)
    gate["reference_pairs_with_no_primary_provider"] = int((ceiling_counts == 0).sum())
    pooled_mean = abundance.mean(axis=1)
    strata, cuts = strata_for(manifest, qc, pooled_mean, primary_pool)
    cuts["gene_count_field"] = "model_graph_qc.model_genes"
    write_json(out / "matching_cutpoints.json", cuts)
    write_table(out / "supplier_matching_strata.tsv", [
        {"catalog_id": k, "primary": bool(primary[i]), "stratum": strata[i],
         "pooled_mean_abundance_percent": pooled_mean[i], "model_genes": qc.loc[k, "model_genes"]}
        for i, k in enumerate(ids)])
    raw_totals = abundance_all.loc[:, metadata.index].sum(axis=0).to_numpy()
    scenarios = [("primary", .01, primary), ("threshold_0.001", .001, primary),
                 ("threshold_0.1", .1, primary), ("without_Enterobacteriaceae", .01, primary & ~enterobacter),
                 ("quality339_fixed_recipients", .01, np.ones(len(ids), bool))]
    donor_rows, metabolite_rows, matching_audit, interventions = [], [], [], []
    primary_masks = []
    null_dir = out / "null_draws"
    null_dir.mkdir(exist_ok=True)
    for si, (sample, meta) in enumerate(metadata.iterrows()):
        print(f"Donor{si + 1}/76 {sample}", flush=True)
        rng = np.random.default_rng(np.random.SeedSequence([SEED, si, 1]))
        for scenario, threshold, eligible in scenarios:
            suppliers = (abundance[:, si] >= threshold) & eligible
            metrics = kernel.detail(suppliers, abundance[:, si], families)
            modelled_abundance = float(abundance[eligible, si].sum())
            row = {"sample": sample, "group": meta.group, "source": meta.source,
                   "scenario": scenario, "threshold": threshold, "recipient_n": len(recipients),
                   "supplier_n": int(suppliers.sum()), **metrics,
                   "catalogue_abundance": raw_totals[si], "modelled_abundance": modelled_abundance,
                   "coverage": modelled_abundance / raw_totals[si]}
            if scenario == "primary":
                primary_masks.append(suppliers)
                random_sets, movable, audit = sample_matched(suppliers, primary_pool, strata, rng)
                for item in audit: matching_audit.append({"sample": sample, "task": "community_null", **item})
                random_scores = kernel.scores(random_sets)["mss_reference"]
                obs = row["mss_reference"]
                null_mean = float(random_scores.mean())
                null_sd = float(random_scores.std(ddof=1))
                row.update({"matched_null_mean": null_mean, "matched_null_sd": null_sd,
                            "matched_null_lo": float(np.quantile(random_scores, .025)),
                            "matched_null_hi": float(np.quantile(random_scores, .975)),
                            "excess_support": float(np.round(obs - null_mean, 12)),
                            "null_p": (1 + int(np.sum(np.round(np.abs(random_scores - null_mean), 12) >=
                                                      np.round(abs(obs - null_mean), 12)))) /
                                      (1 + N_RANDOM),
                            "random_movable_supplier_n": movable,
                            "identical_random_fraction": float(np.mean(np.all(random_sets == suppliers, axis=1))),
                            "null_status": "degenerate_null" if movable == 0 else "matched"})
                np.savez_compressed(null_dir / f"community_{si + 1:02d}.npz",
                                    sets=np.packbits(random_sets, axis=1), scores=random_scores,
                                    catalog_ids=np.asarray(ids), sample=sample)
                counts = kernel.counts(suppliers)[0]
                for met in sorted(set(kernel.em)):
                    mask = kernel.em == met
                    metabolite_rows.append({"sample": sample, "group": meta.group, "source": meta.source,
                                            "metabolite_id": kernel.metabolites[met],
                                            "recipient_n": int(mask.sum()),
                                            "support_fraction": float(np.mean(counts[mask] > 0)),
                                            "redundancy_ge2": float(np.mean(counts[mask] >= 2))})
            donor_rows.append(row)

        # Same locked panel; abundance affects matching, never recipient membership.
        observed = primary_masks[-1]
        baseline = float(kernel.scores(observed)["mss_reference"][0])
        depleted = (truth(ancom.loc[ids, "diff_GroupCase"]).to_numpy() &
                    (pd.to_numeric(ancom.loc[ids, "lfc_GroupCase"]).to_numpy() < 0))
        for ii, operation in enumerate(["remove_depleted", "add_depleted"]):
            removal = operation.startswith("remove")
            pool = primary & (observed if removal else ~observed)
            target = pool & depleted
            irng = np.random.default_rng(np.random.SeedSequence([SEED, si, ii + 2]))
            row = {"sample": sample, "group": meta.group, "source": meta.source,
                   "intervention": operation, "target_label": "ANCOM_depleted",
                   "target_n": int(target.sum()), "candidate_n": int(pool.sum()),
                   "baseline_score": baseline, "recipient_n": len(recipients),
                   "delta": np.nan, "random_delta_mean": np.nan, "random_delta_lo": np.nan,
                   "random_delta_hi": np.nan, "excess_delta": np.nan, "n_random_valid": 0,
                   "random_movable_target_n": 0, "identical_random_fraction": np.nan}
            if not target.any():
                row.update({"match_status": "no_targets", "n_random_valid": 0})
                interventions.append(row)
                continue
            match_abundance = abundance[:, si] if removal else pooled_mean
            istrata, icuts = strata_for(manifest, qc, match_abundance, pool, include_phylum=True)
            random_targets, movable, audit = sample_matched(target, pool, istrata, irng)
            for item in audit: matching_audit.append({"sample": sample, "task": operation, **icuts, **item})
            altered = observed & ~target if removal else observed | target
            random_altered = observed[None, :] & ~random_targets if removal else observed[None, :] | random_targets
            delta = float(np.round(kernel.scores(altered)["mss_reference"][0] - baseline, 12))
            null_delta = np.round(kernel.scores(random_altered)["mss_reference"] - baseline, 12)
            require(delta <= 1e-12 if removal else delta >= -1e-12, "Union intervention monotonicity violated")
            require(np.all(null_delta <= 1e-12) if removal else np.all(null_delta >= -1e-12),
                    "Random union intervention monotonicity violated")
            row.update({"delta": delta, "random_delta_mean": float(null_delta.mean()),
                        "random_delta_lo": float(np.quantile(null_delta, .025)),
                        "random_delta_hi": float(np.quantile(null_delta, .975)),
                        "excess_delta": float(np.round(delta - float(null_delta.mean()), 12)), "n_random_valid": N_RANDOM,
                        "random_movable_target_n": movable,
                        "identical_random_fraction": float(np.mean(np.all(random_targets == target, axis=1))),
                        "match_status": "degenerate_null" if movable == 0 else "matched"})
            np.savez_compressed(null_dir / f"{operation}_{si + 1:02d}.npz",
                                sets=np.packbits(random_targets, axis=1), delta=null_delta,
                                catalog_ids=np.asarray(ids), sample=sample)
            interventions.append(row)

    # Explicit method sensitivities. The stock comparison uses the SAME common
    # recipients on both sides; stock has different requirements by definition.
    # Never compare changed recipient denominators as if only suppliers changed.
    expanded = [k for i, k in enumerate(ids) if detected[i].sum() >= 16 and support[k]["requirements"]]
    stock_common = [k for k in recipients if support[k]["stock_requirements"]]
    additional = [("topology_only_same_recipients", kernel.recipients, "topology_only", primary),
                  ("quality339_expanded_recipients", expanded, "strict", np.ones(len(ids), bool))]
    if stock_common:
        additional += [("primary_stock_common_recipients", stock_common, "strict", primary),
                       ("stock_phyloMint_common_recipients", stock_common, "stock", primary)]
    for scenario, these_recipients, kind, pool in additional:
        skernel = SupportKernel(ids, support, these_recipients, kind=kind)
        for si, (sample, meta) in enumerate(metadata.iterrows()):
            suppliers = detected[:, si] & pool
            modelled = float(abundance[pool, si].sum())
            donor_rows.append({"sample": sample, "group": meta.group, "source": meta.source,
                               "scenario": scenario, "threshold": .01, "recipient_n": len(these_recipients),
                               "supplier_n": int(suppliers.sum()),
                               **skernel.detail(suppliers, abundance[:, si], families),
                               "catalogue_abundance": raw_totals[si], "modelled_abundance": modelled,
                               "coverage": modelled / raw_totals[si]})
    donor = write_table(out / "donor_support.tsv", donor_rows)
    metabolite = write_table(out / "metabolite_support.tsv", metabolite_rows)
    int_frame = write_table(out / "intervention_donors.tsv", interventions)
    write_table(out / "matching_stratum_audit.tsv", matching_audit)
    contrasts = []
    crng = np.random.default_rng(np.random.SeedSequence([SEED, 999]))
    outcomes = ["mss_reference", "redundancy_mean", "redundancy_ge2", "single_family_share",
                "mss_resident", "mss_resident_abundance_weighted"]
    control_sources = sorted(metadata.loc[metadata.group == "Control", "source"].unique())
    for scenario in donor.scenario.unique():
        frame = donor.loc[donor.scenario == scenario]
        for outcome in outcomes + (["excess_support"] if scenario == "primary" else []):
            for comparison in (["all_controls"] + control_sources if scenario == "primary" else ["all_controls"]):
                item = contrast(frame, outcome, scenario, comparison, crng)
                if item: contrasts.append(item)
    prim_frame = donor.loc[donor.scenario == "primary"].copy()
    low = max(prim_frame.loc[prim_frame.group == g, "coverage"].min() for g in ["Case", "Control"])
    high = min(prim_frame.loc[prim_frame.group == g, "coverage"].max() for g in ["Case", "Control"])
    overlap = prim_frame.loc[prim_frame.coverage.between(low, high)].copy()
    write_table(out / "coverage_overlap_donors.tsv", overlap)
    for outcome in outcomes + ["excess_support"]:
        item = contrast(overlap, outcome, "coverage_overlap", "all_controls", crng)
        if item: contrasts.append(item)
    write_table(out / "contrast_summary.tsv", contrasts)
    met_contrasts = []
    for met in metabolite.metabolite_id.unique():
        item = contrast(metabolite.loc[metabolite.metabolite_id == met], "support_fraction",
                        "primary_metabolite", "all_controls", crng)
        if item: met_contrasts.append({"metabolite_id": met, **item})
    met_test = pd.DataFrame(met_contrasts)
    met_test["q"] = bh(met_test.p)
    write_table(out / "metabolite_contrasts.tsv", met_test)
    int_summary = []
    for operation in int_frame.intervention.unique():
        for group in ["Case", "Control", "All"]:
            frame = int_frame.loc[(int_frame.intervention == operation) & int_frame.match_status.eq("matched")]
            if group != "All": frame = frame.loc[frame.group == group]
            if not len(frame): continue
            for outcome in ["delta", "random_delta_mean", "excess_delta"]:
                f = frame.loc[np.isfinite(frame[outcome])]
                if not len(f): continue
                values = np.round(f[outcome].to_numpy(), 12)
                boot = bootstrap_mean(values, f.source, crng)
                # Raw union deltas are directionally automatic, not specificity
                # tests. Only targeted-minus-matched-null deltas receive a test.
                p = (float(wilcoxon(values, alternative="two-sided", method="asymptotic").pvalue)
                     if np.any(values != 0) else 1.) if outcome == "excess_delta" else np.nan
                int_summary.append({"intervention": operation, "group": group, "outcome": outcome,
                                    "donor_n": len(values), "mean": float(values.mean()),
                                    "ci_low": float(np.quantile(boot, .025)), "ci_high": float(np.quantile(boot, .975)),
                                    "p": p, "status": "exploratory_donor_level_matched_null"})
    int_summary = pd.DataFrame(int_summary, columns=["intervention", "group", "outcome", "donor_n", "mean",
                                                     "ci_low", "ci_high", "p", "status"])
    int_summary["q"] = bh(int_summary.p)
    write_table(out / "intervention_summary.tsv", int_summary)
    gate.update({"stage": "complete", "donor_n": len(metadata), "random_sets_per_donor": N_RANDOM,
                 "bootstrap_iterations": N_BOOT, "seed": SEED,
                 "recipient_metabolite_n": len(set(kernel.em)),
                 "stock_common_recipient_n": len(stock_common), "expanded_recipient_n": len(expanded),
                 "coverage_overlap": {"low": float(low), "high": float(high),
                                      "case_n": int((overlap.group == "Case").sum()),
                                      "control_n": int((overlap.group == "Control").sum())},
                 "output_sha256": {p.name: sha(p) for p in sorted(out.glob("*.tsv"))}})
    write_json(out / "analysis_complete.json", gate)
    print(json.dumps(gate, indent=2), flush=True)


if __name__ == "__main__":
    main()
