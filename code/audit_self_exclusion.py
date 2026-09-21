"""Audit of the background self-exclusion correction: baseline, corrected, and the difference.

`control_sets` used to build its candidate pool as the gene universe minus the frozen exclusion set,
which does not contain the 200 background genes. Matching is by absolute distance in reference rank,
so a background gene sat at distance zero from itself and was taken as its own control. The
analysed genes were never affected, because `group_table` adds the scored list to the exclusion set.
The background scores define the reference bands, and the bands are the denominator of the
normalised contrasts and the criterion of the outside-band frames, so the defect reached the figure
through the bands alone.

This script rebuilds the rank matrices from the primary inputs and scores the same 200 background
genes twice in one process: with the selection rule as released, and with the gene excluded from its
own controls. The baseline pass must reproduce the deposited background table before the difference
is read as the effect of self-exclusion, and the script fails if it does not.

Outputs (in the output directory):
  input_manifest.tsv                     every file the loaders open, with size and SHA-256
  baseline_reconstruction_report.json    cohort checks and the reconciliation with the release
  self_exclusion_audit.tsv               gene by organ: self-inclusion before and after, what moved
  control_assignments_baseline.tsv.gz    the ordered controls and their rank gaps, as released
  control_assignments_corrected.tsv.gz   the same after the correction
  background_score_changes.tsv.gz        per profile: released score, corrected score, difference
  analysed_gene_control_sets.tsv          the 67 analysed genes: their control sets under both rules

Usage: python code/audit_self_exclusion.py <pipeline.conf> --tables tables --output-dir DIR
"""
import argparse
import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_layer_scores as B
import level_matched_scores as L
from table_io import read_table

ORGANS, LAYERS = B.ORGANS, B.LAYERS


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def inventory(conf):
    """One row per file the loaders actually open, addressed by a portable logical path.

    The enumeration follows the selection rules of the calculation rather than a directory glob:
    the healthy layer is selected through the pseudobulk manifest, the malignant layer through the
    study-selection table, and each selected dataset contributes its matrix, its gene axis and its
    group table. A required component that selects nothing is an error, not an empty digest. Paths
    are logical: the local directory of an input is configuration, not provenance, and an absolute
    path from the preparing machine belongs in no deposited file.
    """
    rows = []

    def add(role, path, logical):
        path = Path(path)
        if not path.is_file():
            raise SystemExit(f"{role}: required input missing ({logical})")
        rows.append({"role": role, "logical_path": logical, "bytes": path.stat().st_size,
                     "sha256": digest(path)})

    add("gene_sets", conf["MEMBERSHIP_TSV"], "tables/IKM_gene_membership_by_group.tsv")
    add("candidate_panel", conf["PANEL_TSV"], "tables/panoptosis_panel_2026-09-18.tsv")

    depmap = Path(conf["DEPMAP_DIR"])
    for role, name in [("depmap_expression", "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"),
                       ("depmap_models", "Model.csv"),
                       ("depmap_nextgen", "DepMap_NextGen_3D_core_files.zip")]:
        add(role, depmap / name, f"depmap/{name}")

    pb = Path(conf["PSEUDOBULK_DIR"])
    healthy_manifest = pb / "tabula" / "manifest.csv"
    add("healthy_selection", healthy_manifest, "pseudobulk/tabula/manifest.csv")
    selected = pd.read_csv(healthy_manifest)
    if not len(selected):
        raise SystemExit("healthy_selection: the manifest selects no organ")
    for prefix in selected["prefix"].astype(str):
        for role, suffix in [("healthy_matrix", "_counts.npy"), ("healthy_gene_axis", "_genes.csv"),
                             ("healthy_groups", "_groups.tsv")]:
            add(role, pb / "tabula" / f"{prefix}{suffix}", f"pseudobulk/tabula/{prefix}{suffix}")

    studies = Path(conf["KANG_STUDIES_CSV"])
    add("malignant_selection", studies, f"selection/{studies.name}")
    datasets = pd.read_csv(studies)["Dataset"].astype(str).unique()
    if not len(datasets):
        raise SystemExit("malignant_selection: the study table selects no dataset")
    for dataset in datasets:
        # patient_layer opens <dataset>_mean.npy, not _counts.npy: the malignant pseudobulks are
        # per-sample means, and globbing the healthy suffix here selected nothing at all.
        for role, suffix in [("malignant_matrix", "_mean.npy"), ("malignant_gene_axis", "_genes.csv"),
                             ("malignant_groups", "_groups.tsv")]:
            add(role, pb / "kang" / f"{dataset}{suffix}", f"pseudobulk/kang/{dataset}{suffix}")

    frame = pd.DataFrame(rows)
    for role, least in [("healthy_matrix", 6), ("malignant_matrix", 1)]:
        found = int((frame.role == role).sum())
        if found < least:
            raise SystemExit(f"{role}: {found} files selected, at least {least} required")
    return frame


def released_pool(tables):
    """The 200 background genes in the order the released table carries them."""
    background = read_table(Path(tables) / "IKM_random_gene_scores_by_sample.csv.gz")
    genes = list(dict.fromkeys(background["gene"]))
    assert len(genes) == 200, len(genes)
    return genes, background


def nearest_without_self(ref_rank, genes, gene_index, exclude, n_ctrl=400):
    """The released rule: the pool is the universe minus `exclude`, the gene itself included."""
    pool = np.array(sorted(i for g, i in gene_index.items() if g not in exclude))
    names = {i: g for g, i in gene_index.items()}
    pool_names = np.array([names[i] for i in pool])
    ref_pool = ref_rank[pool]
    out = {}
    for g in genes:
        j = gene_index[g]
        d = np.abs(ref_pool - ref_rank[j])
        out[g] = pool[np.lexsort((pool_names, d))[:n_ctrl]]
    return out


def main(conf_path, tables, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    conf = B.read_conf(conf_path)
    manifest = inventory(conf)
    pd.DataFrame(manifest).to_csv(out / "input_manifest.tsv", sep="\t", index=False)

    layers, organ_of, id_of, gene_index, universe = B.build_layers(conf, "universe")
    members = read_table(conf["MEMBERSHIP_TSV"], sep="\t")
    panel = read_table(conf["PANEL_TSV"], sep="\t")
    exclude = set(panel["gene"]) | set(members["gene"])
    genes, released = released_pool(tables)
    missing = [g for g in genes if g not in gene_index]
    assert not missing, missing

    # The analysed genes are scored through group_table, which adds the analysed list to the
    # exclusion set; the claim that their control sets are untouched is therefore checked here
    # directly rather than argued from the code.
    analysed = sorted({g for g in members["gene"].astype(str) if g in gene_index})

    audit, assign_b, assign_c, score_rows, target_rows = [], [], [], [], []
    for organ in ORGANS:
        ref_rows = np.where(organ_of["2D"] == organ)[0]
        ref_rank = np.median(layers["2D"][ref_rows, :], axis=0)
        base = nearest_without_self(ref_rank, genes, gene_index, exclude, 400)
        corr = L.control_sets(ref_rank, genes, gene_index, exclude, 400)
        target_base = nearest_without_self(ref_rank, analysed, gene_index, exclude, 400)
        target_corr = L.control_sets(ref_rank, analysed, gene_index, exclude, 400)
        for g in analysed:
            b_set, c_set = target_base[g].tolist(), target_corr[g].tolist()
            target_rows.append({
                "gene": g, "organ": organ,
                "identical_control_set": bool(b_set == c_set),
                "self_included_baseline": int(gene_index[g] in set(b_set)),
                "self_included_corrected": int(gene_index[g] in set(c_set)),
                "controls": len(c_set),
                "control_set_sha256_prefix_baseline": hashlib.sha256(
                    " ".join(map(str, b_set)).encode()).hexdigest()[:16],
                "control_set_sha256_prefix_corrected": hashlib.sha256(
                    " ".join(map(str, c_set)).encode()).hexdigest()[:16]})
        inv = {i: g for g, i in gene_index.items()}
        for g in genes:
            j = gene_index[g]
            b, c = base[g], corr[g]
            removed = sorted(set(b.tolist()) - set(c.tolist()))
            added = sorted(set(c.tolist()) - set(b.tolist()))
            audit.append({
                "gene": g, "organ": organ,
                "self_included_baseline": int(j in set(b.tolist())),
                "self_included_corrected": int(j in set(c.tolist())),
                "controls_baseline": int(len(b)), "controls_corrected": int(len(c)),
                "removed_gene": ", ".join(inv[i] for i in removed),
                "added_gene": ", ".join(inv[i] for i in added),
                "retained_order_preserved": bool([i for i in b.tolist() if i != j] == [i for i in c.tolist() if i in set(b.tolist())]),
                "reference_rank_target": float(ref_rank[j]),
                "max_control_gap_baseline": float(np.abs(ref_rank[b] - ref_rank[j]).max()),
                "max_control_gap_corrected": float(np.abs(ref_rank[c] - ref_rank[j]).max()),
            })
            for tag, s, sink in (("baseline", b, assign_b), ("corrected", c, assign_c)):
                sink.append({"gene": g, "organ": organ, "selection": tag,
                             "controls": " ".join(inv[i] for i in s),
                             "rank_gap_min": float(np.abs(ref_rank[s] - ref_rank[j]).min()),
                             "rank_gap_max": float(np.abs(ref_rank[s] - ref_rank[j]).max())})
            for lay, mat in layers.items():
                rows = np.where(organ_of[lay] == organ)[0]
                if not len(rows):
                    continue
                v_b = L.sample_scores(mat, rows, j, base[g])
                v_c = L.sample_scores(mat, rows, j, corr[g])
                for n, s_i in enumerate(rows):
                    score_rows.append({"layer": lay, "organ": organ, "gene": g,
                                       "sample_id": int(s_i), "cluster_id": id_of[lay][s_i],
                                       "score_baseline": float(v_b[n]), "score_corrected": float(v_c[n]),
                                       "difference": float(v_c[n] - v_b[n])})
        print("audited organ:", organ, flush=True)

    A = pd.DataFrame(audit)
    A.to_csv(out / "self_exclusion_audit.tsv", sep="\t", index=False)
    T = pd.DataFrame(target_rows)
    T.to_csv(out / "analysed_gene_control_sets.tsv", sep="\t", index=False)
    for sink, name in ((assign_b, "control_assignments_baseline.tsv.gz"),
                       (assign_c, "control_assignments_corrected.tsv.gz")):
        pd.DataFrame(sink).to_csv(out / name, sep="\t", index=False, compression="gzip")
    S = pd.DataFrame(score_rows)
    S.to_csv(out / "background_score_changes.tsv.gz", sep="\t", index=False, compression="gzip")

    # The reconciliation keys on the profile index, which is unique within a layer; the biological
    # cluster identifier is not (11 donors hold 20 healthy profiles, 202 patients hold 231), so it
    # is checked separately instead of being used as the key.
    key = ["layer", "organ", "gene", "sample_id"]
    left = S.set_index(key)["score_baseline"].sort_index()
    right = released.set_index(key)["score"].sort_index()
    assert left.index.is_unique and right.index.is_unique, "profile keys are not unique"
    assert left.index.equals(right.index), "released and rebuilt rows do not correspond"
    dev = float((left - right).abs().max())
    clusters = S.set_index(key)["cluster_id"].sort_index()
    cl_released = released.set_index(key)["cluster_id"].sort_index()
    report = {
        "universe": len(universe),
        "profiles": {lay: int(layers[lay].shape[0]) for lay in LAYERS},
        "unique_cluster_ids": {lay: int(pd.Series(id_of[lay]).nunique()) for lay in LAYERS},
        "rows_compared": int(len(left)),
        "max_abs_difference_to_released": dev,
        "cluster_ids_match": bool((clusters == cl_released).all()),
        "tolerance": 1e-9,
        "baseline_reproduced": bool(dev < 1e-9),
        "self_inclusion_baseline": int(A.self_included_baseline.sum()),
        "self_inclusion_corrected": int(A.self_included_corrected.sum()),
        "selections": int(len(A)),
        "control_count_always_400": bool((A.controls_corrected == 400).all()),
        "analysed_genes_checked": int(T.gene.nunique()),
        "analysed_gene_selections": int(len(T)),
        "analysed_control_sets_identical": bool(T.identical_control_set.all()),
        "analysed_self_inclusion_baseline": int(T.self_included_baseline.sum()),
    }
    json.dump(report, open(out / "baseline_reconstruction_report.json", "w"), indent=1)
    print(json.dumps(report, indent=1), flush=True)
    if not report["baseline_reproduced"]:
        print("FAILED: the released background table was not reproduced", flush=True)
        return 1
    if not report["cluster_ids_match"]:
        print("FAILED: cluster identifiers do not match the released table", flush=True)
        return 1
    if not report["analysed_control_sets_identical"]:
        print("FAILED: an analysed gene's control set changed", flush=True)
        return 1
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("configuration", type=Path)
    parser.add_argument("--tables", type=Path, required=True,
                        help="directory holding the RELEASED background table")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    sys.exit(main(args.configuration, args.tables, args.output_dir))
