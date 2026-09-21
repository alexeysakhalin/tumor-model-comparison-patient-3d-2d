"""Calculate expression scores from separately obtained primary expression inputs.

This stage requires prepared epithelial profiles, model expression and metadata.
It produces profile-level and gene-level scores. build_summary_tables.py then
calculates group estimates and intervals. See docs/PRIMARY_INPUTS.md for scope.
"""

import argparse
import os
import sys
import zipfile
from pathlib import Path
import numpy as np
import pandas as pd
from table_io import read_table

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import level_matched_scores as L

ORGANS = ["stomach", "large_intestine", "lung", "breast", "pancreas", "ovary"]
TABULA = {
    "stomach": "Stomach_ep",
    "large_intestine": "Large_Intestine_ep",
    "lung": "Lung_ep",
    "breast": "Mammary_ep",
    "pancreas": "Pancreas_ep",
    "ovary": "Ovary_ep",
}
LINEAGE = {
    "stomach": "Esophagus/Stomach",
    "large_intestine": "Bowel",
    "lung": "Lung",
    "breast": "Breast",
    "pancreas": "Pancreas",
    "ovary": "Ovary/Fallopian Tube",
}
LAYERS = ["healthy", "malignant", "dome", "2D"]


def read_conf(path):
    conf = {}
    for line in open(path, encoding="utf-8"):
        line = line.split("#", 1)[0].strip()
        if "=" in line:
            k, v = line.split("=", 1)
            conf[k.strip()] = v.strip().strip('"')
    return conf


def healthy_layer(pb_dir):
    """Tabula Sapiens epithelium: one column per donor, summed counts per gene.

    A manifest supplies each organ and its pseudobulk output prefix. Group
    metadata must be present and aligned with the matrix columns.
    """
    mats, organs, syms, ids = ([], [], None, [])
    man_path = f"{pb_dir}/manifest.csv"
    man = read_table(man_path)
    if not {"organ", "prefix"} <= set(man.columns) or set(man["organ"]) != set(ORGANS):
        raise ValueError(
            "Healthy manifest must identify all six organs and output prefixes"
        )
    pairs = list(zip(man["organ"], man["prefix"]))
    for organ, stem in pairs:
        stem = os.path.basename(str(stem))
        M = np.load(f"{pb_dir}/{stem}_counts.npy")
        g = read_table(f"{pb_dir}/{stem}_genes.csv")
        col = next((c for c in ("symbol", "gene") if c in g.columns), None)
        if col is None or g[col].isna().any():
            raise ValueError(f"{stem}: gene symbols are required")
        s = g[col].astype(str).to_numpy()
        if syms is None:
            syms = s
        elif not np.array_equal(syms, s):
            raise SystemExit(f"{stem}: gene axis differs from the first organ")
        mats.append(M)
        organs += [organ] * M.shape[1]
        gp = f"{pb_dir}/{stem}_groups.tsv"
        groups = read_table(gp, sep="\t")
        if M.shape != (len(g), len(groups)) or groups["group"].duplicated().any():
            raise ValueError(f"{stem}: gene or group axis is inconsistent")
        if not np.isfinite(M).all() or np.any(M < 0):
            raise ValueError(f"{stem}: invalid count matrix")
        ids += list(groups["group"].astype(str))
    return (np.hstack(mats).T, np.array(organs), syms, np.array(ids, dtype=str))


ORGAN_OF_ORIGIN = {
    "Stomach": "stomach",
    "Colon": "large_intestine",
    "Colorectum": "large_intestine",
    "Rectum": "large_intestine",
    "Lung": "lung",
    "Breast": "breast",
    "Pancreas": "pancreas",
    "Ovary": "ovary",
}
MIN_CELLS = 50


def patient_layer(pb_dir, studies_csv):
    """Pan-cancer atlas: malignant epithelial pseudobulk samples of the six organs.

    A sample enters the layer when its malignant-epithelium group holds at least MIN_CELLS cells;
    below that a pseudobulk profile is dominated by the few cells that were captured. This rule
    yields the cohort stated in PROV_kang_primary_studies.csv.
    """
    st = read_table(studies_csv)
    organ_of_ds = dict(zip(st["Dataset"], st["organ"]))
    mats, organs, syms, ids = ([], [], None, [])
    for ds, organ in organ_of_ds.items():
        f = f"{pb_dir}/{ds}_mean.npy"
        M = np.load(f)
        g = read_table(f"{pb_dir}/{ds}_genes.csv")
        col = next((c for c in ("symbol", "gene", "gene_id") if c in g.columns))
        s = g[col].astype(str).to_numpy()
        if M.shape[1] == len(s):
            M = M.T
        elif M.shape[0] != len(s):
            raise SystemExit(
                f"{ds}: matrix {M.shape} matches neither orientation of {len(s)} genes"
            )
        if syms is None:
            syms, ref = (s, ds)
        elif not np.array_equal(syms, s):
            raise SystemExit(f"{ds}: gene axis differs from {ref}")
        grp = read_table(f"{pb_dir}/{ds}_groups.tsv", sep="\t")
        if M.shape != (len(s), len(grp)) or grp["group"].duplicated().any():
            raise ValueError(f"{ds}: gene or group axis is inconsistent")
        if not np.isfinite(M).all():
            raise ValueError(f"{ds}: non-finite expression values")
        take = (
            grp["cell_class"].str.startswith("Epithelial_", na=False)
            & grp["cell_class"].str.contains("tumor|malignant", case=False, na=False)
            & (grp["n_cells"] >= MIN_CELLS)
            & grp["Organ_origin"].map(ORGAN_OF_ORIGIN).notna()
        ).to_numpy()
        if not take.any():
            continue
        mats.append(M[:, take])
        organs += list(grp.loc[take, "Organ_origin"].map(ORGAN_OF_ORIGIN))
        if "Patient" not in grp or grp.loc[take, "Patient"].isna().any():
            raise ValueError(f"{ds}: patient identifiers are required")
        pid = grp.loc[take, "Patient"].astype(str)
        ids += [f"{ds}|{p}" for p in pid]
    return (np.hstack(mats).T, np.array(organs), syms, np.array(ids, dtype=str))


def depmap_2d(depmap_dir):
    """DepMap traditional models: log2(TPM+1) per cell line, lineage from Model.csv."""
    expr = read_table(
        f"{depmap_dir}/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv", index_col=0
    )
    if not {"ModelID", "is_default_entry"} <= set(expr.columns):
        raise ValueError(
            "The recorded expression schema requires ModelID and is_default_entry"
        )
    meta_cols = [
        c for c in ("ProfileID", "is_default_entry", "ModelID") if c in expr.columns
    ]
    info = expr[meta_cols]
    expr = expr.drop(columns=meta_cols)
    expr.columns = [c.split(" (")[0] for c in expr.columns]
    mod_meta = read_table(f"{depmap_dir}/Model.csv").set_index("ModelID")
    if mod_meta.index.duplicated().any():
        raise ValueError("Duplicate model identifiers in Model.csv")
    mod = mod_meta["OncotreeLineage"]
    mod_growth = mod_meta["GrowthPattern"]
    lin = info["ModelID"].map(mod)
    keep = lin.isin(LINEAGE.values()).to_numpy()
    if "is_default_entry" in info.columns:
        keep &= (
            info["is_default_entry"]
            .astype(str)
            .str.lower()
            .isin(["true", "1"])
            .to_numpy()
        )
    keep &= (info["ModelID"].map(mod_growth) == "Adherent").to_numpy()
    org_all = (
        pd.Series(lin.to_numpy()).map({v: k for k, v in LINEAGE.items()}).to_numpy()
    )
    code_all = info["ModelID"].map(mod_meta["OncotreeCode"]).to_numpy()
    dis_all = info["ModelID"].map(mod_meta["OncotreePrimaryDisease"]).to_numpy()
    keep &= site_mask(code_all, dis_all, org_all)
    expr = expr[keep]
    if info.loc[keep, "ModelID"].duplicated().any():
        raise ValueError("More than one default expression profile per 2D model")
    organs = (
        pd.Series(lin[keep].to_numpy())
        .map({v: k for k, v in LINEAGE.items()})
        .to_numpy()
    )
    return (
        expr.to_numpy(dtype="float32"),
        organs,
        np.array(expr.columns),
        info.loc[keep, "ModelID"].astype(str).to_numpy(),
    )


SITE_CODES = {
    "stomach": {"STAD", "TSTAD", "DSTAD", "STSC", "MSTAD", "PSTAD", "SPDAC", "GRC"},
    "large_intestine": {"COAD", "READ", "COADREAD", "MACR", "CMC", "SRCCR"},
    "lung": {
        "LUAD",
        "LUSC",
        "SCLC",
        "LCLC",
        "NSCLC",
        "LUAS",
        "LUNE",
        "LUPC",
        "NSCLCPD",
        "SARCL",
    },
    "breast": {
        "IDC",
        "BRCA",
        "ILC",
        "DCIS",
        "BRCNOS",
        "BRCANOS",
        "MDLC",
        "IMMC",
        "BRCANOS",
    },
    "pancreas": {"PAAD", "PAASC", "PANET", "PAAC", "UCP"},
    "ovary": {"HGSOC", "SOC", "CCOV", "MOV", "EOV", "LGSOC", "OVT", "OCS", "SCCO"},
}


def site_mask(codes, diseases, organs):
    """True where the model's Oncotree code belongs to the organ it was assigned to."""
    ok = np.array([str(c) in SITE_CODES.get(o, set()) for c, o in zip(codes, organs)])
    return (
        ok
        & ~pd.Series(diseases)
        .astype(str)
        .str.contains("Non-Cancerous", na=False)
        .to_numpy()
    )


def depmap_3d(depmap_dir):
    """DepMap NextGen organoid models grown as matrigel domes."""
    with zipfile.ZipFile(f"{depmap_dir}/DepMap_NextGen_3D_core_files.zip") as z:
        with z.open("next_gen_expression.csv") as fh:
            expr = read_table(fh, index_col=0)
        with z.open("model_metadata.csv") as fh:
            meta = read_table(fh)
    expr.columns = [c.split(" (")[0] for c in expr.columns]
    col = next(
        (c for c in ("OncotreeLineage", "lineage", "Lineage") if c in meta.columns)
    )
    idc = next((c for c in ("ModelID", "model_id", "ModelId") if c in meta.columns))
    m = meta.set_index(idc)
    if m.index.duplicated().any() or expr.index.duplicated().any():
        raise ValueError("Duplicate model identifiers in the 3D input")
    lin = m[col].reindex(expr.index)
    growth = m["GrowthPattern"].reindex(expr.index)
    keep = (lin.isin(LINEAGE.values()) & (growth == "Dome")).to_numpy()
    org_all = (
        pd.Series(lin.to_numpy()).map({v: k for k, v in LINEAGE.items()}).to_numpy()
    )
    code_all = m["OncotreeCode"].reindex(expr.index).to_numpy()
    dis_all = m["OncotreePrimaryDisease"].reindex(expr.index).to_numpy()
    keep &= site_mask(code_all, dis_all, org_all)
    expr = expr[keep]
    organs = (
        pd.Series(lin[keep].to_numpy())
        .map({v: k for k, v in LINEAGE.items()})
        .to_numpy()
    )
    return (
        expr.to_numpy(dtype="float32"),
        organs,
        np.array(expr.columns),
        np.array(expr.index, dtype=str),
    )


def build_layers(conf, rank_scope="universe"):
    raw = {
        "healthy": healthy_layer(f"{conf['PSEUDOBULK_DIR']}/tabula"),
        "malignant": patient_layer(
            f"{conf['PSEUDOBULK_DIR']}/kang", conf["KANG_STUDIES_CSV"]
        ),
        "dome": depmap_3d(conf["DEPMAP_DIR"]),
        "2D": depmap_2d(conf["DEPMAP_DIR"]),
    }
    universe = None
    for lay, (_, _, syms, _) in raw.items():
        s = pd.Index(syms)
        s = s[~s.duplicated()]
        universe = s if universe is None else universe.intersection(s)
    universe = pd.Index(sorted(universe))
    if len(universe) != 18612:
        raise ValueError(
            f"The Figure 6 input universe has 18,612 genes, obtained {len(universe)}"
        )
    print("gene universe:", len(universe), flush=True)
    layers, organ_of, id_of = ({}, {}, {})
    for lay, (M, organs, syms, sample_ids) in raw.items():
        if (
            M.shape[0] != len(organs)
            or M.shape[1] != len(syms)
            or len(sample_ids) != len(organs)
        ):
            raise ValueError(f"{lay}: inconsistent matrix and annotation axes")
        if set(organs) != set(ORGANS) or not np.isfinite(M).all():
            raise ValueError(f"{lay}: incomplete organ coverage or non-finite values")
        expected = {"healthy": 20, "malignant": 231, "dome": 152, "2D": 401}
        if M.shape[0] != expected[lay]:
            raise ValueError(
                f"{lay}: {M.shape[0]} profiles differ from the Figure 6 cohort"
            )
        first = ~pd.Index(syms).duplicated()
        M, syms = (M[:, first], np.asarray(syms)[first])
        take = pd.Index(syms).get_indexer(universe)
        if (take < 0).any():
            raise SystemExit(f"{lay}: universe gene missing after intersection")
        if rank_scope == "layer":
            layers[lay] = L.percentile_ranks(M.astype("float64"))[:, take]
        else:
            layers[lay] = L.percentile_ranks(M[:, take].astype("float64"))
        organ_of[lay] = organs
        id_of[lay] = np.asarray(sample_ids, dtype=str)
        print(
            f"{lay}: {layers[lay].shape[0]} samples, rank scope {rank_scope} ({(M.shape[1] if rank_scope == 'layer' else len(universe))} genes)",
            flush=True,
        )
    return (layers, organ_of, id_of, {g: i for i, g in enumerate(universe)}, universe)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("configuration", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--rank-scope", choices=["universe", "layer"], default="universe"
    )
    args = parser.parse_args()
    conf = read_conf(args.configuration)
    for key in [
        "PSEUDOBULK_DIR",
        "KANG_STUDIES_CSV",
        "DEPMAP_DIR",
        "MEMBERSHIP_TSV",
        "PANEL_TSV",
    ]:
        if key not in conf or not Path(conf[key]).exists():
            parser.error(f"Missing input for {key}")
    out = args.output_dir
    if out.resolve() == Path(conf["MEMBERSHIP_TSV"]).resolve().parent:
        parser.error("Output must not overwrite the deposited tables")
    out.mkdir(parents=True, exist_ok=True)
    members = read_table(conf["MEMBERSHIP_TSV"], sep="\t")
    panel = read_table(conf["PANEL_TSV"], sep="\t")
    rank_scope = args.rank_scope
    layers, organ_of, id_of, gene_index, universe = build_layers(conf, rank_scope)
    genes = [g for g in members["gene"] if g in gene_index]
    if len(genes) != 67 or len(set(genes)) != 67 or members["group"].nunique() != 12:
        raise ValueError("All 67 unique genes in 12 groups must be available")
    exclude = set(panel["gene"]) | set(members["gene"])
    per_gene = L.group_table(
        genes, ORGANS, layers, "2D", gene_index, organ_of, exclude=exclude
    )
    per_gene = per_gene.merge(
        members[["gene", "group"]].rename(columns={"group": "group"}), on="gene"
    )
    per_gene.rename(columns={"gene": "gene", "organ": "organ"}).to_csv(
        f"{out}/IKM_final_per_gene_by_organ.csv", index=False
    )
    print("per-gene table:", per_gene.shape, flush=True)
    rng = np.random.default_rng(11)
    pool = [g for g in universe if g not in exclude]
    rand = list(rng.choice(pool, 200, replace=False))
    grp_of = dict(zip(members["gene"], members["group"]))
    ps_rows, rand_rows = ([], [])
    for organ in ORGANS:
        ref_rows = np.where(organ_of["2D"] == organ)[0]
        if len(ref_rows) < 10:
            raise ValueError(f"{organ}: at least 10 reference profiles are required")
        ref_rank = np.median(layers["2D"][ref_rows, :], axis=0)
        csets = L.control_sets(ref_rank, genes, gene_index, exclude, 400)
        cs_rand = L.control_sets(ref_rank, rand, gene_index, exclude, 400)
        for lay, mat in layers.items():
            rows = np.where(organ_of[lay] == organ)[0]
            if not len(rows):
                continue
            for grp in sorted(set(grp_of.values())):
                gl = [
                    g
                    for g in members.loc[members["group"] == grp, "gene"]
                    if g in csets
                ]
                if not gl:
                    continue
                med = np.median(
                    np.vstack(
                        [
                            L.sample_scores(mat, rows, gene_index[g], csets[g])
                            for g in gl
                        ]
                    ),
                    axis=0,
                )
                for i, s_i in enumerate(rows):
                    ps_rows.append(
                        (lay, organ, int(s_i), id_of[lay][s_i], grp, float(med[i]))
                    )
            for g in rand:
                if g not in cs_rand:
                    continue
                v = L.sample_scores(mat, rows, gene_index[g], cs_rand[g])
                for i, s_i in enumerate(rows):
                    rand_rows.append(
                        (lay, organ, g, int(s_i), id_of[lay][s_i], float(v[i]))
                    )
        print("scored organ:", organ, flush=True)
    pd.DataFrame(
        ps_rows, columns=["layer", "organ", "sample_id", "cluster_id", "group", "score"]
    ).to_csv(f"{out}/IKM_group_scores_by_sample.csv", index=False)
    pd.DataFrame(
        rand_rows,
        columns=["layer", "organ", "gene", "sample_id", "cluster_id", "score"],
    ).to_csv(
        f"{out}/IKM_random_gene_scores_by_sample.csv.gz",
        index=False,
        compression="gzip",
    )
    prov = []
    for lay in LAYERS:
        for organ, ident in zip(organ_of[lay], id_of[lay]):
            prov.append({"layer": lay, "organ": organ, "cluster_id": ident})
    pv = pd.DataFrame(prov)
    pv.to_csv(f"{out}/PROV_sample_identities.csv", index=False)
    summ = (
        pv.groupby("layer")
        .agg(profile_count=("cluster_id", "size"), unique_ids=("cluster_id", "nunique"))
        .reset_index()
    )
    summ.to_csv(f"{out}/PROV_model_layers.csv", index=False)
    print(summ.to_string(index=False), flush=True)
    print(
        "per-sample rows:",
        len(ps_rows),
        "| random-gene rows:",
        len(rand_rows),
        flush=True,
    )


if __name__ == "__main__":
    main()
