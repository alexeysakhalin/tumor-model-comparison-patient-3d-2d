"""Recalculate Figure 6 expression estimates from deposited per-profile scores.

Scores are aggregated over genes within a profile, profiles within an organ,
and the six organs. Intervals use 2,000 accepted cluster-bootstrap draws, with
donors, studies or models resampled jointly across organs within each layer.
Draws missing an organ are rejected. Reference bands use size-matched random
gene sets and the same estimator. The auxiliary gene-first table is separate.
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from table_io import read_table

LAYERS = ["healthy", "malignant", "dome", "2D"]
ORGANS = ["stomach", "large_intestine", "lung", "breast", "pancreas", "ovary"]
RECOGNITION = ("MHC-I presentation", "NK activating ligands", "Inhibitory signals")


def by_organ_layer(ps, grp):
    """{(organ, layer): (per-sample values, cluster label per value)} for one group.

    Values and cluster labels come from the same rows in the same order, so a bootstrap over
    clusters cannot line up a label with another sample's value.
    """
    sub = ps[ps["group"] == grp]
    out = {}
    for o in ORGANS:
        for l in LAYERS:
            d = sub[(sub["organ"] == o) & (sub["layer"] == l)]
            ident = d["cluster_id"].astype(str)
            cl = (ident.str.split("|").str[0] if l == "malignant" else ident).to_numpy()
            out[o, l] = (d["score"].to_numpy(), cl)
    return out


def point_and_ci(by, lay=None, diff=None, n=2000, rng=None, reject_log=None):
    """Median over samples within an organ, then over the SIX organs; cluster bootstrap.

    Inferential target. The estimate summarises a fixed set of six organs, each represented by the
    samples available for it. The resampling unit is the cluster - the study for tumour samples, the
    donor for healthy tissue, the model for the two culture layers - because a donor contributes
    epithelium to several organs and a study contributes patients to several organs, so sample rows
    are not independent.

    Procedure. Each replicate draws the layer's clusters once and applies that draw to every organ,
    which keeps a donor either inside or outside the whole replicate. A draw that leaves an organ
    with no sample is not a smaller version of the estimand - it is a different estimand over five
    organs - so such draws are rejected and redrawn, and the rejection rate is recorded in
    `IKM_bootstrap_rejection.csv` rather than left implicit. A stratum served by a single cluster
    keeps that cluster in every replicate and therefore contributes no variability of its own.
    """
    rng = rng or np.random.default_rng(7)
    layers_used = list(diff) if diff else [lay]
    if n < 1 or any(l not in LAYERS for l in layers_used):
        raise ValueError(
            "A valid layer or contrast and positive draw count are required"
        )
    for layer in layers_used:
        for organ in ORGANS:
            values, clusters = by[organ, layer]
            if (
                len(values) == 0
                or len(values) != len(clusters)
                or not np.isfinite(values).all()
            ):
                raise ValueError(
                    f"Missing, misaligned or non-finite scores: {organ}, {layer}"
                )
    per_organ = []
    for o in ORGANS:
        if diff:
            a, b = diff
            va, vb = (by[o, a][0], by[o, b][0])
            if len(va) and len(vb):
                per_organ.append(np.median(va) - np.median(vb))
        elif len(by[o, lay][0]):
            per_organ.append(np.median(by[o, lay][0]))
    if not per_organ:
        return (np.nan, np.nan, np.nan, 0, 0)
    pt = float(np.median(per_organ))
    layers_used = [diff[0], diff[1]] if diff else [lay]
    present = {l: [o for o in ORGANS if len(by[o, l][0])] for l in layers_used}
    pools = {
        l: np.unique(np.concatenate([by[o, l][1] for o in present[l]]))
        if present[l]
        else np.array([])
        for l in layers_used
    }

    ordered = {}
    for layer in layers_used:
        for organ in ORGANS:
            values, labels = by[organ, layer]
            order = np.argsort(values, kind="stable")
            ordered[organ, layer] = (
                values[order],
                np.searchsorted(pools[layer], labels[order]),
            )

    def organ_median(o, l, multiplicities):
        values, cluster_positions = ordered[o, l]
        cumulative = np.cumsum(multiplicities[cluster_positions])
        count = int(cumulative[-1])
        if count == 0:
            return None
        # Integer weights reproduce explicit duplication of each sampled cluster.
        low = np.searchsorted(cumulative, (count - 1) // 2, side="right")
        high = np.searchsorted(cumulative, count // 2, side="right")
        return float((values[low] + values[high]) / 2)

    boots, rejected, attempts = ([], 0, 0)
    while len(boots) < n and attempts < n * 50:
        attempts += 1
        drawn = {
            l: np.bincount(rng.choice(len(p), len(p), replace=True), minlength=len(p))
            for l, p in pools.items()
        }
        vals, complete = ([], True)
        for o in ORGANS:
            if diff:
                a, b = diff
                ma, mb = (organ_median(o, a, drawn[a]), organ_median(o, b, drawn[b]))
                if o in present[a] and o in present[b] and (ma is None or mb is None):
                    complete = False
                    break
                if ma is not None and mb is not None:
                    vals.append(ma - mb)
            else:
                m = organ_median(o, lay, drawn[lay])
                if o in present[lay] and m is None:
                    complete = False
                    break
                if m is not None:
                    vals.append(m)
        if not complete:
            rejected += 1
            continue
        if vals:
            boots.append(float(np.median(vals)))
    if reject_log is not None:
        reject_log.append(
            {
                "quantity": f"{lay or diff}",
                "accepted_replicates": len(boots),
                "rejected_draws": rejected,
                "rejection_fraction": round(rejected / max(attempts, 1), 4),
            }
        )
    same = int(np.sum(np.sign(per_organ) == np.sign(pt)))
    if len(boots) != n:
        raise RuntimeError(f"Obtained {len(boots)} accepted replicates; {n} required")
    return (
        pt,
        float(np.percentile(boots, 2.5)),
        float(np.percentile(boots, 97.5)),
        len(per_organ),
        same,
    )


def single_cluster_strata(ps):
    """(layer, organ) pairs served by one cluster only: they carry no resampling variability."""
    rows = []
    for (l, o), sub in ps.groupby(["layer", "organ"]):
        ident = sub["cluster_id"].astype(str)
        cl = ident.str.split("|").str[0] if l == "malignant" else ident
        if cl.nunique() == 1:
            rows.append(
                {
                    "layer": l,
                    "organ": o,
                    "cluster_count": 1,
                    "sample_count": int(sub["sample_id"].nunique()),
                }
            )
    return pd.DataFrame(rows)


DIFF_KEYS = [("malignant", "healthy"), ("dome", "2D")]


def null_bands(rs, sizes, n_rep=2000, seed=3):
    """95% band of the same statistic for k random genes, by the sample route."""
    rng = np.random.default_rng(seed)
    genes = rs["gene"].unique()
    if rs.duplicated(["layer", "organ", "gene", "sample_id"]).any():
        raise ValueError("Duplicate background-gene profile keys")
    cube = {
        (o, l): rs[(rs["organ"] == o) & (rs["layer"] == l)]
        .pivot_table(index="gene", columns="sample_id", values="score")
        .reindex(genes)
        .to_numpy()
        for o in ORGANS
        for l in LAYERS
        if len(rs[(rs["organ"] == o) & (rs["layer"] == l)])
    }
    if (
        len(genes) != 200
        or len(cube) != 24
        or any(not np.isfinite(a).all() for a in cube.values())
    ):
        raise ValueError(
            "Reference data must cover 200 genes and all 24 organ-layer strata"
        )
    out = []
    for k in sizes:
        draws = np.array(
            [rng.choice(len(genes), k, replace=False) for _ in range(n_rep)]
        )
        med = {
            key: np.median(np.median(array[draws], axis=1), axis=1)
            for key, array in cube.items()
        }
        stats = {
            layer: np.median([med[organ, layer] for organ in ORGANS], axis=0)
            for layer in LAYERS
        }
        stats.update(
            {
                f"{a}-{b}": np.median([med[o, a] - med[o, b] for o in ORGANS], axis=0)
                for a, b in DIFF_KEYS
            }
        )
        for key, v in stats.items():
            out.append(
                {
                    "group_size": int(k),
                    "layer": key,
                    "band_lower": round(float(np.percentile(v, 2.5)), 1),
                    "band_upper": round(float(np.percentile(v, 97.5)), 1),
                }
            )
    return pd.DataFrame(out)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tables", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.tables.resolve() == args.output_dir.resolve():
        parser.error(
            "The output directory must differ from the deposited input directory"
        )
    src, out = args.tables, args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    ps = read_table(src / "IKM_group_scores_by_sample.csv")
    rs = read_table(src / "IKM_random_gene_scores_by_sample.csv.gz")
    members = read_table(src / "IKM_gene_membership_by_group.tsv", sep="\t")
    if ps.duplicated(["layer", "organ", "sample_id", "group"]).any():
        raise ValueError("Duplicate sample-group keys")
    sizes = sorted(members.groupby("group").size().unique())
    nb = null_bands(rs, sizes)
    nb.to_csv(f"{out}/IKM_null_bands.csv", index=False)
    band = {
        (int(r["group_size"]), r["layer"]): (r["band_lower"], r["band_upper"])
        for _, r in nb.iterrows()
    }
    rng = np.random.default_rng(7)
    reject_log = []
    mt, ci_rows, rep = ([], [], [])
    for grp, sub in members.groupby("group"):
        print(f"Calculating {grp}", flush=True)
        k = int(sub["gene"].nunique())
        by = by_organ_layer(ps, grp)
        rec, ci = ({"group": grp, "gene_count": k}, {"group": grp, "gene_count": k})
        for lay in LAYERS:
            v, lo, hi, _, _ = point_and_ci(by, lay=lay, rng=rng, reject_log=reject_log)
            b = band.get((k, lay), (np.nan, np.nan))
            outside = bool(not b[0] <= v <= b[1]) if np.isfinite(b[0]) else False
            rec[lay], rec[f"{lay}_outside_reference_band"] = (round(v, 1), outside)
            ci.update(
                {
                    lay: round(v, 1),
                    f"{lay}_lo": round(lo, 1),
                    f"{lay}_hi": round(hi, 1),
                    f"{lay}_outside_band": outside,
                }
            )
        for (a, b_), tag, col in [
            (
                ("malignant", "healthy"),
                "malignant_minus_healthy",
                "malignant_minus_healthy",
            ),
            (("dome", "2D"), "dome_minus_2D", "dome_minus_2D"),
        ]:
            v, lo, hi, n_o, same = point_and_ci(
                by, diff=(a, b_), rng=rng, reject_log=reject_log
            )
            bb = band.get((k, f"{a}-{b_}"), (np.nan, np.nan))
            outside = bool(not bb[0] <= v <= bb[1]) if np.isfinite(bb[0]) else False
            half = (bb[1] - bb[0]) / 2 if np.isfinite(bb[0]) else np.nan
            rec[col] = round(v, 1)
            ci.update(
                {
                    tag: round(v, 1),
                    f"{tag}_lo": round(lo, 1),
                    f"{tag}_hi": round(hi, 1),
                    f"{tag}_CI_excludes_zero": bool(lo > 0 or hi < 0),
                    f"{tag}_band_lo": bb[0],
                    f"{tag}_band_hi": bb[1],
                    f"{tag}_outside_band": outside,
                    f"{tag}_organ_count": f"{same}/{n_o}",
                }
            )
            rep.append(
                {
                    "group": grp,
                    "column": "patient − healthy"
                    if tag == "malignant_minus_healthy"
                    else "organoid − 2D",
                    "block": "representation",
                    "score": round(v, 1),
                    "reference_band_units": round(v / half, 2) if half else np.nan,
                    "outside_band": outside,
                    "CI_excludes_zero": bool(lo > 0 or hi < 0),
                    "organs_with_same_sign": f"{same}/{n_o}",
                    "all_sources_agree": same == n_o,
                }
            )
        mt.append(rec)
        ci_rows.append(ci)
    pg = read_table(src / "IKM_final_per_gene_by_organ.csv")
    gene_rows = []
    for grp, sub in pg.groupby("group"):
        rec = {"group": grp, "gene_count": int(sub["gene"].nunique())}
        per_organ = sub.groupby("organ")[LAYERS].median()
        for lay in LAYERS:
            rec[lay] = round(float(per_organ[lay].median()), 1)
        for (a, b_), col in [
            (("malignant", "healthy"), "malignant_minus_healthy"),
            (("dome", "2D"), "dome_minus_2D"),
        ]:
            d = sub.assign(_d=sub[a] - sub[b_]).groupby("organ")["_d"].median()
            rec[col] = round(float(d.median()), 1)
        rec["estimator"] = "gene route (methods 3a)"
        gene_rows.append(rec)
    pd.DataFrame(gene_rows).to_csv(f"{out}/IKM_final_master_table.csv", index=False)
    pd.DataFrame(ci_rows).to_csv(f"{out}/IKM_final_master_table_CI.csv", index=False)
    pd.DataFrame(rep).to_csv(f"{out}/IKM_representation_in_band_units.csv", index=False)
    mv = []
    for grp in sorted((g for g in members["group"].unique() if g not in RECOGNITION)):
        by = by_organ_layer(ps, grp)
        v = {lay: point_and_ci(by, lay=lay, n=1, rng=rng)[0] for lay in LAYERS}
        d_org, d_2d = (abs(v["dome"] - v["malignant"]), abs(v["2D"] - v["malignant"]))
        mv.append(
            {
                "group": grp,
                "gene_count": int(members[members["group"] == grp]["gene"].nunique()),
                "healthy": round(v["healthy"], 1),
                "malignant_tumor": round(v["malignant"], 1),
                "dome": round(v["dome"], 1),
                "2D": round(v["2D"], 1),
                "healthy_minus_malignant": round(v["healthy"] - v["malignant"], 1),
                "dome_minus_malignant": round(v["dome"] - v["malignant"], 1),
                "2D_minus_malignant": round(v["2D"] - v["malignant"], 1),
                "closer_to_malignant": "dome" if d_org < d_2d else "2D cell line",
                "distance_reduction_pp": round(d_2d - d_org, 1),
            }
        )
    pd.DataFrame(mv).to_csv(f"{out}/IKM_models_vs_patient.csv", index=False)
    pd.DataFrame(reject_log).to_csv(f"{out}/IKM_bootstrap_rejection.csv", index=False)
    sc = single_cluster_strata(ps)
    sc.to_csv(f"{out}/IKM_single_cluster_strata.csv", index=False)
    print("single-cluster strata:", len(sc), flush=True)
    print(
        "written:",
        len(mt),
        "groups |",
        len(mv),
        "death groups |",
        len(nb),
        "band rows",
        flush=True,
    )


if __name__ == "__main__":
    main()
