"""Reconstruct the MLKL/PGAM5 sensitivity analysis from deposited tables.

The original rounded sensitivity table is retained. This script also enumerates
all 200 reference genes to assess Monte Carlo variation in one-gene bands.
Usage: python code/recompute_pgam5_sensitivity.py tables --check
       python code/recompute_pgam5_sensitivity.py tables --output-dir results
"""

from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
from table_io import read_table

LAYERS = ["healthy", "malignant", "dome", "2D"]
NAMES = dict(
    zip(
        LAYERS,
        [
            "Healthy epithelium",
            "Malignant epithelium",
            "3D dome models",
            "2D adherent models",
        ],
    )
)


def reference_bands(tables):
    """One-gene reference bands from the deposited background scores.

    Returned as {layer: ((mc_low, mc_high), (exact_low, exact_high))} together with the per-gene
    reference values. Monte Carlo: 2,000 indices drawn with replacement, NumPy default_rng(3),
    sorted gene symbols; exact: each of the 200 genes once. Factored out of `compute` so that a
    correction which changes the background scores regenerates the deposited bands with this same
    implementation instead of a second copy of it.
    """
    tables = Path(tables)
    background = read_table(tables / "IKM_random_gene_scores_by_sample.csv.gz")
    per_site = background.groupby(["layer", "gene", "organ"])["score"].median()
    reference = per_site.groupby(["layer", "gene"]).median()
    genes = sorted(background["gene"].unique())
    assert len(genes) == 200
    assert per_site.groupby(["layer", "gene"]).size().eq(6).all()
    draws = np.random.default_rng(3).choice(len(genes), size=2000, replace=True)
    distributions = []
    bands = {}
    for layer in LAYERS:
        values = reference.loc[layer].reindex(genes).to_numpy()
        assert np.isfinite(values).all()
        bands[layer] = (
            np.quantile(values[draws], [0.025, 0.975], method="linear"),
            np.quantile(values, [0.025, 0.975], method="linear"),
        )
        distributions.extend(
            (
                {"layer": NAMES[layer], "gene": g, "score_pp": float(v)}
                for g, v in zip(genes, values)
            )
        )
    return bands, distributions, genes


def observed_scores(tables):
    """Full-precision statistics of the analysed genes: the singles and the plotted pair.

    Single genes: the median across the six sites of the site-specific gene scores, which for one
    gene is the same number under either aggregation order. The pair: the sample-first estimator of
    the figure, taken from the group scores, never the mean of the two final single-gene values,
    because medians are not associative. Factored out so that a correction regenerating the
    deposited table uses these numbers rather than the one-decimal display values.
    """
    tables = Path(tables)
    target = read_table(tables / "IKM_final_per_gene_by_organ.csv")
    observed = (
        target[target["gene"].isin(["MLKL", "PGAM5"])].groupby("gene")[LAYERS].median()
    )
    sample = read_table(tables / "IKM_group_scores_by_sample.csv")
    pair = (
        sample[sample["group"] == "Necroptosis executioners"]
        .groupby(["layer", "organ"])["score"]
        .median()
        .groupby("layer")
        .median()
    )
    return observed, pair


def compute(tables):
    tables = Path(tables)
    bands, distributions, genes = reference_bands(tables)
    observed, pair = observed_scores(tables)
    original = read_table(tables / "IKM_pgam5_sensitivity.csv")
    two_gene = read_table(tables / "IKM_null_bands.csv")
    rows = []
    checks = []
    for composition, label, size in [
        ("MLKL + PGAM5", "MLKL + PGAM5 (figure)", 2),
        ("MLKL", "MLKL only", 1),
        ("PGAM5", "PGAM5 only", 1),
    ]:
        for layer in LAYERS:
            expected = original[
                (original["composition"] == label) & (original["layer"] == layer)
            ].iloc[0]
            score = float(
                pair.loc[layer] if size == 2 else observed.loc[composition, layer]
            )
            if size == 1:
                mc, exact = bands[layer]
                exact_low, exact_high = map(float, exact)
                exact_out = bool(score < exact_low or score > exact_high)
            else:
                b = two_gene[
                    (two_gene["group_size"] == 2) & (two_gene["layer"] == layer)
                ].iloc[0]
                mc = [float(b["band_lower"]), float(b["band_upper"])]
                exact_low = exact_high = np.nan
                exact_out = np.nan
            low, high = map(float, mc)
            outside = bool(score < low or score > high)
            checks.append(
                {
                    "composition": composition,
                    "layer": NAMES[layer],
                    "score_matches_1dp": abs(
                        round(score, 1) - float(expected["score_pp"])
                    )
                    < 1e-08,
                    "band_matches_1dp": abs(
                        round(low, 1) - float(expected["band_lower"])
                    )
                    < 1e-08
                    and abs(round(high, 1) - float(expected["band_upper"])) < 1e-08,
                    "outside_flag_matches": outside == bool(expected["outside_band"]),
                }
            )
            rows.append(
                {
                    "composition": composition,
                    "gene_count": size,
                    "layer": NAMES[layer],
                    "score_pp": score,
                    "reported_MC_band_lower_pp": low,
                    "reported_MC_band_upper_pp": high,
                    "outside_reported_MC_band": outside,
                    "exact_200_gene_band_lower_pp": exact_low,
                    "exact_200_gene_band_upper_pp": exact_high,
                    "outside_exact_200_gene_band": exact_out,
                    "reference_method": "Size-two band retained from Figure 6"
                    if size == 2
                    else "MC: 2000 draws, NumPy default_rng(3), sorted gene symbols; exact: each of 200 genes once",
                }
            )
    raw = read_table(tables / "IKM_pgam5_screen_gene_values.csv")
    supplied = read_table(tables / "IKM_pgam5_sensitivity_arms.csv")
    screen_checks = []
    for row in supplied.to_dict("records"):
        chosen = ["MLKL", "PGAM5"] if row["composition"] == "MLKL + PGAM5" else ["MLKL"]
        src = raw[(raw["model"] == row["arm"]) & raw["gene"].isin(chosen)]
        assert len(src) == len(chosen)
        for condition in ["baseline", "IFNg", "IL6"]:
            calculated = float(src[condition].median())
            error = abs(calculated - float(row[condition]))
            screen_checks.append(
                {
                    "arm": row["arm"],
                    "composition": row["composition"],
                    "condition": condition,
                    "recomputed_from_rounded_genes": calculated,
                    "supplied": float(row[condition]),
                    "absolute_difference": error,
                    "within_input_rounding": error <= 0.000500001,
                }
            )
    assert all(
        (
            all(
                (
                    c[k]
                    for k in [
                        "score_matches_1dp",
                        "band_matches_1dp",
                        "outside_flag_matches",
                    ]
                )
            )
            for c in checks
        )
    ), checks
    assert all((c["within_input_rounding"] for c in screen_checks))
    result = pd.DataFrame(rows)
    reference_table = pd.DataFrame(distributions)
    return (
        result,
        reference_table,
        {
            "expression_checks": checks,
            "screen_checks": screen_checks,
            "scope": "Derived-table reconstruction. Individual-gene cluster-bootstrap intervals and raw-data reanalysis are not performed.",
        },
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tables", type=Path)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    result, reference, checks = compute(args.tables)
    if args.check:
        for name, frame in [
            ("IKM_pgam5_sensitivity_recomputed.csv", result),
            ("IKM_pgam5_reference_gene_scores.csv", reference),
        ]:
            expected = read_table(args.tables / name)
            pd.testing.assert_frame_equal(
                frame,
                expected,
                check_dtype=False,
                check_exact=False,
                rtol=1e-08,
                atol=1e-10,
            )
        print(
            "PASS 12 expression estimates and reported bands; 24 screen values within source rounding; exact one-gene reference distribution"
        )
    if args.output_dir:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        result.to_csv(
            args.output_dir / "IKM_pgam5_sensitivity_recomputed.csv",
            index=False,
            float_format="%.12g",
        )
        reference.to_csv(
            args.output_dir / "IKM_pgam5_reference_gene_scores.csv",
            index=False,
            float_format="%.12g",
        )
        (args.output_dir / "PGAM5_RECONSTRUCTION_CHECKS.json").write_text(
            json.dumps(checks, indent=2) + "\n"
        )


if __name__ == "__main__":
    main()
