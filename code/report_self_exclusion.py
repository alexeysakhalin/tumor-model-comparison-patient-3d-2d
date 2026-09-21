"""What the background self-exclusion correction changed, table by table.

Compares the released tables with the corrected ones and writes the differences at full precision,
in the convention each table uses: reference bands are stored to one decimal place, and a normalised
contrast is the full-precision statistic divided by half the width of the stored band, rounded to two
decimals. The comparison therefore measures the correction against the released numbers and not
against a reconstruction in a different convention.

Usage: python code/report_self_exclusion.py --released DIR --corrected DIR --output-dir DIR
"""
import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_summary_tables import DIFF_KEYS, null_bands
from table_io import read_table

LAYERS = ["healthy", "malignant", "dome", "2D"]


def full_precision_bands(background, sizes):
    """Band limits before the release rounds them, from the same routine the pipeline uses.

    `null_bands(..., decimals=None)` is the released routine with its rounding switched off, so the
    comparison of band edges is not quantised to the stored 0.1 grid.
    """
    return null_bands(read_table(background), sizes, decimals=None)


def main(released, corrected, out):
    released, corrected, out = Path(released), Path(corrected), Path(out)
    out.mkdir(parents=True, exist_ok=True)
    summary = {}

    members = read_table(corrected / "IKM_gene_membership_by_group.tsv", sep="\t")
    sizes = sorted(members.groupby("group").size().unique())
    nb_r = full_precision_bands(released / "IKM_random_gene_scores_by_sample.csv.gz", sizes)
    nb_c = full_precision_bands(corrected / "IKM_random_gene_scores_by_sample.csv.gz", sizes)
    key = ["group_size", "layer"]
    bands = nb_r.merge(nb_c, on=key, suffixes=("_released", "_corrected"))
    bands["width_released"] = bands.band_upper_released - bands.band_lower_released
    bands["width_corrected"] = bands.band_upper_corrected - bands.band_lower_corrected
    bands["width_change"] = bands.width_corrected - bands.width_released
    bands["lower_change"] = bands.band_lower_corrected - bands.band_lower_released
    bands["upper_change"] = bands.band_upper_corrected - bands.band_upper_released
    bands.to_csv(out / "reference_band_changes.tsv", sep="\t", index=False)
    summary["band_rows"] = int(len(bands))
    summary["all_bands_wider"] = bool((bands.width_change > 0).all())
    summary["width_change_range"] = [float(bands.width_change.min()), float(bands.width_change.max())]
    summary["max_abs_edge_change"] = float(pd.concat([bands.lower_change.abs(), bands.upper_change.abs()]).max())
    summary["median_abs_edge_change"] = float(pd.concat([bands.lower_change.abs(), bands.upper_change.abs()]).median())

    rep_r = read_table(released / "IKM_representation_in_band_units.csv")
    rep_c = read_table(corrected / "IKM_representation_in_band_units.csv")
    rep = rep_r.merge(rep_c, on=["group", "column", "block"], suffixes=("_released", "_corrected"))
    rep["band_units_change"] = rep.reference_band_units_corrected - rep.reference_band_units_released
    rep["classification_changed"] = rep.outside_band_released != rep.outside_band_corrected
    rep.to_csv(out / "representation_changes.tsv", sep="\t", index=False)
    changed = rep[rep.band_units_change.abs() > 1e-12]
    summary["contrasts"] = int(len(rep))
    summary["contrasts_changed"] = int(len(changed))
    summary["max_abs_band_unit_change"] = float(rep.band_units_change.abs().max())
    summary["contrast_classification_changes"] = int(rep.classification_changed.sum())

    ci_r = read_table(released / "IKM_final_master_table_CI.csv")
    ci_c = read_table(corrected / "IKM_final_master_table_CI.csv")
    rows = []
    for lay in LAYERS + [f"{a}_minus_{b}" for a, b in DIFF_KEYS]:
        flag = f"{lay}_outside_band"
        if flag not in ci_r.columns:
            continue
        # Column names such as "2D" are not valid attribute names, so the merged frame is read by
        # label rather than through itertuples attributes.
        m = ci_r[["group", lay, flag]].merge(ci_c[["group", lay, flag]], on="group",
                                             suffixes=("_released", "_corrected"))
        for _, r in m.iterrows():
            rows.append({"group": r["group"], "column": lay,
                         "score_released": r[f"{lay}_released"],
                         "score_corrected": r[f"{lay}_corrected"],
                         "outside_released": r[f"{flag}_released"],
                         "outside_corrected": r[f"{flag}_corrected"]})
    cls = pd.DataFrame(rows)
    cls["classification_changed"] = cls.outside_released != cls.outside_corrected
    cls["score_changed"] = cls.score_released != cls.score_corrected
    cls.to_csv(out / "layer_classification_changes.tsv", sep="\t", index=False)
    summary["classification_rows"] = int(len(cls))
    summary["layer_classification_changes"] = int(cls.classification_changed.sum())
    summary["layer_score_changes"] = int(cls.score_changed.sum())

    s_r = read_table(released / "IKM_pgam5_sensitivity.csv")
    s_c = read_table(corrected / "IKM_pgam5_sensitivity.csv")
    sens = s_r.merge(s_c, on=["composition", "gene_count", "layer"], suffixes=("_released", "_corrected"))
    sens["classification_changed"] = sens.outside_band_released != sens.outside_band_corrected
    # The two sensitivity files use different label vocabularies for the same rows: the deposited
    # table says "MLKL only" and "healthy", the recomputed one "MLKL" and "Healthy epithelium".
    # They are mapped explicitly and the join is asserted, because a silent mismatch would report a
    # missing comparison as an unchanged one.
    COMPOSITION = {"MLKL + PGAM5 (figure)": "MLKL + PGAM5", "MLKL only": "MLKL", "PGAM5 only": "PGAM5"}
    LAYER_NAME = {"healthy": "Healthy epithelium", "malignant": "Malignant epithelium",
                  "dome": "3D dome models", "2D": "2D adherent models"}
    rec_r = read_table(released / "IKM_pgam5_sensitivity_recomputed.csv")
    rec_c = read_table(corrected / "IKM_pgam5_sensitivity_recomputed.csv")
    exact = rec_r.merge(rec_c, on=["composition", "gene_count", "layer"], suffixes=("_released", "_corrected"))
    assert len(exact) == len(rec_r) == 12, (len(exact), len(rec_r))
    sens["composition_recomputed"] = sens.composition.map(COMPOSITION)
    sens["layer_recomputed"] = sens.layer.map(LAYER_NAME)
    assert sens.composition_recomputed.notna().all() and sens.layer_recomputed.notna().all()
    exact_cols = ["composition", "gene_count", "layer",
                  "exact_200_gene_band_lower_pp_released", "exact_200_gene_band_upper_pp_released",
                  "outside_exact_200_gene_band_released",
                  "exact_200_gene_band_lower_pp_corrected", "exact_200_gene_band_upper_pp_corrected",
                  "outside_exact_200_gene_band_corrected"]
    sens = sens.merge(exact[[c for c in exact_cols if c in exact.columns]],
                      left_on=["composition_recomputed", "gene_count", "layer_recomputed"],
                      right_on=["composition", "gene_count", "layer"], how="left",
                      suffixes=("", "_drop"))
    sens = sens.drop(columns=[c for c in sens.columns if c.endswith("_drop")])
    one = sens.gene_count == 1
    assert sens.loc[one, "exact_200_gene_band_lower_pp_corrected"].notna().all(), "exact band join failed"
    assert sens.loc[~one, "exact_200_gene_band_lower_pp_corrected"].isna().all(), "two-gene rows have no exact band"
    # The exact enumeration is defined for one-gene sets only, so the two-gene rows carry no value
    # there; a missing value is not a changed classification.
    both = sens.outside_exact_200_gene_band_released.notna() & sens.outside_exact_200_gene_band_corrected.notna()
    sens["exact_classification_changed"] = both & (sens.outside_exact_200_gene_band_released
                                                   != sens.outside_exact_200_gene_band_corrected)
    sens["monte_carlo_and_exact_disagree"] = (sens.outside_exact_200_gene_band_corrected.notna()
                                              & (sens.outside_band_corrected
                                                 != sens.outside_exact_200_gene_band_corrected))
    sens.to_csv(out / "sensitivity_changes.tsv", sep="\t", index=False)
    summary["sensitivity_rows"] = int(len(sens))
    summary["sensitivity_classification_changes"] = int(sens.classification_changed.sum())
    summary["sensitivity_exact_classification_changes"] = int(sens.exact_classification_changed.sum())
    summary["sensitivity_exact_rows"] = int(both.sum())
    summary["monte_carlo_and_exact_disagreements"] = int(sens.monte_carlo_and_exact_disagree.sum())

    claims_r = json.loads((released.parent / "validation/manuscript_values.json").read_text()) \
        if (released.parent / "validation/manuscript_values.json").exists() else None
    claims_c = json.loads((corrected.parent / "validation/manuscript_values.json").read_text())
    if claims_r:
        mrows = []
        for a, b in zip(claims_r["expression"], claims_c["expression"]):
            assert (a["group"], a["column"]) == (b["group"], b["column"])
            mrows.append({"group": a["group"], "column": a["column"],
                          "released": a["value"], "corrected": b["value"],
                          "changed": abs(a["value"] - b["value"]) > 1e-9})
        M = pd.DataFrame(mrows)
        M.to_csv(out / "manuscript_value_changes.tsv", sep="\t", index=False)
        summary["manuscript_values"] = int(len(M))
        summary["manuscript_values_changed"] = int(M.changed.sum())

    json.dump(summary, open(out / "correction_effect_summary.json", "w"), indent=1)
    print(json.dumps(summary, indent=1), flush=True)
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--released", type=Path, required=True)
    parser.add_argument("--corrected", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    sys.exit(main(args.released, args.corrected, args.output_dir))
