"""Verify release integrity and Figure 6 numerical and graphical consistency.

Use --build-dir after run_all.sh has recalculated expression statistics and
rendered the figure. Missing or inconsistent outputs produce a nonzero exit.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

from recompute_pgam5_sensitivity import compute
from table_io import read_table

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition, label):
    if not condition:
        raise AssertionError(label)
    print(f"PASS {label}", flush=True)


def verify(build_dir=None):
    manifest = pd.read_csv(ROOT / "MANIFEST_SHA256.tsv", sep="\t")
    require(not manifest.file.duplicated().any(), "Unique manifest paths")
    require(
        all(
            (ROOT / r.file).is_file() and sha(ROOT / r.file) == r.sha256
            for r in manifest.itertuples()
        ),
        "File checksums",
    )
    snapshot = json.loads((ROOT / "ZENODO_TABLE_SNAPSHOT.json").read_text())
    require(
        len(snapshot["files"]) == 51
        and all(
            sha(ROOT / name) == checksum for name, checksum in snapshot["files"].items()
        ),
        "All 51 tables match the source-data archive",
    )
    mapping = pd.read_csv(ROOT / "FIGURE6_PANEL_MAP.tsv", sep="\t")
    require(
        len(mapping) == 6
        and all(
            sha(ROOT / r.source_data_file)
            == sha(ROOT / r.analysis_table)
            == r.source_data_sha256
            for r in mapping.itertuples()
        ),
        "Six plotted source tables match analysis tables",
    )

    tables = ROOT / "tables"
    members = read_table(tables / "IKM_gene_membership_by_group.tsv", sep="\t")
    require(
        len(members) == members.gene.nunique() == 67 and members.group.nunique() == 12,
        "67 unique genes in 12 non-overlapping groups",
    )
    require(
        members.groupby("group").figure_label.nunique().eq(1).all(),
        "One scientific display label per group",
    )
    panel = read_table(tables / "panoptosis_panel_2026-09-18.tsv", sep="\t")
    require(
        len(set(panel.gene) & set(members.gene)) == 41
        and len(set(members.gene) - set(panel.gene)) == 26,
        "41 panel genes and 26 additional genes",
    )
    ps = read_table(tables / "IKM_group_scores_by_sample.csv")
    require(
        ps.groupby("group").size().eq(804).all() and np.isfinite(ps.score).all(),
        "804 finite profiles per group",
    )
    require(
        not ps.duplicated(["layer", "organ", "sample_id", "group"]).any(),
        "Unique sample-group keys",
    )
    unique = ps.drop_duplicates(["layer", "organ", "sample_id"])
    require(
        unique.groupby("layer").size().to_dict()
        == {"healthy": 20, "malignant": 231, "dome": 152, "2D": 401},
        "Four cohort sizes",
    )
    require(
        unique.groupby("layer").organ.nunique().eq(6).all(), "Six organs in each cohort"
    )

    screen = read_table(tables / "IKM_models_panoptosis_layers.csv")
    require(
        len(screen) == 84
        and screen.model.nunique() == 7
        and screen.groupby("model").group.nunique().eq(12).all()
        and not screen.duplicated(["model", "group"]).any(),
        "Seven complete CRISPR arms",
    )
    half_width = (screen.band_hi - screen.band_lo) / 2
    require(
        (half_width > 0).all()
        and np.isfinite(screen[["score", "reference_band_units"]]).all().all(),
        "Finite functional effects and positive reference widths",
    )
    require(
        np.allclose(
            screen.score / half_width, screen.reference_band_units, atol=0.025, rtol=0
        ),
        "84 normalized effects agree within source rounding",
    )
    require(
        np.array_equal(
            (screen.score < screen.band_lo) | (screen.score > screen.band_hi),
            screen.outside_band,
        ),
        "Functional reference-band classifications",
    )
    hela = screen[screen.model.str.startswith("HeLa")]
    require(
        hela.replicate_rho.isna().all() and hela.screen_count.eq(1).all(),
        "HeLa library halves are not counted as biological replicates",
    )
    weights = read_table(tables / "IKM_group_weight_all12.csv").set_index("group")
    require(
        len(weights) == 12
        and all(
            abs(
                round(float(d.reference_band_units.abs().median()), 2)
                - weights.loc[g, "magnitude"]
            )
            < 0.011
            and int(d.outside_band.sum()) == int(weights.loc[g, "arms_outside_band"])
            for g, d in screen.groupby("group")
        ),
        "Twelve functional magnitude summaries",
    )
    ranked = read_table(tables / "IKM_patel2017_ranked_genes.csv")
    overlap = read_table(tables / "IKM_patel2017_group_overrepresentation.csv")
    tail = set(ranked.gene_id)
    require(
        len(overlap) == 12
        and all(
            len(set(members.loc[members.group == r.group, "gene"]) & tail)
            == r.published_list_overlap
            for r in overlap.itertuples()
        ),
        "Twelve published-list overlaps",
    )
    pgam5, reference, _ = compute(tables)
    for name, frame in [
        ("IKM_pgam5_sensitivity_recomputed.csv", pgam5),
        ("IKM_pgam5_reference_gene_scores.csv", reference),
    ]:
        pd.testing.assert_frame_equal(
            frame, read_table(tables / name), check_dtype=False, rtol=1e-8, atol=1e-10
        )
    require(True, "MLKL/PGAM5 estimates, reference distributions and 24 screen values")

    ci = read_table(tables / "IKM_final_master_table_CI.csv").set_index("group")
    claims = json.loads((ROOT / "validation/manuscript_values.json").read_text())
    for row in claims["expression"]:
        require(
            abs(float(ci.loc[row["group"], row["column"]]) - row["value"]) < 1e-9,
            f"Manuscript expression value: {row['group']}, {row['column']}",
        )

    if build_dir is not None:
        for name in [
            "IKM_null_bands.csv",
            "IKM_final_master_table_CI.csv",
            "IKM_representation_in_band_units.csv",
            "IKM_final_master_table.csv",
            "IKM_models_vs_patient.csv",
            "IKM_bootstrap_rejection.csv",
            "IKM_single_cluster_strata.csv",
        ]:
            rebuilt = read_table(build_dir / "tables" / name)
            expected = read_table(tables / name)
            pd.testing.assert_frame_equal(
                rebuilt,
                expected,
                check_dtype=False,
                check_exact=False,
                rtol=1e-10,
                atol=1e-10,
            )
            require(True, f"Recalculated {name}")
        require(
            (build_dir / "Figure6.pdf").is_file()
            and (build_dir / "Figure6.pdf").stat().st_size > 0,
            "Rebuilt PDF exists",
        )
        with (
            Image.open(ROOT / "figures/Figure6.png") as source,
            Image.open(build_dir / "Figure6.png") as rebuilt,
        ):
            require(
                np.array_equal(
                    np.asarray(source.convert("RGBA")),
                    np.asarray(rebuilt.convert("RGBA")),
                ),
                "Rebuilt figure is pixel-identical to the accepted Figure 6",
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path)
    arguments = parser.parse_args()
    verify(arguments.build_dir)
    print("All requested checks passed.")
