"""Build the per-group metadata table for an already-computed pseudobulk matrix.

Only obs is read, and every carried column is taken at the first cell of each sample rather than
by building a per-cell frame, so memory stays at a few tens of megabytes instead of gigabytes.
"""

import argparse

import h5py
import numpy as np
import pandas as pd

CARRY = [
    "dataset",
    "origin",
    "ever_smoker",
    "EGFR_mutation",
    "KRAS_mutation",
    "TP53_mutation",
    "study_id",
    "donor_id",
    "sample_type",
    "disease",
    "tissue",
    "assay",
    "platform",
    "medical_condition",
    "cancer_type",
    "tumor_stage",
    "microsatellite_status",
    "CMS_type",
    "immune_infiltration_type",
    "treatment_status_before_resection",
    "treatment_response",
    "GEO_sample_accession",
    "SRA_sample_accession",
    "NCBI_BioProject_accession",
    "matrix_type",
    "tissue_cell_state",
    "anatomic_location",
    "sex",
    "age",
    "tumor_source",
]


def read_obs_col(f, name):
    """Return an obs column as a string array, handling categorical and plain encodings."""
    node = f["obs"][name]
    if isinstance(node, h5py.Group):
        cats = node["categories"][:]
        cats = np.array([c.decode() if isinstance(c, bytes) else str(c) for c in cats])
        codes = node["codes"][:]
        return np.where(codes >= 0, cats[np.clip(codes, 0, len(cats) - 1)], "nan")
    vals = node[:]
    if vals.dtype.kind in "SO":
        return np.array([v.decode() if isinstance(v, bytes) else str(v) for v in vals])
    return vals.astype(str)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("infile")
    p.add_argument("out_prefix")
    p.add_argument("--group-col", default="sample_id")
    p.add_argument("--celltype-col", default="cell_type")
    p.add_argument("--keep-celltypes", default=None)
    p.add_argument(
        "--require",
        action="append",
        default=[],
        help="obs filter as col=value, repeatable. MUST be given exactly as in "
        "pseudobulk_stream.py: this script decides which groups exist, and a filter "
        "applied to the matrix but not here produces a metadata table with more rows "
        "than the matrix has columns.",
    )
    p.add_argument(
        "--expect-columns",
        type=int,
        default=None,
        help="number of columns of the matrix this table describes; the run fails if the "
        "group count differs.",
    )
    args = p.parse_args()

    f = h5py.File(args.infile, "r")
    samples = read_obs_col(f, args.group_col)
    ctypes = read_obs_col(f, args.celltype_col)
    mask = (
        np.isin(ctypes, args.keep_celltypes.split(";"))
        if args.keep_celltypes
        else np.ones(len(ctypes), bool)
    )
    for req in args.require:
        col, val = req.split("=", 1)
        mask &= read_obs_col(f, col) == val

    key = np.char.add(np.char.add(samples, "||"), ctypes)
    groups, codes_g = np.unique(key[mask], return_inverse=True)
    cells = np.bincount(codes_g, minlength=len(groups))

    s_uniq, s_inv = np.unique(samples, return_inverse=True)
    pos = np.arange(len(samples))
    first = np.full(len(s_uniq), len(samples), dtype="int64")
    np.minimum.at(first, s_inv[mask], pos[mask])
    valid = first < len(samples)
    per_sample = {"sample": s_uniq[valid]}
    for c in CARRY:
        if c in f["obs"]:
            per_sample[c] = read_obs_col(f, c)[first[valid]]
    meta = pd.DataFrame({"group": groups, "n_cells": cells})
    meta[["sample", "cell_type"]] = meta.group.str.split(r"\|\|", n=1, expand=True)
    meta = meta.merge(pd.DataFrame(per_sample), on="sample", how="left")
    if args.expect_columns is not None and len(meta) != args.expect_columns:
        raise SystemExit(
            f"groups table has {len(meta)} rows but the matrix has "
            f"{args.expect_columns} columns: the filters of the two scripts differ"
        )
    meta.to_csv(f"{args.out_prefix}_groups.tsv", sep="\t", index=False)
    print(
        "groups table:",
        meta.shape,
        "| samples:",
        meta["sample"].nunique(),
        "| cell types:",
        meta.cell_type.nunique(),
        flush=True,
    )

    # An X-only file has no raw group; the comparison is simply not available there, and raising
    # a KeyError on a valid file would stop the metadata stage for no reason.
    raw_ids = (
        np.array(
            [
                v.decode() if isinstance(v, bytes) else str(v)
                for v in f["raw/var/_index"][:]
            ]
        )
        if "raw" in f and "var" in f["raw"]
        else None
    )
    x_ids = np.array(
        [v.decode() if isinstance(v, bytes) else str(v) for v in f["var/_index"][:]]
    )
    print(
        "raw/var identical to var:",
        "no raw group" if raw_ids is None else bool(np.array_equal(raw_ids, x_ids)),
        flush=True,
    )


if __name__ == "__main__":
    main()
