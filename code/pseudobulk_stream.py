"""Pseudobulk a CELLxGENE .h5ad by streaming its CSR blocks straight out of HDF5.

Sparse rows are read in blocks; memory also depends on the number of genes,
groups and cells in the metadata. Raw counts are read from raw/X when present.
Using X requires an explicit assertion that it contains counts. Graphs are not read.

Outputs: <prefix>_counts.npy (genes x groups, float32), <prefix>_genes.csv, <prefix>_groups.tsv.

Usage: python pseudobulk_stream.py <in.h5ad> <out_prefix> [--keep-celltypes "a;b"] [--block N]
"""

import argparse

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

CARRY = [
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
        out = np.where(codes >= 0, cats[np.clip(codes, 0, len(cats) - 1)], "nan")
        return out
    vals = node[:]
    if vals.dtype.kind in "SO":
        return np.array([v.decode() if isinstance(v, bytes) else str(v) for v in vals])
    return vals.astype(str)


def decode(arr):
    return np.array([v.decode() if isinstance(v, bytes) else str(v) for v in arr])


def main():
    p = argparse.ArgumentParser()
    p.add_argument("infile")
    p.add_argument("out_prefix")
    p.add_argument("--group-col", default="sample_id")
    p.add_argument("--celltype-col", default="cell_type")
    p.add_argument("--keep-celltypes", default=None)
    p.add_argument("--block", type=int, default=20000)
    p.add_argument(
        "--matrix-is-counts",
        action="store_true",
        help="Confirm from source metadata that X contains counts if raw/X is absent",
    )
    p.add_argument(
        "--require",
        action="append",
        default=[],
        help="obs filter as col=value, repeatable. Cells failing any filter are dropped, "
        "so two sequencing platforms are never summed into one group.",
    )
    args = p.parse_args()

    f = h5py.File(args.infile, "r")
    # The path is decided first and the dataset is opened from it. Asking h5py for the same object
    # twice returns two different Python objects, so an identity test (`xg is f.get("raw/X")`) is
    # always False and would label raw counts with the annotation of the other var group.
    source = "raw/X" if ("raw" in f and "X" in f["raw"]) else "X"
    if source == "X" and not args.matrix_is_counts:
        raise ValueError(
            "X cannot be assumed to contain counts; confirm the source modality"
        )
    if args.block < 1:
        raise ValueError("Block size must be positive")
    xg = f[source]
    encoding = xg.attrs.get("encoding-type", "")
    if isinstance(encoding, bytes):
        encoding = encoding.decode()
    if not isinstance(xg, h5py.Group) or encoding != "csr_matrix":
        raise ValueError("This count aggregator requires a CSR-encoded matrix")
    n_obs, n_v = (int(x) for x in xg.attrs["shape"])
    indptr = xg["indptr"][:]

    samples = read_obs_col(f, args.group_col)
    ctypes = read_obs_col(f, args.celltype_col)
    if len(samples) != n_obs or len(ctypes) != n_obs:
        raise ValueError("Observation metadata do not match the matrix row axis")
    if args.keep_celltypes:
        keep = set(args.keep_celltypes.split(";"))
        mask = np.isin(ctypes, list(keep))
    else:
        mask = np.ones(n_obs, dtype=bool)
    for req in args.require:
        col, val = req.split("=", 1)
        mask &= read_obs_col(f, col) == val
    key = np.char.add(np.char.add(samples, "||"), ctypes)
    groups = np.unique(key[mask])
    if not len(groups) or np.isin(samples[mask], ["", "nan", "None"]).any():
        raise ValueError("No selected cells or missing grouping identifiers")
    gpos = {g: i for i, g in enumerate(groups)}
    codes = np.array([gpos.get(k, -1) for k in key], dtype="int64")
    codes[~mask] = -1
    cells = np.bincount(codes[codes >= 0], minlength=len(groups))
    print(
        f"cells {n_obs} | selected {int(mask.sum())} | groups {len(groups)} | genes {n_v} | matrix {source}",
        flush=True,
    )

    acc = np.zeros((len(groups), n_v), dtype="float32")
    for start in range(0, n_obs, args.block):
        stop = min(start + args.block, n_obs)
        sel = codes[start:stop]
        rows = np.where(sel >= 0)[0]
        if rows.size == 0:
            continue
        lo, hi = int(indptr[start]), int(indptr[stop])
        values = xg["data"][lo:hi].astype("float32")
        if (
            not np.isfinite(values).all()
            or np.any(values < 0)
            or not np.allclose(values, np.rint(values), atol=1e-6, rtol=0)
        ):
            raise ValueError(
                "Count aggregation encountered invalid or non-integer values"
            )
        block = sp.csr_matrix(
            (values, xg["indices"][lo:hi], indptr[start : stop + 1] - lo),
            shape=(stop - start, n_v),
        )
        ind = sp.csr_matrix(
            (np.ones(rows.size, dtype="float32"), (sel[rows], rows)),
            shape=(len(groups), stop - start),
        )
        acc += (ind @ block).toarray()
        if (start // args.block) % 25 == 0:
            print(f"rows {stop}/{n_obs}", flush=True)

    # Identifiers and symbols must both come from the var of the matrix that was actually summed.
    # In a CELLxGENE file raw/var and var can differ in length and in order, so taking ids from one
    # and symbols from the other would label a gene's counts with another gene's name. Both are read
    # from the same group and the length is asserted against the matrix.
    vgrp = "raw/var" if source == "raw/X" else "var"
    gene_ids = decode(f[f"{vgrp}/_index"][:])
    if len(gene_ids) != n_v:
        raise SystemExit(
            f"{vgrp}/_index has {len(gene_ids)} entries but {source} has {n_v} columns"
        )
    genes = pd.DataFrame({"gene_id": gene_ids})
    if "feature_name" in f[vgrp]:
        fn = f[vgrp]["feature_name"]
        # Two encodings occur: a categorical group (categories + codes) and a plain string dataset.
        # Reading only the categorical form silently drops the symbols of a file that uses the other,
        # leaving a gene table with identifiers alone.
        if isinstance(fn, h5py.Group):
            cats = decode(fn["categories"][:])
            codes = fn["codes"][:]
            if np.any(codes < 0):
                raise ValueError("Missing feature names in the selected gene axis")
            sym = cats[codes]
        else:
            sym = decode(fn[:])
        if len(sym) != n_v:
            raise SystemExit(
                f"{vgrp}/feature_name has {len(sym)} entries but {source} has {n_v} columns"
            )
        genes["symbol"] = sym
    elif (
        source == "raw/X"
        and "var" in f
        and "feature_name" in f["var"]
        and len(decode(f["var/_index"][:])) == n_v
    ):
        # Same gene axis in both var groups: the symbols of var can be used without misalignment.
        fn = f["var"]["feature_name"]
        if (
            isinstance(fn, h5py.Group)
            and (decode(f["var/_index"][:]) == gene_ids).all()
        ):
            cats = decode(fn["categories"][:])
            genes["symbol"] = cats[fn["codes"][:]]
    np.save(f"{args.out_prefix}_counts.npy", acc.T)
    genes.to_csv(f"{args.out_prefix}_genes.csv", index=False)

    meta = pd.DataFrame({"group": groups, "n_cells": cells})
    meta[["sample", "cell_type"]] = meta.group.str.split(r"\|\|", n=1, expand=True)
    carry = {c: read_obs_col(f, c) for c in CARRY if c in f["obs"]}
    if carry:
        per_cell = pd.DataFrame({"sample": samples, **carry})[mask]
        for field in ["donor_id", "assay", "tissue"]:
            if (
                field in per_cell
                and per_cell.groupby("sample")[field].nunique().gt(1).any()
            ):
                raise ValueError(f"Conflicting {field} within a pseudobulk sample")
        meta = meta.merge(
            per_cell.groupby("sample").first(),
            left_on="sample",
            right_index=True,
            how="left",
        )
    meta.to_csv(f"{args.out_prefix}_groups.tsv", sep="\t", index=False)
    print("written:", acc.T.shape, "| groups table:", meta.shape, flush=True)


if __name__ == "__main__":
    main()
