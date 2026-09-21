"""Aggregate the Kang pan-cancer atlas files to per (sample, cell class) profiles.

These files carry log-normalised values, not counts, so groups are summarised by the MEAN of the
stored values rather than a sum: summing log-normalised values is not a count and would scale with
the number of cells. The mean is written as its own modality and must not be mixed with the
count-based pseudobulk of the CELLxGENE atlases.

Malignant epithelium is identified by the file's own cnv_status field rather than by a cell-type
label, because this atlas has no malignant label. Cells are read in row blocks straight out of
HDF5, so peak memory stays at one block.
"""

import os
import sys

import h5py
import numpy as np
import pandas as pd

KEEP = {"Epithelial", "T cell", "NK cell", "Macrophage"}
CARRY = ["Cancer type", "Dataset", "Organ_origin", "Patient", "Tissue"]
BLOCK = 20000


def obs_col(f, name):
    node = f["obs"][name]
    if isinstance(node, h5py.Group) and "categories" in node:
        cats = np.array(
            [
                c.decode() if isinstance(c, bytes) else str(c)
                for c in node["categories"][:]
            ]
        )
        codes = node["codes"][:]
        return np.where(codes >= 0, cats[np.clip(codes, 0, len(cats) - 1)], "nan")
    vals = node[:]
    if vals.dtype.kind in "SO":
        return np.array([v.decode() if isinstance(v, bytes) else str(v) for v in vals])
    return vals.astype(str)


def var_names(f):
    for path in ("var/_index", "var/index", "var/gene_symbols", "var/features"):
        if path in f:
            raw = f[path][:]
            return np.array(
                [v.decode() if isinstance(v, bytes) else str(v) for v in raw]
            )
    raise KeyError("no var index")


def one_file(path, out_dir):
    tag = os.path.basename(path)[:-5]
    with h5py.File(path, "r") as f:
        genes = var_names(f)
        ct = obs_col(f, "Celltype")
        cnv = obs_col(f, "cnv_status")
        samp = obs_col(f, "Sample")
        label = np.where(ct == "Epithelial", np.char.add("Epithelial_", cnv), ct)
        keep = np.isin(ct, list(KEEP))
        if keep.sum() == 0:
            return None
        group = np.char.add(np.char.add(samp.astype(str), "|"), label.astype(str))
        uniq, inv = np.unique(group[keep], return_inverse=True)
        idx_keep = np.flatnonzero(keep)
        acc = np.zeros((len(uniq), len(genes)), dtype="float32")
        cnt = np.bincount(inv, minlength=len(uniq)).astype("float32")
        g = f["X"]
        n_cells = len(ct)
        pos_of = np.full(n_cells, -1, dtype="int64")
        pos_of[idx_keep] = inv
        if isinstance(g, h5py.Group):
            encoding = g.attrs.get("encoding-type", "")
            if isinstance(encoding, bytes):
                encoding = encoding.decode()
            if encoding != "csr_matrix" or tuple(g.attrs["shape"]) != (
                n_cells,
                len(genes),
            ):
                raise ValueError(
                    "A CSR matrix aligned with the observation and gene axes is required"
                )
            indptr = g["indptr"][:]
            for start in range(0, n_cells, BLOCK):
                stop = min(start + BLOCK, n_cells)
                lo, hi = indptr[start], indptr[stop]
                if hi == lo:
                    continue
                data = g["data"][lo:hi]
                if not np.isfinite(data).all():
                    raise ValueError("Non-finite expression values")
                indices = g["indices"][lo:hi]
                rows = np.repeat(
                    np.arange(start, stop), np.diff(indptr[start : stop + 1])
                )
                sel = pos_of[rows] >= 0
                if not sel.any():
                    continue
                np.add.at(acc, (pos_of[rows[sel]], indices[sel]), data[sel])
        else:
            if g.shape != (n_cells, len(genes)):
                raise ValueError("Matrix and annotation axes do not match")
            for start in range(0, n_cells, BLOCK):
                stop = min(start + BLOCK, n_cells)
                blk = np.asarray(g[start:stop], dtype="float32")
                if not np.isfinite(blk).all():
                    raise ValueError("Non-finite expression values")
                loc = pos_of[start:stop]
                for i in np.flatnonzero(loc >= 0):
                    acc[loc[i]] += blk[i]
        acc /= cnt[:, None]
        meta = pd.DataFrame({"group": uniq, "n_cells": cnt.astype(int)})
        meta[["sample", "cell_class"]] = meta.group.str.split("|", n=1, expand=True)
        first = pd.DataFrame({c: obs_col(f, c) for c in CARRY if c in f["obs"]})
        first["sample"] = samp
        first = first.groupby("sample").first()
        meta = meta.merge(first, left_on="sample", right_index=True, how="left")
        meta["study_file"] = tag
    np.save(os.path.join(out_dir, f"{tag}_mean.npy"), acc)
    pd.DataFrame({"gene": genes}).to_csv(
        os.path.join(out_dir, f"{tag}_genes.csv"), index=False
    )
    meta.to_csv(os.path.join(out_dir, f"{tag}_groups.tsv"), sep="\t", index=False)
    return tag, acc.shape, len(meta)


if __name__ == "__main__":
    src, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    files = sorted(f for f in os.listdir(src) if f.endswith(".h5ad"))
    if not files:
        raise SystemExit("No H5AD inputs found")
    done = 0
    errors = []
    for fn in files:
        tag = fn[:-5]
        try:
            res = one_file(os.path.join(src, fn), out)
            done += 1
            print(
                f"[{done}/{len(files)}] {tag}: {res[1] if res else 'no target cells'} groups={res[2] if res else 0}",
                flush=True,
            )
        except Exception as exc:
            errors.append(tag)
            print(f"[fail] {tag}: {type(exc).__name__}: {exc}", flush=True)
    print("finished:", done, "of", len(files), flush=True)
    if errors:
        raise SystemExit(f"Failed inputs: {', '.join(errors)}")
