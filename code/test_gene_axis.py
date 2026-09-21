"""Regression tests for the gene axis of pseudobulk_stream.py.

The defect these guard against: the matrix source was chosen with an identity test on an h5py
object, which is always False, so counts read from raw/X could be written out under the annotation
of the other var group. Each test builds a small H5AD in a temporary directory, runs the script and
checks the identifier and the count that come out.

Run: python test_gene_axis.py
"""

import os
import subprocess
import sys
import tempfile

import h5py
import numpy as np
from table_io import read_table
import scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__))


def write_h5ad(
    path, counts, var_index, raw_counts=None, raw_var_index=None, cell_types=None
):
    """Minimal H5AD: X (and optionally raw/X), var/_index, obs with two columns."""
    n_obs = counts.shape[0]
    with h5py.File(path, "w") as f:
        for grp, mat, idx in [
            ("X", counts, var_index),
            ("raw/X", raw_counts, raw_var_index),
        ]:
            if mat is None:
                continue
            m = sp.csr_matrix(mat)
            g = f.create_group(grp)
            g.attrs["shape"] = np.array(m.shape, dtype="int64")
            g.attrs["encoding-type"] = "csr_matrix"
            g.create_dataset("data", data=m.data.astype("float32"))
            g.create_dataset("indices", data=m.indices.astype("int32"))
            g.create_dataset("indptr", data=m.indptr.astype("int64"))
            vg = f.create_group("var" if grp == "X" else "raw/var")
            vg.create_dataset("_index", data=np.array(idx, dtype="S32"))
        o = f.create_group("obs")
        o.create_dataset(
            "_index", data=np.array([f"c{i}" for i in range(n_obs)], dtype="S8")
        )
        o.create_dataset("donor_id", data=np.array(["d1"] * n_obs, dtype="S8"))
        o.create_dataset(
            "compartment",
            data=np.array(cell_types or ["Epithelium"] * n_obs, dtype="S24"),
        )


def run(infile, prefix):
    r = subprocess.run(
        [
            sys.executable,
            os.path.join(HERE, "pseudobulk_stream.py"),
            infile,
            prefix,
            "--group-col",
            "donor_id",
            "--celltype-col",
            "compartment",
            "--keep-celltypes",
            "Epithelium",
            "--block",
            "8",
            "--matrix-is-counts",
        ],
        capture_output=True,
        text=True,
    )
    if r.returncode != 0:
        raise AssertionError(
            r.stderr.strip().splitlines()[-1] if r.stderr else "run failed"
        )
    return np.load(prefix + "_counts.npy"), read_table(prefix + "_genes.csv")


def case_permuted_axes(tmp):
    """raw/var and var hold the same genes in different order: counts must follow raw/var."""
    p = os.path.join(tmp, "permuted.h5ad")
    write_h5ad(
        p,
        counts=np.array([[0.0, 20.0, 10.0]]),
        var_index=["C", "B", "A"],
        raw_counts=np.array([[10.0, 20.0, 0.0]]),
        raw_var_index=["A", "B", "C"],
    )
    counts, genes = run(p, os.path.join(tmp, "permuted_out"))
    got = dict(zip(genes["gene_id"], counts[:, 0]))
    assert got == {"A": 10.0, "B": 20.0, "C": 0.0}, (
        f"counts landed under the wrong genes: {got}"
    )
    return "counts follow raw/var when the two var groups are permuted"


def case_shorter_var(tmp):
    """var holds fewer genes than raw/var: the run must not mix the two axes."""
    p = os.path.join(tmp, "shorter.h5ad")
    write_h5ad(
        p,
        counts=np.array([[1.0, 2.0]]),
        var_index=["A", "B"],
        raw_counts=np.array([[10.0, 20.0, 30.0]]),
        raw_var_index=["A", "B", "C"],
    )
    counts, genes = run(p, os.path.join(tmp, "shorter_out"))
    assert len(genes) == counts.shape[0] == 3, (
        f"{len(genes)} names for {counts.shape[0]} rows"
    )
    assert list(genes["gene_id"]) == ["A", "B", "C"], list(genes["gene_id"])
    return "a shorter var group cannot be used to label the raw matrix"


def case_no_raw(tmp):
    """No raw group at all: X is summed and labelled by var."""
    p = os.path.join(tmp, "noraw.h5ad")
    write_h5ad(p, counts=np.array([[5.0, 0.0, 7.0]]), var_index=["A", "B", "C"])
    counts, genes = run(p, os.path.join(tmp, "noraw_out"))
    got = dict(zip(genes["gene_id"], counts[:, 0]))
    assert got == {"A": 5.0, "B": 0.0, "C": 7.0}, got
    return "X is used and labelled by var when no raw group exists"


def case_plain_feature_name(tmp):
    """feature_name stored as a plain string dataset: symbols must survive into the gene table."""
    p = os.path.join(tmp, "plain.h5ad")
    write_h5ad(p, counts=np.array([[1.0, 2.0]]), var_index=["ENSG1", "ENSG2"])
    with h5py.File(p, "a") as f:
        f["var"].create_dataset(
            "feature_name", data=np.array(["AAA", "BBB"], dtype="S8")
        )
    counts, genes = run(p, os.path.join(tmp, "plain_out"))
    assert "symbol" in genes.columns, f"symbols dropped: {list(genes.columns)}"
    assert list(genes["symbol"]) == ["AAA", "BBB"], list(genes["symbol"])
    return "a plain string feature_name is read as symbols"


def case_groups_x_only(tmp):
    """pseudobulk_groups.py must run on a file without a raw group."""
    p = os.path.join(tmp, "grouponly.h5ad")
    write_h5ad(p, counts=np.array([[1.0, 2.0], [3.0, 4.0]]), var_index=["A", "B"])
    r = subprocess.run(
        [
            sys.executable,
            os.path.join(HERE, "pseudobulk_groups.py"),
            p,
            os.path.join(tmp, "grouponly_out"),
            "--group-col",
            "donor_id",
            "--celltype-col",
            "compartment",
        ],
        capture_output=True,
        text=True,
    )
    assert r.returncode == 0, (r.stderr or "").strip().splitlines()[-1]
    return "the metadata stage runs on an X-only file"


if __name__ == "__main__":
    fails = []
    with tempfile.TemporaryDirectory() as tmp:
        for case in (
            case_permuted_axes,
            case_shorter_var,
            case_no_raw,
            case_plain_feature_name,
            case_groups_x_only,
        ):
            try:
                print("[ok  ]", case(tmp))
            except AssertionError as e:
                fails.append(f"{case.__name__}: {e}")
                print("[FAIL]", case.__name__, "—", e)
    print("\nFAILED:", ", ".join(fails) if fails else "none")
    sys.exit(1 if fails else 0)
