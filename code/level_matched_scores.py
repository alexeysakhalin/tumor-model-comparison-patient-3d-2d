"""Expression-rank differences relative to matched control genes.

Main Figure 6 ranks use the shared 18,612-gene universe. Each target gene is
compared with 400 controls matched to its organ-specific median rank in 2D
models. These scores describe relative expression, not pathway activity or
absolute expression. The 2D distribution is approximately centered on zero.
Reference bands and cluster-bootstrap intervals are calculated separately by
build_summary_tables.py.
"""

from __future__ import annotations

import numpy as np
from scipy.stats import rankdata
import pandas as pd


def percentile_ranks(mat: np.ndarray) -> np.ndarray:
    """Within-sample percentile ranks, 0-100, one row per sample.

    Ties share the average of the ranks they span. This matters here and is not a detail: a
    pseudobulk row holds thousands of zeros, and an ordinal rank would hand each of them a
    different percentile decided by the order of the gene axis, so a gene that is simply absent
    could land anywhere inside the zero block. With average ranks every absent gene in a sample
    gets the same percentile, and the value no longer depends on how the matrix is sorted.
    """
    if mat.ndim != 2 or mat.shape[1] == 0 or not np.isfinite(mat).all():
        raise ValueError("Ranks require a finite, nonempty two-dimensional matrix")
    out = np.empty(mat.shape, dtype="float32")
    for i in range(mat.shape[0]):
        out[i] = 100.0 * rankdata(mat[i], method="average") / mat.shape[1]
    return out


def ranks_are_permutation_invariant(mat: np.ndarray, seed: int = 0) -> bool:
    """Shuffle the gene axis, rank, unshuffle: the result must be identical."""
    rng = np.random.default_rng(seed)
    perm = rng.permutation(mat.shape[1])
    back = np.empty_like(perm)
    back[perm] = np.arange(len(perm))
    return np.allclose(percentile_ranks(mat), percentile_ranks(mat[:, perm])[:, back])


def control_sets(
    ref_rank: np.ndarray,
    genes: list[str],
    gene_index: dict,
    exclude: set,
    n_ctrl: int = 400,
) -> dict:
    """For every gene, the indices of the `n_ctrl` genes nearest to it in reference rank.

    Selection is by absolute distance in reference rank, not by a window on the sorted axis: a
    window centred on the insertion position takes whatever lies to the left and right of it, which
    at the edges of the axis and inside large tie blocks is not the nearest set. Ties in distance
    are broken by gene symbol, so the chosen controls do not depend on the order of the gene axis.
    """
    pool = np.array(sorted((i for g, i in gene_index.items() if g not in exclude)))
    if len(pool) < n_ctrl or n_ctrl < 1 or not np.isfinite(ref_rank).all():
        raise ValueError("Insufficient eligible controls or invalid reference ranks")
    missing = set(genes) - set(gene_index)
    if missing:
        raise ValueError(
            f"Target genes missing from the shared axis: {sorted(missing)}"
        )
    names = {i: g for g, i in gene_index.items()}
    pool_names = np.array([names[i] for i in pool])
    ref_pool = ref_rank[pool]
    sets = {}
    for g in genes:
        j = gene_index.get(g)
        if j is None:
            continue
        d = np.abs(ref_pool - ref_rank[j])
        order = np.lexsort((pool_names, d))  # distance first, symbol as the tie-break
        sets[g] = pool[order[:n_ctrl]]
    return sets


def controls_are_nearest(ref_rank, gene_index, exclude, gene, n_ctrl=50) -> bool:
    """The selected controls are the n_ctrl smallest distances available in the pool."""
    sets = control_sets(ref_rank, [gene], gene_index, exclude, n_ctrl)
    pool = np.array(sorted((i for g, i in gene_index.items() if g not in exclude)))
    d_all = np.sort(np.abs(ref_rank[pool] - ref_rank[gene_index[gene]]))[:n_ctrl]
    d_sel = np.sort(np.abs(ref_rank[sets[gene]] - ref_rank[gene_index[gene]]))
    return bool(np.allclose(d_all, d_sel))


def control_match_report(
    ref_rank: np.ndarray, sets: dict, gene_index: dict
) -> "pd.DataFrame":
    """How closely each control set matches its gene in reference rank."""
    import pandas as pd

    rows = []
    for g, cols in sets.items():
        j = gene_index[g]
        d = np.abs(ref_rank[cols] - ref_rank[j])
        rows.append(
            {
                "gene": g,
                "ref_rank": float(ref_rank[j]),
                "n_controls": len(cols),
                "median_abs_rank_gap": float(np.median(d)),
                "max_abs_rank_gap": float(d.max()),
            }
        )
    return pd.DataFrame(rows)


def sample_scores(
    rank_mat: np.ndarray, rows: np.ndarray, gene_col: int, ctrl_cols: np.ndarray
) -> np.ndarray:
    """Per-sample score of one gene against its control set."""
    return rank_mat[rows, gene_col] - np.median(
        rank_mat[np.ix_(rows, ctrl_cols)], axis=1
    )


def group_table(
    genes: list[str],
    organs: list[str],
    layers: dict,
    ref_layer: str,
    gene_index: dict,
    organ_of: dict,
    n_ctrl: int = 400,
    exclude: set | None = None,
) -> pd.DataFrame:
    """Per-gene, per-organ score in every layer.

    layers: {layer_name: rank matrix}. ref_layer: the layer used for level matching (flat lines).
    gene_index: {symbol: column} shared by the matrices. organ_of: {layer_name: organ vector}.
    Every layer must use the same gene axis. This function summarizes individual
    genes over profiles; group-level sample-first estimates are built separately.
    """
    excl = set(genes) | (exclude or set())
    out = []
    for organ in organs:
        ref_rows = np.where(organ_of[ref_layer] == organ)[0]
        if len(ref_rows) < 10:
            raise ValueError(f"{organ}: fewer than 10 reference profiles")
        ref_rank = np.median(layers[ref_layer][ref_rows, :], axis=0)
        csets = control_sets(ref_rank, genes, gene_index, excl, n_ctrl)
        for g, ctrl in csets.items():
            rec = {"gene": g, "organ": organ}
            for lay, mat in layers.items():
                rows = np.where(organ_of[lay] == organ)[0]
                rec[lay] = (
                    float(np.median(sample_scores(mat, rows, gene_index[g], ctrl)))
                    if len(rows)
                    else np.nan
                )
            out.append(rec)
    return pd.DataFrame(out)
