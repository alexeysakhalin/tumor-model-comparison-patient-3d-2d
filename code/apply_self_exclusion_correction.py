"""Refresh the deposited tables that the background self-exclusion correction changes.

`control_sets` excludes every scored gene from its own control set. Only the 200 background genes
were affected, because `group_table` already adds the analysed genes to the exclusion set, so the
analysed estimates, their intervals and the functional-screen inputs keep their values. What has to
be regenerated is the chain that starts at the background scores: the reference bands, the tables
built from them by `build_summary_tables.py`, and the MLKL/PGAM5 sensitivity table, whose one-gene
bands come from the background and whose two-gene band is the size-two band of the main figure.

This script writes only the sensitivity table; the rest is produced by the normal pipeline stages and
copied into `tables/`. It keeps the deposited column names and row order, edits numbers in place, and
prints what changed. Band limits are stored to one decimal place and the normalisation divides the
full-precision score by half the width of the stored band, which is the convention of the release.

Usage: python code/apply_self_exclusion_correction.py tables [--baseline DIR]
"""
import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from recompute_pgam5_sensitivity import LAYERS, observed_scores, reference_bands
from table_io import aliases, english_label, read_table

def deposited_columns(frame_columns):
    """Map the English schema names onto the deposited header, using the versioned aliases.

    The deposited tables keep the column names of the data record, and `table_io` translates them on
    read. Writing has to address the same columns without naming them literally, so the alias file is
    inverted here: a script stays in English and the data record keeps its bytes.
    """
    inverse = {}
    for source, target in aliases().items():
        inverse.setdefault(target, source)
    lookup = {}
    for column in frame_columns:
        lookup[english_label(column)] = column
    return lookup, inverse


def main(tables, baseline=None):
    tables = Path(tables)
    bands, _, genes = reference_bands(tables)
    observed, pair = observed_scores(tables)
    two = read_table(tables / "IKM_null_bands.csv")
    two = two[two["group_size"] == 2].set_index("layer")

    path = tables / "IKM_pgam5_sensitivity.csv"
    raw = pd.read_csv(path)
    before = raw.copy()
    col, _ = deposited_columns(raw.columns)
    english = read_table(path)
    single = {"MLKL only": "MLKL", "PGAM5 only": "PGAM5"}
    for i in raw.index:
        layer = english.loc[i, "layer"]
        if int(english.loc[i, "gene_count"]) == 2:
            low, high = float(two.loc[layer, "band_lower"]), float(two.loc[layer, "band_upper"])
        else:
            low, high = (float(x) for x in bands[layer][0])
        low, high = round(low, 1), round(high, 1)
        # The stored value is rounded to one decimal place; the normalisation and the outside-band
        # test use the full-precision statistic, as the release convention does.
        composition = english.loc[i, "composition"]
        score = (float(pair.loc[layer]) if composition not in single
                 else float(observed.loc[single[composition], layer]))
        assert round(score, 1) == float(english.loc[i, "score_pp"]), (composition, layer, score)
        raw.loc[i, col["band_lower"]] = low
        raw.loc[i, col["band_upper"]] = high
        raw.loc[i, col["outside_band"]] = bool(score < low or score > high)
        raw.loc[i, col["reference_band_units"]] = round(score / ((high - low) / 2), 2)
    raw.to_csv(path, index=False)

    ref = pd.read_csv(Path(baseline) / "IKM_pgam5_sensitivity.csv") if baseline else before
    changed = []
    for i in raw.index:
        for name in ("band_lower", "band_upper", "outside_band", "reference_band_units"):
            column = col[name]
            if str(ref.loc[i, column]) != str(raw.loc[i, column]):
                changed.append((english.loc[i, "composition"], english.loc[i, "layer"], name,
                                ref.loc[i, column], raw.loc[i, column]))
    print(f"sensitivity rows: {len(raw)} | fields changed: {len(changed)}", flush=True)
    for comp, layer, name, a, b in changed:
        print(f"  {comp} | {layer} | {name}: {a} -> {b}", flush=True)
    flips = [c for c in changed if c[2] == "outside_band"]
    print("outside-band classifications changed:", len(flips), flush=True)
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tables", type=Path)
    parser.add_argument("--baseline", type=Path, default=None,
                        help="directory of the released tables, for the change report")
    args = parser.parse_args()
    sys.exit(main(args.tables, args.baseline))
