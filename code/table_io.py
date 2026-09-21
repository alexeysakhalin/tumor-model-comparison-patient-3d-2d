"""Read deposited tables through the documented English schema.

The source files retain their original bytes to match the Zenodo data deposit.
The translation affects column names and categorical labels, never numbers.
"""

from functools import lru_cache
from pathlib import Path
import json

import pandas as pd


@lru_cache(maxsize=1)
def aliases():
    path = Path(__file__).resolve().parents[1] / "schema/table_aliases.json"
    return json.loads(path.read_text(encoding="utf-8"))


def english_label(value):
    if isinstance(value, str):
        for source, target in aliases().items():
            value = value.replace(source, target)
    return value


def read_table(path, **kwargs):
    frame = pd.read_csv(path, **kwargs)
    frame.columns = [english_label(c) for c in frame.columns]
    if frame.columns.duplicated().any():
        raise ValueError(f"Column aliases are not unique: {path}")
    for column in frame.select_dtypes(include=["object", "string"]).columns:
        frame[column] = frame[column].map(english_label)
    return frame
