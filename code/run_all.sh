#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
export MPLBACKEND=Agg
python code/test_numerics.py
python code/test_gene_axis.py
python code/verify_package.py
python code/build_summary_tables.py tables --output-dir build/tables
python code/recompute_pgam5_sensitivity.py tables --check --output-dir build/sensitivity
python code/make_full_figure.py tables build/Figure6
python code/verify_package.py --build-dir build
