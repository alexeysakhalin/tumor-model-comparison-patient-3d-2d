# Primary-data processing

The default workflow recalculates expression summaries from deposited profile-level scores. The primary-data utilities require separately obtained input files and have not been rerun on the full original matrices in this release. Tests verify gene-axis alignment, ranking and selection behavior on small numerical fixtures.

## Healthy epithelial profiles

`code/pseudobulk_stream.py` streams a CSR-encoded H5AD matrix and writes counts, gene symbols and aligned group metadata. Select `donor_id` as the grouping column and `compartment` as the cell-class column, retain `Epithelium`, and restrict the assay to the recorded 10x 3-prime v3 value. Use exactly the same filters in `pseudobulk_groups.py` if metadata are rebuilt separately. The count script uses `raw/X`; an X-only input requires `--matrix-is-counts` and integer, nonnegative values. Do not apply count summation to log-normalized expression.

The prepared Tabula directory requires `manifest.csv` with `organ` and `prefix` columns. Organ identifiers are `stomach`, `large_intestine`, `lung`, `breast`, `pancreas` and `ovary`. Each prefix must have `_counts.npy` (genes by groups), `_genes.csv` (gene symbols) and `_groups.tsv` (the same ordered groups). No donor identifiers are inferred when metadata are missing.

## Malignant epithelial profiles

`code/kang_pseudobulk.py` takes an input directory of atlas H5AD files and an output directory. It averages the stored log-normalized expression within sample and cell class. Epithelial cells are classified using the source `cnv_status` field. The retained epithelial groups must have at least 50 cells. All selected source studies must have their matrix, genes and group metadata available; failed files terminate the run with a nonzero status.

## Culture models

The recorded 2D inputs are `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` and `Model.csv`. The expression schema must provide `ModelID` and `is_default_entry`. Exactly one default profile per retained model is required. The 3D archive is `DepMap_NextGen_3D_core_files.zip`, containing `next_gen_expression.csv` and `model_metadata.csv`. Only models with the recorded `Adherent` or `Dome` growth pattern and organ-specific OncoTree codes are selected. Non-cancerous models are excluded.

Obtain the NextGen inputs from the DepMap file set **NextGen Model Manuscript 2026**. The associated Figshare deposition contains other analytical results and must not be substituted for these matrices. See [the authors' data instructions](https://github.com/broadinstitute/DepMap-NextGen-Public) and `SOURCE_DATASETS.tsv`.

## Configuration and score calculation

`code/build_layer_scores.py` reads an ordinary `KEY=value` configuration file; it does not execute the file as shell code. The required keys are `PSEUDOBULK_DIR` (with `tabula` and `kang` subdirectories), `KANG_STUDIES_CSV`, `DEPMAP_DIR`, `MEMBERSHIP_TSV` and `PANEL_TSV`. The last three table inputs correspond to the supplied provenance, membership and frozen candidate-panel tables. Give the configuration file as the positional argument and a separate output directory with `--output-dir`.

The main calculation uses `--rank-scope universe`. The alternative `layer` scope is a sensitivity analysis and must not be substituted for the main figure. The code checks the recorded cohort sizes and 18,612-gene intersection. These checks identify incompatible input versions; they do not replace raw-file checksum verification.

## Functional screens

The supplied archive does not contain a complete script and input set for rebuilding all seven functional arms from raw guide counts. The default workflow checks the deposited 84 arm-by-group summaries, their normalization, 12 magnitude summaries, published-list overlaps and the available MLKL/PGAM5 sensitivity inputs. It does not establish a new raw-count reconstruction of those screens. Raw sources and the analytical procedure are specified in `SOURCE_DATASETS.tsv` and `METHODS.md`.

## Provenance limits

The exact quarterly release of the traditional DepMap expression input remains unresolved. Recorded local modification dates and a separately cited essential-gene list cannot establish that release. The original ontology query underlying the broader candidate panel is also unavailable; the deposited membership is a frozen, inspectable input.

The control-selection utility preserves the supplied rule for the 200 background genes, which can include the background target itself among its controls. This differs from analyzed genes, all of which are explicitly excluded from control pools. A sensitivity analysis excluding each background target from its own controls requires the original expression matrices; it has not been reconstructed from the deposited group scores.
