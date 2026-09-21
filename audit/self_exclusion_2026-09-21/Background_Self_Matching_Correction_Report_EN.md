# Correction of background self-matching in Figure 6

Release 1.1.0, 21 September 2026. Baseline: commit `9702153c7564c2fa48bf2d58f052066c0e4bb442` of this
repository, which remains the reference for every comparison below.

## What was wrong

A scored gene's controls are the 400 genes nearest to it in reference rank within an anatomical site.
The candidate pool was the shared gene universe minus the frozen exclusion set, which holds the
analysed genes and the broader candidate panel. `group_table` adds the analysed genes to that set, so
they were never their own controls. The 200 background genes are scored by a separate call that passes
the same set, which does not contain them; matching is by absolute distance, so each background gene
sat at distance zero from itself and was taken as its own control in 1,199 of the 1,200 gene-by-organ
selections. The background scores define the
reference bands, and the bands are the denominator of the normalised contrasts and the criterion of the
outside-band indicators, so the defect reached Figure 6 through the bands and through nothing else.

## What changed in the code

`control_sets` in `code/level_matched_scores.py` drops the scored gene from the ordered candidates
before the nearest 400 are taken, so the set holds 400 distinct controls, and raises when fewer than
400 eligible genes remain. `controls_are_nearest` validates against the same eligible pool.
`code/test_numerics.py` fixes the rule on explicit fixtures: a gene outside the exclusion set is not
its own control, the property holds inside a full tie block, and an insufficient pool raises. Removing
the exclusion line fails three of those tests. Nothing else in the selection was altered: the pool is
still the universe minus the frozen exclusions, the background pool is not excluded as a whole, ties
are still broken by gene symbol, and the count, seeds and draw sequences are unchanged.

## Reproduction of the released baseline

`code/audit_self_exclusion.py` rebuilds the four rank matrices from the primary inputs
(18,612 shared symbols; 20 healthy,
231 malignant, 152 dome and
401 2D profiles) and scores the same 200 background genes under both rules in one
process. All 160,800 released values were matched one-to-one on layer, organ, gene
and profile index — the biological cluster identifier is not unique within a layer, so it is verified
separately rather than used as the key — with a maximum absolute difference of
1.421e-14 percentile points against a declared tolerance of 1e-9.
The script fails if rows are missing, if keys are not unique, or if that tolerance is exceeded.

## How often self-inclusion occurred

1,199 of 1,200 gene-by-organ selections included the
scored gene. The exception is SMIM40 in stomach: its reference rank lies inside a tie
block where 400 other genes at distance zero precede it by symbol. After the correction no selection
contains its own gene, and every corrected set holds 400 distinct controls.

## Measured effect

| Quantity | Result |
|---|---|
| Background scores (160,800 gene-by-profile) | 59.7 per cent unchanged; 99th percentile of \|change\| 0.19611; largest 11.98555 for ZNF222 in one large_intestine profile of the malignant layer |
| Direction | 32,981 increase, 31,756 decrease, median 0.00 in every layer |
| Reference bands (48) | every band wider, by 0.000537 to 0.263899; edges move by a median of 0.026734 and at most 0.264554 |
| Normalised contrasts (24) | 5 change, each by at most 0.04 reference-band units; 0 classification changes |
| Group-by-layer classifications (72 comparisons: 48 layer rows and 24 contrasts) | 0 changes; 0 score changes |
| MLKL/PGAM5 sensitivity (12) | 0 Monte Carlo and 0 exact-enumeration classification changes; the one method disagreement (PGAM5 in healthy epithelium) existed before the correction and remains |
| Manuscript values (22) | 2 band limits move; no point estimate or interval moves |

The direction of the band change is the one the mechanism predicts: a control at distance zero carries
the target's own rank, which pulled the control median toward the target and made the background
statistics slightly too tightly distributed.

## Convention of the comparison

Band limits are stored to one decimal place and a normalised contrast is the full-precision statistic
divided by half the width of the stored band, rounded to two decimals — the convention of the released
tables, kept unchanged so that this comparison isolates self-exclusion. Band-edge changes above are
quoted at full precision, obtained from the released routine with its rounding switched off
(`null_bands(..., decimals=None)`), which does not affect what the tables carry.

## What was regenerated

`tables/IKM_random_gene_scores_by_sample.csv.gz`, `tables/IKM_null_bands.csv`,
`tables/IKM_final_master_table.csv`, `tables/IKM_final_master_table_CI.csv`,
`tables/IKM_representation_in_band_units.csv`, `tables/IKM_models_vs_patient.csv`,
`tables/IKM_bootstrap_rejection.csv`, `tables/IKM_single_cluster_strata.csv`,
`tables/IKM_pgam5_sensitivity.csv`, `tables/IKM_pgam5_sensitivity_recomputed.csv`,
`tables/IKM_pgam5_reference_gene_scores.csv`, `tables/PGAM5_RECONSTRUCTION_CHECKS.json`, the two
affected `source_data/` exports, `figures/Figure6.png`, `figures/Figure6.pdf`,
`validation/manuscript_values.json`, `FIGURE6_PANEL_MAP.tsv`, `ZENODO_TABLE_SNAPSHOT.json` and
`MANIFEST_SHA256.tsv`. Functional-screen inputs, the published-list overlaps, the 72 group estimates
and their conditional intervals are unchanged.

## Limits

The primary expression matrices and the assembled rank matrices are not redistributed here; the audit
identifies them by path, size and checksum in `input_manifest.tsv`, and `docs/PRIMARY_INPUTS.md`
describes how they are prepared. Control-set identity for all 67 analysed genes under both rules is recorded in
`analysed_gene_control_sets.tsv`: the ordered lists are identical in every one of the 402 gene-by-site selections, which is
what makes the analysed estimates independent of this correction; the ordered control sets themselves are
deposited for the 200 background genes only, since those are the genes the correction moves. The quarterly release of the traditional DepMap expression input remains unresolved and is not
inferred from a checksum. No manuscript file is part of this repository: the values a manuscript would
quote are in `validation/manuscript_values.json` and the two that moved are listed in
`manuscript_value_changes.tsv`, but the manuscript text itself has not been checked here.

## The malignant selection input and the deposited provenance table

The loader reads the study-selection table named in the configuration: 2,269 bytes, seven columns, SHA-256 prefix `1803b234cab3`, recorded in `input_manifest.tsv` under `malignant_selection`. The repository ships an expanded version of the same table, `tables/PROV_kang_primary_studies.csv` (4,918 bytes, nine columns, prefix `a6c15a0ec346`), which adds `source_reference` and `source_url`. The two are not byte-identical and are not presented as such. They carry the same 36 study-by-organ rows and the same 34 datasets, and six of the seven shared columns agree row for row, including the dataset-to-organ mapping, which with the `Dataset` column is all the loader reads from this file. The seventh column, the accession, differs in six rows: where the selection input records that the accession is not recoverable from the file name, the deposited table carries the study citation instead (`br_pan_blueprint`, `crc_pan_blueprint`, `lu_nsclc_pan_blueprint_lung`, `ov_pan_blueprint_ovary` — Qian et al., Cell Res 2020; `lu_adc_codeocean` — Bischoff et al., Oncogene 2021; `ov_hgstoc_pan_blueprint_ovary` — Olbrecht et al., Genome Med 2021). The deposited table is therefore the same selection with citations attached and six accession fields filled in from the literature rather than from a file name. The selection input is the smaller file; the deposited table is the same selection with citations attached.

## Files

| File | Content |
|---|---|
| `input_manifest.tsv` | every file the loaders open, one row each: 127 files in 13 roles, with a portable logical path, size and SHA-256. The enumeration follows the selection rules of the calculation — the healthy layer through `pseudobulk/tabula/manifest.csv`, the malignant layer through the study-selection table, each dataset contributing its matrix, gene axis and group table — and fails if a required component selects nothing. |
| `analysed_gene_control_sets.tsv` | the 67 analysed genes in all six sites (402 selections): control set identical before and after, self-inclusion count, 16-character SHA-256 prefixes of the ordered control-index lists |
| `baseline_reconstruction_report.json` | cohort sizes, reconciliation, tolerance, self-inclusion counts |
| `self_exclusion_audit.tsv` | 1,200 rows: gene by organ, self-inclusion before and after, removed and added control, rank gaps |
| `control_assignments_baseline.tsv.gz`, `control_assignments_corrected.tsv.gz` | the ordered control sets under both rules |
| `background_score_changes.tsv.gz` | 160,800 rows: released score, corrected score, difference |
| `reference_band_changes.tsv` | 48 bands at full precision, both rules, widths and changes |
| `representation_changes.tsv` | the 24 normalised contrasts, released and corrected |
| `layer_classification_changes.tsv` | group-by-layer estimates and outside-band flags, both rules |
| `sensitivity_changes.tsv` | the 12 MLKL/PGAM5 rows in both reference schemes |
| `manuscript_value_changes.tsv` | the 22 validated manuscript values, released and corrected |
| `correction_effect_summary.json` | the measured effects quoted above |

Reproduce with `python code/audit_self_exclusion.py <pipeline.conf> --tables <released tables>
--output-dir DIR`, then `python code/report_self_exclusion.py --released <released tables> --corrected
tables --output-dir DIR`. The selection-rule tests in `code/test_numerics.py` need no primary data.
