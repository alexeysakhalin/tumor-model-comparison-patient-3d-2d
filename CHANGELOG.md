# Change record

## 1.1.0 — 21 September 2026

**Every scored gene is excluded from its own control pool.** Background genes could previously select
themselves at distance zero, which affected the reference bands and everything normalised by them.
Analysed estimates and intervals are unaffected and were verified unchanged.

The released baseline was reproduced first: all 160,800 background values matched
one-to-one within 1.4e-14 percentile points. Self-inclusion had
occurred in 1,199 of 1,200 gene-by-organ selections.
5 of 24 normalised contrasts change, each by at most
0.04 reference-band units; all 48 bands widen; no classification in
Figure 6 changes; two quoted band limits move.

Changed: `code/level_matched_scores.py`, `code/build_summary_tables.py` (rounding made an argument,
default unchanged), `code/recompute_pgam5_sensitivity.py` (band and observation helpers factored out,
results unchanged for identical inputs), `code/test_numerics.py`, the eleven regenerated tables, two
`source_data/` exports, `figures/Figure6.*`, `validation/manuscript_values.json`,
`FIGURE6_PANEL_MAP.tsv`, `ZENODO_TABLE_SNAPSHOT.json`, `MANIFEST_SHA256.tsv`, `METHODS.md`,
`PGAM5_SENSITIVITY_METHODS.md`, `docs/PRIMARY_INPUTS.md`.
New: `code/audit_self_exclusion.py`, `code/report_self_exclusion.py`,
`code/apply_self_exclusion_correction.py`, `audit/self_exclusion_2026-09-21/`.

## 1.0.0 — 19 September 2026

First public release: the Figure 6 analysis, its derived tables, source-data exports, provenance
records, reference metadata and verification workflow.

### Export consistency

The composite PDF, PNG and three standalone panels are regenerated in the recorded rendering environment. The exact pixel comparison remains part of the verification workflow. PDF creation timestamps are omitted to make repeated exports deterministic. Ordered-control digests are labeled explicitly as 16-character SHA-256 prefixes.
