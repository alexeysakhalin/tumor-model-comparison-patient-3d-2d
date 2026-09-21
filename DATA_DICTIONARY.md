# Table schema conventions

Tables use UTF-8 CSV or TSV. The compressed background-score file contains ordinary CSV. CSV and TSV do not contain workbook authors, hidden worksheets, comments or custom Office properties. Empty fields denote unavailable or inapplicable values and must not be replaced by zero.

The original column identifiers are retained for compatibility. The principal identifiers are: `группа` = analysis group; `ген` = gene symbol; `слой` = expression layer; `орган` = organ; `образец` = sample index within layer and organ; `идентификатор` = original donor, study-qualified patient or model identifier; `значение` = score; `генов` = number of group genes; `модель` = screen arm.

Expression layers are `здоровая` (healthy), `пациент` (malignant), `органоид` (3D dome) and `2D` (2D adherent). Organ identifiers are `желудок` (stomach), `кишка` (large intestine), `лёгкое` (lung), `молочная` (breast), `панкреас` (pancreas) and `яичник` (ovary). Prostate and skin occur in source-inventory tables but are excluded from Figure 6.

In the interval table, `_lo` and `_hi` identify the lower and upper conditional bootstrap limits; `_полоса_lo` and `_полоса_hi` identify the reference-band limits. `пац−здор` is tumor minus healthy and `орг−2D` is 3D minus 2D. `_CI_исключает_0` records whether the interval excludes zero. `_вне_полосы` records whether the unrounded estimate lies outside its reference band. `_органов` is the count of organs with the same contrast sign as the aggregate, followed by the number of represented organs. Scores are in percentile points, not expression fold changes.

In the CRISPR table, `значение` is the group log2 enrichment statistic relative to the genomic background; `в_единицах_полосы` is the score divided by the reference-band half-width; `повторов` is the number of screen contrasts pooled; `ρ_между_повторами` records the reported screen-agreement summary; `ρ_полубиблиотек` is library-half agreement. The HeLa replicate-correlation field is empty because its halves do not establish biological replication. `вес` in the summary table is the median absolute normalized score across the seven arms, and `плеч_вне_полосы` is the number of arms outside their reference bands.

In the panel-B table, `RNAi_Z` is the published siRNA rescue Z score and `RNAi_выживание` is relative viability. `надёжный_хит_статьи` denotes the higher-confidence published candidates. `в_нашем_скрине` denotes availability in the reanalysed public CRISPR screen, not a new screen generated for this manuscript. `персентиль_WT` and the corresponding knockout columns are CRISPR enrichment percentiles on a 0-100 scale. `символ_в_скрине` records the matched screen symbol.

In the Patel overlap table, `генов_в_группе` is the group size, `в_опубликованном_списке` is the overlap count, `гены` lists overlapping symbols and `доля` is the overlapping fraction. The ranked-list table preserves source RIGER statistics. A source empirical P value of zero reflects the deposited calculation and is not an exact probability of zero.

In the original bootstrap log, group identifiers are implicit in row order: alphabetical group order with healthy, malignant, 3D, 2D, tumor-minus-healthy and 3D-minus-2D statistics within each group. `BOOTSTRAP_STATISTICS.tsv` makes that mapping explicit and records independently recomputed values. The numerical content of the original log is preserved.

Auxiliary model distances are differences between aggregate scores. They are not errors measured against a clinically validated predictor. Historical broader-panel coverage fields, including `in_tcga_matrix`, are annotation records and are not Figure 6 sample provenance.

The analysis scripts read these tables through code/table_io.py and schema/table_aliases.json. This boundary translates labels into English and checks for column-name collisions. The source files and their numerical values remain unchanged.
