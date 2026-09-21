# Figure 6 methods

## Expression cohorts and anatomical matching

The analysis compares healthy epithelial profiles from Tabula Sapiens 2.0, malignant epithelial profiles from the atlas of Kang et al., 3D dome cultures from the DepMap NextGen collection and traditional 2D adherent DepMap models. All four groups are represented in stomach, large intestine, lung, breast, pancreas and ovary. The retained cohorts contain 20 healthy profiles from 11 donors, 231 malignant profiles from 202 study-qualified patients in 34 studies, 152 dome models and 401 adherent models. The patient provenance table contains 36 organ-by-study records because some studies contribute to more than one anatomical site.

Culture models are restricted by growth pattern and anatomical origin. One default expression profile is retained per model. Model identifiers, OncoTree codes and growth patterns are listed in tables/PROV_model_sites.csv. Cohorts are matched by anatomical site; they are not paired patient tumors and derived cultures. Culture format and cohort provenance are therefore not separable as causal effects in this comparison.

Healthy epithelial pseudobulk profiles use summed counts from cells annotated as Epithelium and measured with the 10x 3-prime v3 assay. Malignant epithelial profiles use the mean of the stored log-normalized values within each sample, with at least 50 malignant epithelial cells per retained sample. The deposited model expression values are used for the two culture groups. These different modalities are compared through within-profile ranks rather than directly on an absolute expression scale.

## Gene sets

The plotted collection contains 67 genes assigned to 12 non-overlapping computational groups. Forty-one genes derive from the broader 108-gene candidate panel, and 26 additional genes describe antigen presentation and immune interactions. Full membership and display labels are in tables/IKM_gene_membership_by_group.tsv; GROUP_ANNOTATIONS.tsv documents functional scope and references. The organizational domains do not imply biological independence or sequential barriers to immune killing.

MLKL/PGAM5 is a descriptive pair. PGAM5 is not an obligatory terminal necroptotic effector, and the combined score is not interpreted as a canonical necroptosis-execution module. The added single-gene sensitivity analysis is described in PGAM5_SENSITIVITY_METHODS.md. MLKL alone remains inside the corresponding reference band in all four layers. The interpretation remains descriptive even after this check. Death-receptor mediators comprise the adaptors FADD and TRADD, the initiator caspases CASP8 and CASP10, and BID. NK receptor ligands and immune-regulatory ligands include genes with receptor-dependent effects. IFNGR1, IFNGR2, JAK1, JAK2, STAT1 and IRF1 are not a separate expression group; part of this axis is assessed in the functional comparisons. The groups are selected mechanistic sets, not exhaustive pathways or validated resistance signatures.

## Expression scores and aggregation

Within-profile percentile ranks are calculated over the common universe of 18,612 gene symbols. Tied expression values receive their average rank. For each target gene and anatomical site, 400 control genes are selected by proximity of their median rank in the 2D reference models. The analyzed genes and broader candidate panel are excluded from the control pool. A gene score is its rank in a profile minus the median rank of its matched controls in that profile. In the supplied primary-processing implementation, the control pool excludes the analyzed genes and the broader candidate panel. A randomly selected background gene is not additionally excluded from its own control pool. The deposited reference-gene scores preserve this rule. The size of any change after excluding such self-matches cannot be assessed without the primary rank matrices and control assignments.

The plotted group estimator first takes the median across genes within a profile, then the median across profiles within an anatomical site, and finally the median across the six sites. For a contrast, the difference between the two site-specific group scores is calculated before taking the median across sites. The 2D reference is approximately centered on zero by construction; individual genes need not have a score of exactly zero.

tables/IKM_final_master_table_CI.csv contains the plotted estimator. The auxiliary gene-first estimator in tables/IKM_final_master_table.csv takes the median across profiles for each gene before aggregating genes within a site. The estimators differ because medians are not associative. The gene-first table is retained as an auxiliary analysis and is not substituted for the plotted table.

## Cluster-bootstrap intervals

The 95% intervals use the 2.5th and 97.5th percentiles of 2,000 accepted cluster-bootstrap replicates. Clusters are donors for healthy epithelium, studies for malignant epithelium and model identifiers for the culture groups. A cluster is sampled once per layer and its multiplicity is applied jointly across all sites to which it contributes. Gene membership and the six-site target are fixed.

Draws that leave any of the six sites without observations are rejected and resampled. The intervals are therefore conditional on representation of all six sites. Rejection counts and fractions are reported in tables/IKM_bootstrap_rejection.csv. The gastric malignant stratum is supplied by one study and cannot estimate between-study variability within that stratum. An interval containing zero does not establish equivalence. Filled expression markers in panel c denote intervals excluding zero; hollow markers denote intervals containing zero.

## Descriptive expression reference bands

Reference bands are the central 95% ranges obtained from 2,000 size-matched random gene sets drawn from a fixed pool of 200 background genes. The same aggregation is applied to these sets within each layer or contrast. The bands are descriptive comparisons with the selected background pool. They do not constitute multiplicity-adjusted significance tests, and they are not a calibrated null model for every possible co-regulated pathway.

Panel a divides each expression contrast by half the width of its own reference band. A band need not be symmetric around zero, so a normalized magnitude above one is not equivalent to falling outside the band. Frame status is computed from the original band limits. The bootstrap interval and the reference band answer different questions and determine different graphical elements.

## Functional CRISPR comparisons

Seven selection arms from five models are summarized. Six derive from Watterson et al.: SK-MEL-2, A375, HT-29, HT-29 JAK1 knockout and CRC-9 under IFN-gamma, and CRC-9 under autologous T-cell selection. One additional arm derives from the HeLa TNF-alpha plus IFN-gamma screen of Zhou et al. The arms aggregate 19 screen contrasts: 18 Watterson contrasts and one pooled HeLa screen. The two HeLa guide-library halves are combined as parts of one library and are not biological replicates.

For Watterson data, the supplied processing record specifies counts-per-million normalization, log2 transformation, exclusion of guides with fewer than 30 baseline counts, treatment-minus-matched-control contrasts and centering on non-targeting guides. Guide effects are averaged within genes and screen effects within arms. For HeLa, guide log2 fold changes are calculated from normalized treated and control counts, each library half is median-centered and guide effects are averaged across both halves within each gene. The deposited derived tables contain the resulting summaries; full raw-count reconstruction is outside the scope of this package.

Within each arm, a group statistic is the median effect of its genes minus the genomic median. It is divided by half the width of the arm-specific 95% reference band from 2,000 size-matched random gene sets. Positive values indicate relative enrichment of knockouts under selection, consistent with resistance; negative values indicate depletion, consistent with sensitization or reduced relative fitness. Recognition, proliferation and general fitness can also contribute to these effects.

Frames identify values outside the reference band. Solid frames denote a consistent sign across the contributing screens; dotted frames denote disagreement. The corresponding expression frames use consistency across anatomical sites. These indicators are descriptive and are not adjusted for the number of cells in the matrix. Replicate correlations are displayed only where replicate screens exist. HeLa has no biological replicate estimate; the correlation between its library halves is recorded separately. Panel c shows the median across seven arms of the absolute normalized group statistic, which conveys magnitude without direction.

## Cross-platform comparison and published-list overlap

Panel b compares the 34 published protective candidates from the HT-29 kinome siRNA screen of Woznicki et al. with their percentiles in the HT-29 CRISPR ranking from Watterson et al. Historical symbols are reconciled through Entrez Gene identifiers. JAK1 and JAK2 meet the original siRNA study's high-confidence criterion and are highlighted. The siRNA experiment used TNF-alpha plus IFN-gamma, whereas the CRISPR arm used IFN-gamma alone. The comparison does not isolate platform effects or establish a common death mechanism.

The final column of panel a reports direct overlap with the published Mel624 TCR-selection candidate list of Patel et al. The supplied list contains 554 identifiers, including 502 non-miRNA identifiers used in this comparison. The displayed value is the number of group genes in that list divided by the group size. A suitable eligible-screen background was not reconstructed, so no new enrichment ratio, hypergeometric test or false-discovery rate is reported. Absence from the published list is not evidence of no functional effect. TCR-mediated recognition is distinct from direct CD19 CAR recognition.

## Sources and reproducibility boundaries

SOURCE_DATASETS.tsv specifies the eight primary resources, publication DOIs and relevant data components. REFERENCE_METADATA.json contains bibliographic metadata. The recorded input filenames, sizes and hashes remain in the provenance tables. The traditional DepMap release identifier is absent from the input records and cannot be inferred from a general CCLE citation. The exact source files for the NextGen expression and model annotation are associated with the DepMap file set named NextGen Model Manuscript 2026 by the publication.

All six plotted source tables match their corresponding analysis tables. The original 72 expression estimates and conditional intervals, 48 expression reference bands, 12 functional magnitude summaries and 12 candidate-list overlap counts are unchanged. The deposited PDF and PNG renderings correspond to these six numerical source tables. The code directory contains the plotting script, expression-summary reconstruction and verification workflow. Primary-processing utility inputs and limits are described in docs/PRIMARY_INPUTS.md. SOFTWARE_VERSIONS.json records the reconstruction environment. Complete primary-expression reconstruction and original ontology-query reconstruction are not established by these checks.

The additional DepMap provenance tables retain reported local timestamps, file sizes, checksums and comparisons of Model.csv copies. Local modification dates do not establish when the portal changed or which quarterly release supplied a file. The comparison of four source downloads was reported by the supplying analyst; the original copies are not included. Truncated SHA-256 values are explicitly labeled as prefixes. No exact traditional DepMap release is inferred.

## MLKL and PGAM5 sensitivity

The eight single-gene expression estimates and one-gene Monte Carlo bands can be reconstructed from the deposited gene-level and random-gene tables. The pair estimates in the auxiliary sensitivity table were corrected from the gene-first estimator to the sample-first estimator used in Figure 6. The reproducible implementation and a comparison with exact enumeration of all 200 background genes are in PGAM5_SENSITIVITY_METHODS.md. Removing PGAM5 places MLKL inside the one-gene band in all four layers under both methods. PGAM5 is outside the band in malignant epithelium under both methods; its healthy-layer classification changes between Monte Carlo sampling and exact enumeration. This distinction is preserved in the full-precision sensitivity table. The main Figure 6 values and bands are unchanged.

The auxiliary HT-29 comparison reconciles with the supplied per-gene screen values within their rounding precision. Exclusion of PGAM5 changes the group statistic by approximately 0.03 to 0.45 log2 units across the 12 configuration-condition combinations. Neither these changes nor the expression-band comparison demonstrate a death modality. Individual-gene cluster-bootstrap intervals are not estimated because the required profile-level target-gene scores are not deposited.

## Data sources and software

The exact expression and screen source components are listed in SOURCE_DATASETS.tsv, with full bibliographic details in REFERENCE_METADATA.json. NextGen expression and model metadata are from the DepMap file set NextGen Model Manuscript 2026 associated with Neiswender et al. The linked Figshare collection contains other analytical results. The shared-gene scoring implementation uses seeds 11 for the background-gene pool, 3 for random-set reference bands and 7 for cluster resampling. Recalculation preserves the deposited background-gene order.
