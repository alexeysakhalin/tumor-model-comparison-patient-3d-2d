# MLKL and PGAM5 sensitivity analysis

The Figure 6 row is a descriptive MLKL/PGAM5 pair. The new analysis compares the pair with each gene separately. It does not define a necroptosis activity score.

## Expression estimates

For MLKL and PGAM5 separately, the median across the six anatomical sites is calculated from the site-specific gene scores in tables/IKM_final_per_gene_by_organ.csv. With one gene, the gene-first and sample-first aggregation orders coincide. The pair is reconstructed from tables/IKM_group_scores_by_sample.csv using the sample-first estimator of Figure 6. The pair is not reconstructed by averaging the two final single-gene scores, because medians are not associative.

The eight single-gene estimates, band limits and outside-band indicators are reported in tables/IKM_pgam5_sensitivity.csv. The four pair estimates use the sample-first estimator of Figure 6: -15.7, -10.7, -5.3 and 0.0 percentile points for healthy epithelium, malignant epithelium, 3D and 2D models, respectively. The auxiliary gene-first estimator gives different values because medians are not associative; it is not substituted for the plotted estimator.

The single-gene reference distribution is obtained by taking the median within each site and then across the six sites for each of the 200 genes in tables/IKM_random_gene_scores_by_sample.csv.gz. The supplied one-gene bands are reproduced by sorting gene symbols, sampling 2000 indices with replacement using NumPy default_rng(3), and taking the 2.5th and 97.5th percentiles with linear interpolation. These settings reproduce the deposited Monte Carlo estimates to their reported precision.

## Exact enumeration of the reference pool

For a one-gene set, all 200 possible reference values can be enumerated without Monte Carlo sampling. The exact empirical central 95% range of those 200 values is therefore provided as an additional sensitivity check. This is exact only for the deposited finite reference pool; it does not remove uncertainty in the choice of background genes or in the biological samples. The two-gene reference bands of the main figure are unchanged and are not replaced by one-gene bands.

MLKL remains inside the one-gene band in all four layers under both methods. In malignant epithelium, PGAM5 is outside the band under both methods. In healthy epithelium, PGAM5 is inside the reconstructed Monte Carlo band but outside the exact empirical band: its score is -24.715 percentile points, the Monte Carlo lower limit is -25.621, and the exact empirical lower limit is -22.142. That healthy-layer classification is sensitive to Monte Carlo sampling of a finite pool and should not be described as a robust inside-band result. Full-precision estimates, both sets of limits and the 800 gene-layer reference values are deposited.

Excluding PGAM5 removes the pair's outside-band classification in healthy and malignant epithelium. Because both gene composition and the size-specific reference distribution change, this comparison does not establish that PGAM5 alone causally determines the pair's displacement. It also does not demonstrate reduced necroptotic competence. PGAM5 was dispensable for necroptosis in the experimental systems examined by Moriwaki et al. (2016), supporting the descriptive interpretation of this pair (doi:10.4049/jimmunol.1501662).

## HT29 screen estimates

tables/IKM_pgam5_screen_gene_values.csv contains the MLKL and PGAM5 rows extracted from the supplied per-gene HT-29 screen table. The four entries are HT-29 screen configurations, including two JAK1 knockout clones; they are not the seven aggregated arms of Figure 6. Baseline, IFNg and IL6 are retained as distinct conditions. The 24 values in tables/IKM_pgam5_sensitivity_arms.csv reconcile with the per-gene medians within 0.0005 log2 units, the precision limit imposed by the three-decimal input values. The absolute difference between the pair and MLKL alone is approximately 0.03 to 0.45 log2 units across all 12 configuration-condition combinations, or 0.04 to 0.45 for IFNg alone. No uncertainty or reference-band classification for these auxiliary screen values is estimated.

## Reproduction and boundaries

The reconstructed estimates and reference limits are in tables/IKM_pgam5_sensitivity_recomputed.csv; the 800 reference gene-layer values are in tables/IKM_pgam5_reference_gene_scores.csv. tables/PGAM5_RECONSTRUCTION_CHECKS.json records the numerical checks. Analysis scripts are maintained separately from this data archive. Project repository: https://github.com/alexeysakhalin/tumor-model-comparison-patient-3d-2d

Individual-gene cluster-bootstrap intervals were not recalculated because individual-gene scores for every biological profile are not included. Such intervals are mathematically defined for a one-gene set; their absence here is a data-availability limitation. The original Figure 6 intervals and its six plotted source tables remain unchanged.
