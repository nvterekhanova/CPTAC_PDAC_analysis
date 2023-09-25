# Analysis for CNV-impact on protein/gene expression, that was part of Cell publication (2021).

## SCNA arm and focal significance.

  * To infer focal-level significant somatic copy number alterations (SCNA) we used GISTIC2 (Mermel et al., 2011) with the default parameters except for increased thresholds for amplifications and deletions (i.e., -ta and -td parameters of GISTIC2), that were set to 0.4, and confidence level set to 0.95. This analysis was performed on the segment-level SCNA data from the autosomes.

## Prioritizing Putative SCNA Drivers.

  * We first filtered all the genes to those with quantifiable copy number, gene expression, and proteomics (N=11,623). Next, we also filtered genes for those occurring in the focal amplified regions identified by GISTIC2 with Q value < 0.25 (N = 543). Finally, we filtered the genes by their CN-mRNA correlation and CN-protein correlation to keep the genes with significant CN cis-effect (FDR < 0.05, Spearman’s correlation). The resulted set of genes (N=23) was used for the over-representation analysis to identify significantly enriched GO-biological processes.


## Analyses:

1. ```Barplot_driver_CNVs_summary``` -- Barplot summary of top frequent CNVs in PDAC cohort.

2. ```CNV_impact_protein_gExpr``` -- CNV impact on gene expression and protein.

3. ```GISTIC``` -- visualization of GISTIC results using ```maftools```.

   * GISTIC was run using GenePattern web-server: https://www.genepattern.org/modules/docs/GISTIC_2.0.