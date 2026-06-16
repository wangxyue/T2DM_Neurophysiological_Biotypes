# T2DM_Neurophysiological_Biotypes
## Overview

## Code
**1. MIND_network**
- run_MIND.py: Calculates the Morphometric INverse Divergence (MIND) network. 
- The MIND networks were constructed using an open-source code (https://github.com/isebenius/MIND).

**2. Deviation scores**
- The normative model of the morphometric brain networks used in this study was derived from our previous work (Liang et al., 2025).
- Regional deviation scores relative to the established normative model were calculated using the out-of-sample estimation framework (Bethlehem et al., 2022).

**3. Biotyping**
- S1_T2DM_Biotyping_NbClust.r: Applies K-means clustering to the deviation score matrix, utilizing the NbClust package to determine the optimal number of T2DM biotypes via comprehensive index voting.
- S2_Cluster_stability_analysis.m: Performs a resampling-based stability analysis across k=2 to 10 (1,000 iterations) to evaluate the robustness of the clustering solutions.
- S3_plot_brain_deviation_map.m: Plots the biotype-specific brain deviation maps.

**4. Case_control_analysis**
- Deviations
  - S1_global_analysis.m: Conducts case-control comparisons of global deviation scores among bioype 1, biotype 2, and healthy controls (HC).
  - S2_Economo7_analysis.m: Evaluates group-level differences in deviation scores across the seven von Economo cytoarchitectonic classes.
  - S3_ROI_analysis.m: Performs ROI-level statistical analyses to identify localized brain deviation differences among the three groups.
- Cognition
  - Group_differences_in_cognition.m: Assesses behavioral phenotype differences among the three groups across multiple cognitive domains, controlling for key demographic and scanning site covariates.
- Clinical
  - Group_differences_in_clinical_10metrics.m: Assesses group-level differences in 10 shared clinical and metabolic metrics across all three groups using covariate-adjusted linear models.
  - Biotype_differences_in_clinical_3metrics.m: Evaluates biotype-specific differences in patient-only clinical metrics (Duration, HOMA2B, HOMA2IR) using covariate-adjusted linear models.

**5. Phenotype_analysis**
- Brain_deviation_cognition_PLSC.m: Performs Partial Least Squares Correlation (PLSC) analysis to identify multivariate coupling patterns between regional brain deviations and cognitive metrics within specific T2DM biotypes. The framework evaluates the statistical significance of latent variables via permutation testing and calculates Bootstrap Ratios (BSR) through resampling to identify stable brain-behavior associations.
  
- Brain_deviation_clinical_PLSC.m: Utilizes PLSC to map multivariate associations between brain deviation patterns and clinical/metabolic profiles, incorporating permutation tests for overall significance and bootstrapping for feature reliability.

**6. Gene_analysis**
- The Allen Human Brain Atlas (AHBA) datasets were preprocessed using the abagen toolbox (https://github.com/rmarkello/abagen) 
- Brain_deviation_gene_PLSR.m: Performs Partial Least Squares Regression (PLSR) to associate biotype-specific spatial patterns of brain deviation with transcriptomic profiles. The script rigorously validates spatial significance using spin tests and assesses individual gene weight stability via bootstrap resampling, outputting Z-scored ranked gene lists.
- The clusterProfiler R package was utilized to perform Gene Set Enrichment Analysis (GSEA) on the ranked gene lists.

**Spintest**
- Generate_SpinTest_Permutations.m: Generates spatially constrained null models by performing spherical rotations of the cortical surface parcellation. 
- Map_SpinTest_To_Parcels.m: Maps the vertex-level spatial permutations into ROI space to construct the empirical null distribution of brain deviation scores.

## Data

## References
