# T2DM_Neurophysiological_Biotypes
## Overview

## Code
**1. MIND_network**
- run_MIND.py : Morphometric INverse Divergence (MIND) network calculation
- Open Python code for the estimation of MIND networks (https://github.com/isebenius/MIND)

**2. Deviation scores**
- The normative model of the morphometric brain networks used in this study was derived from our previous work (Liang et al., 2025).
- Regional deviation scores relative to the established normative model were calculated using the out-of-sample estimation framework (Bethlehem et al., 2022).

**3. Biotyping**
- S1_T2DM_Biotyping_NbClust.r : Applies K-means clustering to the deviation score matrix, utilizing the NbClust package to determine the optimal number of T2DM biotypes via comprehensive index voting.
- S2_Cluster_stability_analysis.m : Performs a resampling-based stability analysis across k=2 to 10 (1,000 iterations) to evaluate the robustness of the clustering solutions.
- S3_plot_brain_deviation_map.m: Plots the biotype-specific brain deviation maps.

**4. Case_control_analysis**
- Deviations
  - S1_global_analysis.m : Conducts case-control comparisons of global deviation scores among bioype 1, biotype 2, and healthy controls (HC).
  - S2_Economo7_analysis.m : Evaluates group-level differences in deviation scores across the seven von Economo cytoarchitectonic classes.
  - S3_ROI_analysis.m : Performs ROI-level statistical analyses to identify localized brain deviation differences among the three groups.
- Cognition
  - Group_differences_in_cognition.m：Assesses behavioral phenotype differences among the three groups across multiple cognitive domains, controlling for key demographic and scanning site covariates.
- Clinical
  - Group_differences_in_clinical_10metrics.m
  - Biotype_differences_in_clinical_3metrics.m

**5. Phenotype_analysis**
- Brain_deviation_cognition_PLSC.m
- Brain_deviation_clinical_PLSC.m

**6. Gene_analysis**
- Brain_deviation_gene_PLSR.m

**Spintest**

## Data

## References
