# T2DM_Neurophysiological_Biotypes
## Overview


## Code
**1. MIND_network**
- run_MIND.py: Calculates the Morphometric INverse Divergence (MIND) network.
  
- The MIND networks were constructed using an open-source code (https://github.com/isebenius/MIND) [1].

**2. Deviation scores**
- The normative model of the morphometric brain networks used in this study was derived from our previous work [2].
  
- Regional deviation scores relative to the established normative model were calculated using the out-of-sample estimation framework (https://github.com/brainchart/Lifespan) [3].

**3. Biotyping**
- S1_T2DM_Biotyping_NbClust.r: Applies K-means clustering to the deviation score matrix, utilizing the NbClust package (https://www.rdocumentation.org/packages/NbClust/versions/3.0.1/topics/NbClust) [4] to determine the optimal number of T2DM biotypes via comprehensive index voting.
  
- S2_Cluster_stability_analysis.m: Performs a resampling-based stability analysis across k=2 to 10 (1,000 iterations) to evaluate the robustness of the clustering solutions.
  
- S3_plot_brain_deviation_map.m: Plots the biotype-specific brain deviation maps.

**4. Case_control_analysis**
- Deviations
  - S1_global_analysis.m: Conducts case-control comparisons of global deviation scores among bioype 1, biotype 2, and healthy controls (HC).
    
  - S2_Economo7_analysis.m: Evaluates group-level differences in deviation scores across the seven von Economo cytoarchitectonic classes [5].
    
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
- The Allen Human Brain Atlas (AHBA) datasets (https://human.brain-map.org/static/download) [6] were preprocessed using the abagen toolbox (https://github.com/rmarkello/abagen) [7].
  
- Brain_deviation_gene_PLSR.m: Performs Partial Least Squares Regression (PLSR) to associate biotype-specific spatial patterns of brain deviation with transcriptomic profiles. The script rigorously validates spatial significance using spin tests and assesses individual gene weight stability via bootstrap resampling, outputting Z-scored ranked gene lists.
  
- The clusterProfiler R package (https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) [8] was utilized to perform gene set enrichment analysis [9] on the ranked gene lists.
  
- The canonical cell classes data (https://github.com/jms290/PolySyn_MSNs/tree/master/Data) [10] was used to perform cell-type enrichment analysis.

**Spintest**

Spin-based spatial permutation test [11]:
- Generate_SpinTest_Permutations.m: Generates spatially constrained null models by performing spherical rotations of the cortical surface parcellation.
  
- Map_SpinTest_To_Parcels.m: Maps the vertex-level spatial permutations into ROI space to construct the empirical null distribution of brain deviation scores.

## Data

1. The spatial distribution maps of the mean brain deviation scores for T2DM biotypes: Biotype1_mean_deviation.mat, Biotype2_mean_deviation.mat
   
2. The mapping index assigning each ROI to one of the seven von Economo cytoarchitectonic classes: parcel_to_Economo.mat
   
3. The cortical sensorimotor-association (S-A) axis mapped onto the custom DK318 atlas: SAaxis_318_grad.mat

## References
[1] Sebenius, Isaac, et al. "Robust estimation of cortical similarity networks from brain MRI." Nature neuroscience 26.8 (2023): 1461-1471.

[2] Liang, Xinyuan, et al. "Dissecting human cortical similarity networks across the lifespan." Neuron 113.19 (2025): 3275-3295.

[3] Bethlehem, Richard AI, et al. "Brain charts for the human lifespan." Nature 604.7906 (2022): 525-533.

[4] Charrad, Malika, et al. "NbClust: an R package for determining the relevant number of clusters in a data set." Journal of statistical software 61 (2014): 1-36.

[5] von Economo, Constantin Freiherr, and Georg N. Koskinas. Die cytoarchitektonik der hirnrinde des erwachsenen menschen. J. Springer, 1925.

[6] Hawrylycz, Michael J., et al. "An anatomically comprehensive atlas of the adult human brain transcriptome." Nature 489.7416 (2012): 391-399.

[7] Markello, Ross D., et al. "Standardizing workflows in imaging transcriptomics with the abagen toolbox." elife 10 (2021): e72129.

[8] Yu, Guangchuang, et al. "clusterProfiler: an R package for comparing biological themes among gene clusters." Omics: a journal of integrative biology 16.5 (2012): 284-287.

[9] Subramanian, Aravind, et al. "Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles." Proceedings of the national academy of sciences 102.43 (2005): 15545-15550.

[10] Seidlitz, Jakob, et al. "Transcriptomic and cellular decoding of regional brain vulnerability to neurogenetic disorders." Nature communications 11.1 (2020): 3358.

[11] Alexander-Bloch, Aaron F., et al. "On testing for spatial correspondence between maps of human brain structure and function." Neuroimage 178 (2018): 540-551.
