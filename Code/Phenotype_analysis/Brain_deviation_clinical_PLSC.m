%% PLSC: Brain Deviation - Clinical/Metabolic Multivariate Coupling Patterns

clear; clc; close all;

% Target Biotype to analyze (1 = Biotype 1; 2 = Biotype 2)
TARGET_BIOTYPE = 2;

% Expected variables in .mat: X_brain, Y_clin, Biotype_Label
CLEANED_DATA_PATH = 'path/to/data/Cleaned_Brain_Clinical_Data.mat';
OUTPUT_DIR = sprintf('path/to/results/PLSC_Output/Clinical/Biotype%d', TARGET_BIOTYPE);

% PLSC Statistical Parameters
nPerm = 10000;      
nBoot = 1000;       
rng_seed = 20260320; 
bsr_thresh = 2.0;    % Threshold for stable brain regions (|BSR| > 2)

rng(rng_seed);
if ~exist(OUTPUT_DIR, 'dir'), mkdir(OUTPUT_DIR); end

%% 1. Load Data & Define Variables
fprintf('Loading data and initializing parameters...\n');
load(CLEANED_DATA_PATH); 

% The 9 clinical variables directly entering PLSC
clinical_vars = {'BMI', 'Duration_Year', 'HbA1c', 'SBP', 'DBP', 'TC', 'TG', 'HDL', 'LDL'};

%% 2. Filter by Biotype & Handle Missing Values
fprintf('Selecting Biotype %d and removing missing data...\n', TARGET_BIOTYPE);

biotype_idx = (Biotype_Label == TARGET_BIOTYPE);
X_sub = X_brain(biotype_idx, :);
Y_sub = Y_clin(biotype_idx, :);

% Keep only complete cases
valid_idx = all(~isnan(X_sub), 2) & all(~isnan(Y_sub), 2);
X_clean = X_sub(valid_idx, :);
Y_clean = Y_sub(valid_idx, :);

[n_subj, p_brain] = size(X_clean);
q_clin = size(Y_clean, 2);
fprintf('Final sample size for PLSC: N = %d\n', n_subj);

%% 3. Standardization (Z-scoring)
fprintf('Standardizing matrices...\n');
Xz = zscore(X_clean);
Yz = zscore(Y_clean);

if any(std(Xz, 0, 1) == 0) || any(std(Yz, 0, 1) == 0)
    error('Zero variance detected in standardized matrices. Check input data.');
end

%% 4. Execute PLSC (SVD)
fprintf('Running PLSC (SVD)...\n');
plsc = run_plsc_svd(Xz, Yz);

U = plsc.U;                  
V = plsc.V;                  
singvals = plsc.singvals;    
explained = plsc.explained;  

fprintf('Extracted %d latent variables (LVs).\n', length(singvals));

%% 5. Permutation Test (Significance)
fprintf('Running Permutation Test (%d permutations)...\n', nPerm);
perm_results = permutation_test_plsc(Xz, Yz, nPerm);

perm_pvals = zeros(length(singvals), 1);
for lv = 1:length(singvals)
    perm_pvals(lv) = mean(perm_results.perm_singvals(:, lv) >= singvals(lv));
end

perm_explained = (perm_results.perm_singvals.^2) ./ sum(perm_results.perm_singvals.^2, 2);

%% 6. Bootstrap Resampling (Robustness/Reliability)
fprintf('Running Bootstrap Resampling (%d iterations)...\n', nBoot);
boot_results = bootstrap_plsc(Xz, Yz, U, V, nBoot);

brain_bsr = U ./ boot_results.brain_se;
clin_bsr  = V ./ boot_results.clin_se;

%% 7. Identify Stable ROIs (LV1)
fprintf('\nIdentifying LV1 stable brain regions...\n');
lv1_bsr = brain_bsr(:, 1);
stable_idx = abs(lv1_bsr) > bsr_thresh;
fprintf('LV1 stable ROIs (|BSR| > %.2f): %d / %d\n', bsr_thresh, sum(stable_idx), length(stable_idx));

%% 8. Print Summary & Save Results
fprintf('\n================== PLSC RESULTS SUMMARY ==================\n');
fprintf('LV\tSingVal\t\tExplained(%%)\tPerm_P\n');
for lv = 1:length(singvals)
    fprintf('%d\t%.6f\t%.2f\t\t%.5f\n', lv, singvals(lv), explained(lv)*100, perm_pvals(lv));
end
fprintf('==========================================================\n');

results = struct();
results.target_biotype = TARGET_BIOTYPE;
results.clinical_vars_used = clinical_vars;
results.biotype_idx_before_missing = biotype_idx;  
results.valid_idx_within_biotype = valid_idx;      
results.subIDs = results.subIDs(valid_idx);

results.Xz = Xz;
results.Yz = Yz;
results.U = U;
results.V = V;
results.singvals = singvals;
results.explained = explained;            
results.crosscov = plsc.crosscov;
results.brain_scores = plsc.brain_scores;
results.clin_scores = plsc.clin_scores;

results.perm_results = perm_results;
results.perm_explained = perm_explained;  
results.perm_pvals = perm_pvals;          

results.boot_results = boot_results;
results.brain_bsr = brain_bsr;
results.clin_bsr = clin_bsr;
results.lv1_bsr_thresh = bsr_thresh;
results.lv1_stable_idx = stable_idx;

save_file = fullfile(OUTPUT_DIR, sprintf('PLSC_Clinical_Results_Biotype%d.mat', TARGET_BIOTYPE));
save(save_file, 'results');
fprintf('\nResults successfully saved to:\n%s\n', save_file);


