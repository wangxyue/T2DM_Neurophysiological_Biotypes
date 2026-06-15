%% Identifies spatial transcriptomic signatures of macroscopic brain deviation.

clear; clc; close all;

% Target Biotype to analyze (1 = Biotype 1; 2 = Biotype 2)
TARGET_BIOTYPE = 2;   

% Statistical Parameters
n_comp = 10;              % Number of PLS components to extract
n_spin = 10000;           % Number of spin permutations for spatial null models
n_boot = 5000;            % Number of bootstraps for gene weight stability
rng_seed = 123;          
lh_rois = 1:156;          

BIOTYPING_FILE = 'path/to/data/Final_Biotyping_Results.mat';
SPIN_NULL_FILE = sprintf('path/to/Spintest/results/Perm_biotype%d_deviation.mat', TARGET_BIOTYPE);
AHBA_GENE_FILE = 'path/to/genedata/DK318_DS01_expression_final.csv';
OUTPUT_DIR     = 'path/to/results/AHBA_Gene_PLS_Results/';

rng(rng_seed);
if ~exist(OUTPUT_DIR, 'dir'), mkdir(OUTPUT_DIR); end

%% 1. Load Data
fprintf('Loading LH spatial data and complete gene expression matrices...\n');

% A. Brain Deviation Data
load(BIOTYPING_FILE); 

% B. Spin Test Null Models
load(SPIN_NULL_FILE); 

% C. AHBA Gene Expression Matrix
gene_table = readtable(AHBA_GENE_FILE);
gene_names = gene_table.Properties.VariableNames(2:end);
G_raw = table2array(gene_table(lh_rois, 2:end));

%% 2. Spatial Alignment & Standardization (LH Only)
% Extract average deviation map for the target biotype
biotype_mask = (final_idx_NbClust == TARGET_BIOTYPE);
Y_obs_full = mean(X(biotype_mask, :), 1)'; 
Y_obs_lh = Y_obs_full(lh_rois);

if TARGET_BIOTYPE == 1
    Y_perm_lh = Perm_biotype1_deviation(lh_rois, :);
else
    Y_perm_lh = Perm_biotype2_deviation(lh_rois, :);
end

% Data Extraction and Z-score Standardization
Y_clean_nozscore = Y_obs_lh;
Y_clean = zscore(Y_clean_nozscore);  

G_clean_nozscore = G_raw;
G_clean = zscore(G_clean_nozscore);

Y_perm_clean = zscore(Y_perm_lh);

%% 3. PLS Regression Analysis 
fprintf('Executing PLS Regression for Biotype %d...\n', TARGET_BIOTYPE);

[XL, YL, XS, YS, BETA, PCTVAR, MSE, stats_orig] = plsregress(G_clean, Y_clean, n_comp);

obs_r = zeros(n_comp, 1);
obs_pct_var = PCTVAR(2, :)'; 

for k = 1:n_comp
    obs_r(k) = corr(XS(:, k), Y_clean);
end

%% 4. Spin Test Significance Validation
fprintf('Validating spatial significance via %d Spin Tests...\n', n_spin);

perm_r = zeros(n_spin, n_comp);
perm_pct_var = zeros(n_spin, n_comp); 

for s = 1:n_spin
    Y_p = Y_perm_clean(:, s);
    [~, ~, XS_p, ~, ~, PCTVAR_p] = plsregress(G_clean, Y_p, n_comp);
    
    perm_pct_var(s, :) = PCTVAR_p(2, :); 
    for k = 1:n_comp
        perm_r(s, k) = corr(XS_p(:, k), Y_p);
    end
end

p_spin_r = zeros(n_comp, 1);
p_spin_varexp = zeros(n_comp, 1);

for k = 1:n_comp
    p_spin_r(k) = sum(abs(perm_r(:, k)) >= abs(obs_r(k))) / n_spin;
    p_spin_varexp(k) = sum(perm_pct_var(:, k) >= obs_pct_var(k)) / n_spin;
end

component_stats = table((1:n_comp)', obs_pct_var, p_spin_varexp, obs_r, p_spin_r, ...
    'VariableNames', {'Component', 'VarExp_Y', 'P_spin_VarExp', 'R_obs', 'P_spin_R'});

disp('PLSC Analysis Overview (First 10 Components):');
disp(component_stats);

%% 5. Save Intermediate Matrix Results
gene_weights = stats_orig.W; 
save_file = fullfile(OUTPUT_DIR, sprintf('Biotype%d_LH_PLSC_Results.mat', TARGET_BIOTYPE));

save(save_file, 'component_stats', 'XS', ...
    'Y_clean_nozscore', 'G_clean_nozscore', ...          
    'Y_clean', 'G_clean', ...                            
    'gene_weights', 'gene_names', 'TARGET_BIOTYPE', ...
    'obs_pct_var', 'p_spin_varexp', 'perm_pct_var', ...  
    'obs_r', 'p_spin_r', 'perm_r');                                 
    
fprintf('Matrix results successfully saved to: %s\n', save_file);

%% 6. Bootstrap Reliability Testing 
fprintf('\nExecuting %d Bootstraps to determine gene stability for PC%d...\n', n_boot, TARGET_PC);

n_regions = size(G_clean, 1);
n_genes = size(G_clean, 2);

boot_weights = zeros(n_boot, n_genes);
W_orig = stats_orig.W(:, TARGET_PC);

for b = 1:n_boot
    boot_indices = randsample(n_regions, n_regions, true);
    G_boot = zscore(G_clean(boot_indices, :));
    Y_boot = zscore(Y_clean(boot_indices));
    
    [~, ~, ~, ~, ~, ~, ~, stats_boot] = plsregress(G_boot, Y_boot, n_comp);
    W_boot = stats_boot.W(:, TARGET_PC);
    
    if corr(W_boot, W_orig) < 0
        W_boot = -W_boot;
    end
    
    boot_weights(b, :) = W_boot';
end

% Calculate Bootstrap Ratio (Z-score)
boot_mean = mean(boot_weights, 1);
boot_std = std(boot_weights, 0, 1);
boot_z = (boot_mean ./ boot_std)';

%% 7. Gene Ranking & Result Export
fprintf('Ranking genes and saving results...\n');
real_weights = stats_orig.W(:, TARGET_PC);
gene_ranking = table(gene_names', boot_z, boot_mean', real_weights, ...
    'VariableNames', {'GeneSymbol', 'Z_score', 'Mean_Weight', 'Weights'});

gene_ranking_descend = sortrows(gene_ranking, 'Z_score', 'descend');
gene_ranking_ascend = sortrows(gene_ranking, 'Z_score', 'ascend');

descend_name = fullfile(OUTPUT_DIR, sprintf('Biotype%d_PC%d_Gene_Ranking_Descend.csv', TARGET_BIOTYPE, TARGET_PC));
ascend_name = fullfile(OUTPUT_DIR, sprintf('Biotype%d_PC%d_Gene_Ranking_Ascend.csv', TARGET_BIOTYPE, TARGET_PC));

writetable(gene_ranking_descend, descend_name);
writetable(gene_ranking_ascend, ascend_name);

fprintf('Gene ranking export complete.\n');
fprintf('Total Genes Z > 1.96 (Sig Positive): %d\n', sum(boot_z > 1.96));
fprintf('Total Genes Z < -1.96 (Sig Negative): %d\n', sum(boot_z < -1.96));