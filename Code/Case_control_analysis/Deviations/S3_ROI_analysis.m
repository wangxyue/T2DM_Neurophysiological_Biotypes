%% Case-Control Analysis:£¨Regional-level£©
% Biotype1 vs HC
% Biotype2 vs HC

clear; clc; close all;

BIOTYPING_RESULTS_PATH = 'path/to/data/Final_Biotyping_Results.mat';
HC_DEVIATION_PATH = 'path/to/data/CN_deviation_z_matrix.mat';
OUTPUT_DIR = 'path/to/results/biotype_vs_HC_output/';

load(BIOTYPING_RESULTS_PATH); 
load(HC_DEVIATION_PATH);

Z_T2DM = X;
Z_HC   = CN_deviation_z_matrix;

Z_biotype1 = Z_T2DM(final_idx_NbClust == 1, :);
Z_biotype2 = Z_T2DM(final_idx_NbClust == 2, :);

if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
end

%% Biotype1 vs HC
results_biotype1 = run_regional_level_case_control_analysis( ...
    Z_biotype1, Z_HC, 'Biotype1', OUTPUT_DIR);

%% Biotype2 vs HC
results_biotype2 = run_regional_level_case_control_analysis( ...
    Z_biotype2, Z_HC, 'Biotype2', OUTPUT_DIR);

fprintf('\nRegional-level analysis completed successfully.\n');

