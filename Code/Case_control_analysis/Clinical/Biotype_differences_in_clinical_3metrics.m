%% T2DM Biotype Comparison for Patient-Specific Clinical Metrics
% Evaluates group differences between Biotype 1 and Biotype 2 for metrics
% exclusive to the patient cohort (Duration, HOMA2B, HOMA2IR).
%
% Statistical Model: Metric ~ Group + Age + Sex + SITE

clear; clc; close all;

% Expected columns: SubID, Group (1=Biotype1, 2=Biotype2), Age, Sex, SITE, Duration_Year, HOMA2B, HOMA2IR
INPUT_CSV_PATH = 'path/to/clindata/PatientOnly_Metrics.csv';
OUTPUT_DIR = 'path/to/results/case_control_clinical/PatientOnly_Comparison_Output';
if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
end

%% 1. Load Cleaned Patient Data
fprintf('Loading patient dataset...\n');
patient_data = readtable(INPUT_CSV_PATH);

% Ensure covariates and Group are categorical
patient_data.Group = categorical(patient_data.Group);
patient_data.Sex = categorical(patient_data.Sex);
patient_data.SITE = categorical(patient_data.SITE);

%% 2. Define Metrics and Covariates
target_metrics = {'Duration_Year', 'HOMA2B', 'HOMA2IR'};
covariates = {'Age', 'Sex', 'SITE'};

%% 3. Main Analysis (Biotype 2 vs Biotype 1)
biotype_comp_results = table();

fprintf('Running statistical analysis for patient-specific profiles...\n');

for i = 1:length(target_metrics)
    var_name = target_metrics{i};
    
    % Extract current metric and covariates, then drop missing rows
    temp_table = patient_data(:, [{'SubID', var_name}, covariates, {'Group'}]);
    temp_table.Properties.VariableNames{var_name} = 'TargetVar';
    temp_table(any(ismissing(temp_table), 2), :) = [];
    temp_table.SITE = removecats(temp_table.SITE);
    
    % Calculate valid sample sizes
    n1 = sum(temp_table.Group == '1');
    n2 = sum(temp_table.Group == '2');
    fprintf('-> Metric %s: Valid N (Sub1=%d, Sub2=%d)\n', var_name, n1, n2);
    
    % Initialize statistical variables
    beta_val = NaN; t_val = NaN; p_val = NaN;
    
    if n1 > 1 && n2 > 1
        % Fit Linear Model
        lm = fitlm(temp_table, 'TargetVar ~ Group + Age + Sex + SITE');
        
        % Extract group effect (Group_2 coefficient = Biotype 2 - Biotype 1)
        coef_names = lm.CoefficientNames;
        g2_idx = find(strcmp(coef_names, 'Group_2'));
        
        if ~isempty(g2_idx)
            beta_val = lm.Coefficients.Estimate(g2_idx); 
            t_val = lm.Coefficients.tStat(g2_idx);       
            p_val = lm.Coefficients.pValue(g2_idx);      
        end
    else
        warning('Not enough subjects in both biotypes for metric %s. Skipping stat test.', var_name);
    end
    
    % Append to results table
    biotype_comp_results = [biotype_comp_results; table({var_name}, n1, n2, beta_val, t_val, p_val, ...
        'VariableNames', {'Metric', 'N_Sub1', 'N_Sub2', 'Estimate', 't_value', 'P_uncorr'})];
end

%% 4. FDR Correction
fprintf('\nApplying Benjamini-Hochberg FDR correction...\n');
if ~isempty(biotype_comp_results)
    biotype_comp_results.P_FDR = bh_fdr(biotype_comp_results.P_uncorr);
end

% Save results
output_file = fullfile(OUTPUT_DIR, 'PatientSpecific_Biotype_Comparison_Results.csv');
writetable(biotype_comp_results, output_file);
fprintf('Analysis complete. Results saved to: %s\n', output_file);


