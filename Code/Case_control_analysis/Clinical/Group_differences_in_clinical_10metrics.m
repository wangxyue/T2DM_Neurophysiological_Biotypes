%% Statistical Comparison of Clinical and Metabolic Profiles Across Biotypes
% This script assesses group differences in 10 clinical/metabolic metrics 
% among HC (Group=0), Biotype1 (Group=1), and Biotype2 (Group=2).
%
% Model: Metric ~ Group + Age + Sex + SITE

clear; clc; close all;

INPUT_CSV_PATH = 'path/to/clindata/Metabolic_dataset_3groups.csv';
OUTPUT_DIR = 'path/to/results/case_control_clinical/Metabolic_Analysis_Output';
if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
end

%% 1. Load Cleaned Data
Clindata = readtable(INPUT_CSV_PATH);

% Ensure covariates and Group are categorical
Clindata.Group = categorical(Clindata.Group);
Clindata.Sex = categorical(Clindata.Sex);
Clindata.SITE = categorical(Clindata.SITE);

%% 2. Define Metrics and Covariates
% Strictly 10 metrics shared across all 3 groups
metrics = {'BMI', 'HbA1c', 'SBP', 'DBP', 'TC', 'TG', 'HDL', 'LDL', 'FPG', 'FINS'};
covariates = {'Age', 'Sex', 'SITE'};

%% 3. Main Analysis
overall_results = table(); 
posthoc_results = table();

fprintf('Running statistical analysis across 10 metabolic profiles...\n');

for i = 1:length(metrics)
    var_name = metrics{i};
    
    % Extract relevant columns and rename target variable for fitlm
    temp_table = Clindata(:, [{var_name}, covariates, {'Group'}]);
    temp_table.Properties.VariableNames{var_name} = 'TargetVar';
    
    % Remove missing data for current metric
    temp_table(any(ismissing(temp_table), 2), :) = [];
    temp_table.SITE = removecats(temp_table.SITE);
    
    % Note: Assuming HC=0, Sub1=1, Sub2=2 based on original categorical logic
    n0 = sum(temp_table.Group == '0'); 
    n1 = sum(temp_table.Group == '1'); 
    n2 = sum(temp_table.Group == '2');
    fprintf('-> Metric %s: Valid N (HC=%d, Sub1=%d, Sub2=%d)\n', var_name, n0, n1, n2);

    % Initialize output variables
    f_val = NaN; p_ov = NaN; 
    d1=NaN; t1=NaN; p1=NaN; 
    d2=NaN; t2=NaN; p2=NaN; 
    d21=NaN; t21=NaN; p21=NaN;

    % Require all 3 groups to be present to run the full model
    if n0 > 1 && n1 > 1 && n2 > 1
        % General Linear Model framework
        lm = fitlm(temp_table, 'TargetVar ~ Group + Age + Sex + SITE');
        
        % Omnibus Effect via ANOVA
        a = anova(lm); 
        f_val = a.F(strcmp(a.Row, 'Group')); 
        p_ov = a.pValue(strcmp(a.Row, 'Group'));
        
        names = lm.CoefficientNames; 
        idx1 = find(strcmp(names, 'Group_1')); 
        idx2 = find(strcmp(names, 'Group_2'));
        
        % 1) HC vs Biotype1
        if ~isempty(idx1) 
            d1 = lm.Coefficients.Estimate(idx1); 
            p1 = lm.Coefficients.pValue(idx1); 
            t1 = lm.Coefficients.tStat(idx1); 
        end
        
        % 2) HC vs Biotype2
        if ~isempty(idx2) 
            d2 = lm.Coefficients.Estimate(idx2); 
            p2 = lm.Coefficients.pValue(idx2); 
            t2 = lm.Coefficients.tStat(idx2); 
        end
        
        % 3) Biotype1 vs Biotype2
        if ~isempty(idx1) && ~isempty(idx2)
            H = zeros(1, lm.NumCoefficients); H(idx1)=-1; H(idx2)=1; 
            [p21, f21, ~] = coefTest(lm, H); 
            d21 = d2 - d1; 
            t21 = sign(d21) * sqrt(f21); 
        end
    else
        warning('Not enough subjects in all groups for metric %s. Skipping stat test.', var_name);
    end
    
    % Append to overall results
    overall_results = [overall_results; table({var_name}, n0+n1+n2, n0, n1, n2, f_val, p_ov, ...
        'VariableNames', {'Metric','Total_N','N_HC','N_Sub1','N_Sub2','F_value','P_uncorr'})];
    
    % Append to post-hoc pairwise results
    posthoc_results = [posthoc_results; table(repmat({var_name},3,1), {'HC_vs_S1';'HC_vs_S2';'S1_vs_S2'}, ...
        [d1;d2;d21], [t1;t2;t21], [p1;p2;p21], ...
        'VariableNames', {'Metric','Comparison','Difference','t_value','P_uncorr'})];
end

%% 4. FDR Correction
fprintf('\nApplying Benjamini-Hochberg FDR correction...\n');

overall_results.P_FDR = bh_fdr(overall_results.P_uncorr);
posthoc_results.P_FDR = bh_fdr(posthoc_results.P_uncorr);

%% 5. Save Results
writetable(overall_results, fullfile(OUTPUT_DIR, 'Metabolic_Overall_Comparison.csv'));
writetable(posthoc_results, fullfile(OUTPUT_DIR, 'Metabolic_PostHoc_Comparison.csv'));
fprintf('Analysis complete. Results saved to %s\n', OUTPUT_DIR);


