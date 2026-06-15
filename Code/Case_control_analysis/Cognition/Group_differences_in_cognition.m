%% Statistical Comparison of Behavioral Phenotypes Across Biotypes
% This script assesses the group differences in 7 cognitive measures 
% among HC, biotype1, and biotype2 using a linear regression / ANCOVA framework.
% 
% Model: Cognition ~ Group + Age + Sex + Education_level + SITE

clear; clc; close all;

INPUT_CSV_PATH = 'path/to/cogdata/Cognition_data.csv';
OUTPUT_DIR = 'path/to/results/case_control_cognition/Cognitive_Analysis_Output';

if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
end

%% 1. Load Data
fprintf('Loading cognitive dataset...\n');
Cogdata = readtable(INPUT_CSV_PATH);

% Ensure Group is a categorical variable and set HC as the reference level
Cogdata.Group = categorical(Cogdata.Group);
Cogdata.Group = reordercats(Cogdata.Group, {'HC', 'Biotype1', 'Biotype2'});

% Define cognitive variables
cog_outnames = {'TrailMaking1', 'TrailMaking2', 'SymbolDigitSubstitution', ...
                'NumericMemoryDigits', 'PairsMatching', ...
                'FluidIntelligenceScore', 'MatrixPatternCompletion'};
            
cog_names = {'Trail Making #1', 'Trail Making #2', 'Symbol Digit Substitution', ...
             'Numeric Memory (Digits)', 'Pairs Matching', ...
             'Fluid Intelligence (Score)', 'Matrix Pattern Completion'};

%% 2. Main Analysis
OmnibusTable = table();
PairwiseTable = table();

for j = 1:numel(cog_outnames)
    yvar = cog_outnames{j};
    yname = cog_names{j};

    % Extract relevant columns and remove rows with missing values
    vars_to_keep = {'Group', 'Age', 'Sex', 'Education_level', 'SITE', yvar};
    D = Cogdata(:, vars_to_keep);
    D = rmmissing(D);

    % Calculate sample sizes
    n_hc   = sum(D.Group == 'HC');
    n_sub1 = sum(D.Group == 'Biotype1');
    n_sub2 = sum(D.Group == 'Biotype2');

    D.Group = removecats(D.Group);
    D.Sex   = removecats(categorical(D.Sex));
    D.SITE  = removecats(categorical(D.SITE));

    fprintf('Analyzing: %s\n', yname);
    fprintf('Valid N = %d (HC=%d, Biotype1=%d, Biotype2=%d)\n', height(D), n_hc, n_sub1, n_sub2);

    if numel(categories(D.Group)) < 3
        warning('Less than 3 groups present for %s. Skipping.', yvar);
        continue;
    end

    % Fit models
    mdl_red  = fitlm(D, sprintf('%s ~ Age + Sex + Education_level + SITE', yvar));
    mdl_full = fitlm(D, sprintf('%s ~ Group + Age + Sex + Education_level + SITE', yvar));

    % Omnibus group effect (Nested F-test)
    cmp = nested_model_compare(mdl_red, mdl_full);

    tmp = table();
    tmp.CognitiveMeasure = string(yname);
    tmp.VarName = string(yvar);
    tmp.N_Total = height(D);
    tmp.N_HC = n_hc;
    tmp.N_Biotype1 = n_sub1;
    tmp.N_Biotype2 = n_sub2;
    tmp.F_stat = cmp.F;
    tmp.df1 = cmp.df1;
    tmp.df2 = cmp.df2;
    tmp.pValue = cmp.pValue;
    tmp.AdjR2_full = mdl_full.Rsquared.Adjusted;
    OmnibusTable = [OmnibusTable; tmp];

    % Pairwise contrasts
    coefTbl = mdl_full.Coefficients;
    coefNames = coefTbl.Properties.RowNames;
    B = coefTbl.Estimate;
    CovB = mdl_full.CoefficientCovariance;
    df = mdl_full.DFE;

    row_sub1 = find(strcmp(coefNames, 'Group_Biotype1'));
    row_sub2 = find(strcmp(coefNames, 'Group_Biotype2'));

    if isempty(row_sub1) || isempty(row_sub2)
        continue;
    end

    [est1, se1, t1, p1] = linear_contrast(B, CovB, row_sub1, df);
    [est2, se2, t2, p2] = linear_contrast(B, CovB, row_sub2, df);
    
    c3 = zeros(length(B), 1);
    c3(row_sub2) = 1;
    c3(row_sub1) = -1;
    [est3, se3, t3, p3] = linear_contrast(B, CovB, c3, df);


    Ptmp = table();
    Ptmp.CognitiveMeasure = repmat(string(yname), 3, 1);
    Ptmp.Contrast = ["HC_vs_Biotype1"; "HC_vs_Biotype2"; "Biotype1_vs_Biotype2"];
    Ptmp.Estimate = [est1; est2; est3];  
    Ptmp.SE = [se1; se2; se3];
    Ptmp.tValue = [t1; t2; t3];
    Ptmp.pValue = [p1; p2; p3];
    
    PairwiseTable = [PairwiseTable; Ptmp];
end

%% 3. FDR Correction across tests
fprintf('\nApplying Benjamini-Hochberg FDR correction...\n');

OmnibusTable.pFDR = bh_fdr(OmnibusTable.pValue);

PairwiseTable.pFDR = bh_fdr(PairwiseTable.pValue);

%% 4. Save Results
writetable(OmnibusTable, fullfile(OUTPUT_DIR, 'Cognition_Omnibus_Results.csv'));
writetable(PairwiseTable, fullfile(OUTPUT_DIR, 'Cognition_Pairwise_Results.csv'));
fprintf('Analysis complete. Results saved to %s\n', OUTPUT_DIR);


