%% Resampling-based Cluster Stability Analysis for Biotype Validation
% Evaluates the robustness of clustering solutions (k) by repeatedly subsampling 
% the dataset, applying K-means, and measuring label agreement on overlapping 
% subjects using the Hungarian matching algorithm.

clear; clc; close all;

INPUT_CSV = 'path/to/data/T2DM_deviation_z_with_SubID.csv';
OUT_DIR   = 'path/to/results/Cluster_analysis/Cluster_Stability_Results';

% Resampling & Clustering Parameters
k_list      = 2:10;             % Range of k to evaluate
nIter       = 1000;             % Number of resampling iterations
subsample_p = 0.80;             % Subsampling ratio (80%)
rng_seed    = 2026;             

% K-means Parameters
nReplicates = 10;               
maxIterKM   = 1000;             
distanceKM  = 'sqeuclidean';   

rng(rng_seed, 'twister'); 
if ~exist(OUT_DIR, 'dir'), mkdir(OUT_DIR); end

%% 1. Load and Prepare Data
fprintf('Loading data...\n');
T = readtable(INPUT_CSV);

if ~ismember('SubID', T.Properties.VariableNames)
    error('The input table must contain a column named "SubID".');
end

SubID = T.SubID;
X = T{:, 2:end}; 

if ~isnumeric(X)
    error('Feature matrix is not numeric. Please check the csv file.');
end

[nSub, nFeat] = size(X);
fprintf('Loaded data: %d subjects, %d features.\n', nSub, nFeat);

% Remove rows with NaN if any exist
nan_rows = any(isnan(X), 2);
if any(nan_rows)
    warning('Found %d subjects with NaN values. They will be removed.', sum(nan_rows));
    X = X(~nan_rows, :);
    SubID = SubID(~nan_rows, :);
    nSub = size(X,1);
end

fprintf('After NaN removal: %d subjects remain.\n', nSub);

%% 2. Feature Standardization (Z-score)
Xz = zscore(X);

%% 3. Main Resampling Stability Analysis
nK = numel(k_list);
stability_all = nan(nIter, nK);

fprintf('\nStarting resampling-based stability analysis (%d iterations)...\n', nIter);

for ik = 1:nK
    k = k_list(ik);
    fprintf('Processing k = %d ...\n', k);
    
    parfor iter = 1:nIter 
        idxA = sort(randperm(nSub, round(subsample_p * nSub)));
        idxB = sort(randperm(nSub, round(subsample_p * nSub)));
        
        XA = Xz(idxA, :);
        XB = Xz(idxB, :);

        labelA = kmeans(XA, k, 'Distance', distanceKM, 'Replicates', nReplicates, ...
            'MaxIter', maxIterKM, 'Display', 'off');
        
        labelB = kmeans(XB, k, 'Distance', distanceKM, 'Replicates', nReplicates, ...
            'MaxIter', maxIterKM, 'Display', 'off');

        [commonIdx, posA, posB] = intersect(idxA, idxB);

        if numel(commonIdx) < k
            stability_all(iter, ik) = NaN;
            continue;
        end
        
        overlapLabelA = labelA(posA);
        overlapLabelB = labelB(posB);
        
        % Align labels and compute best agreement via Hungarian algorithm
        bestAcc = cluster_agreement(overlapLabelA, overlapLabelB, k);
        stability_all(iter, ik) = bestAcc;
    end
end
fprintf('Stability analysis finished.\n');

%% 4. Summarize Results & Save CSVs
mean_stab = nanmean(stability_all, 1);
std_stab  = nanstd(stability_all, 0, 1);
sem_stab  = std_stab ./ sqrt(sum(~isnan(stability_all), 1));

[best_mean_stab, best_idx] = max(mean_stab);
best_k = k_list(best_idx);

fprintf('\n===== Stability Summary =====\n');
for ik = 1:nK
    fprintf('k = %d: mean stability = %.4f, std = %.4f\n', ...
        k_list(ik), mean_stab(ik), std_stab(ik));
end
fprintf('Most stable k by mean stability: %d (%.4f)\n', best_k, best_mean_stab);

summaryTable = table(k_list(:), mean_stab(:), std_stab(:), sem_stab(:), ...
    'VariableNames', {'k', 'MeanStability', 'StdStability', 'SEM'});
writetable(summaryTable, fullfile(OUT_DIR, 'Stability_Summary.csv'));

iterTable = array2table(stability_all, 'VariableNames', strcat('k', string(k_list)));
writetable(iterTable, fullfile(OUT_DIR, 'Stability_All_Iterations.csv'));
