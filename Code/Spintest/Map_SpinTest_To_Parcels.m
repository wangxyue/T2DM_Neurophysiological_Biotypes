%% Spatial Spin Test Mapping: Parcel-Level Null Distribution Generation
clear; clc; close all;

DATA_FILE  = 'path/to/data/Biotype1_mean_deviation.mat';
PERM_LABEL_DIR = 'path/to/Spintest/results/Perm_DK318_Label';

LABEL_GII_LH   = 'lh.DK318.label.gii';
LABEL_GII_RH   = 'rh.DK318.label.gii';

OUTPUT_MAT     = 'path/to/Spintest/results/Perm_biotype1_deviation.mat';
N_PARCELS      = 318;  

%% 1. Load Brain Metric Data
fprintf('Loading brain phenotype data...\n');
load(DATA_FILE);

%% 2. Detect and Verify Permutation Files
file_list = dir(fullfile(PERM_LABEL_DIR, 'Rand_*.mat'));
NumPerm = length(file_list);

if NumPerm == 0
    error('No permutation files (Rand_*.mat) found in directory: %s', PERM_LABEL_DIR);
end
fprintf('Detected %d permutation files. Setting NumPerm = %d\n', NumPerm, NumPerm);

%% 3. Load Template Surface
LV = gifti(LABEL_GII_LH); 
RV = gifti(LABEL_GII_RH);
DK_LR = [LV.cdata; RV.cdata];
roi_vertex_indices = cell(N_PARCELS, 1);
for node = 1:N_PARCELS
    roi_vertex_indices{node} = find(DK_LR == node);
end

%% 4. Map Vertex Rotations to Parcel Space
Perm_biotype1_deviation = zeros(N_PARCELS, NumPerm);
Perm_data = Biotype1_mean_deviation;

fprintf('\nMapping %d permutations from vertex to parcel space...\n', NumPerm);
tic;

for i = 1:NumPerm
    current_file = fullfile(PERM_LABEL_DIR, ['Rand_' sprintf('%.5d', i) '.mat']);
    load(current_file);
    Rand_Label = [LNewLabel; RNewLabel];
    for node = 1:N_PARCELS
        vertex_idx = roi_vertex_indices{node};
        perm_roi = mode(Rand_Label(vertex_idx));
        
        if perm_roi == 0
            Perm_biotype1_deviation(node, i) = 0; 
            continue;
        end
        Perm_biotype1_deviation(node, i) = Perm_data(perm_roi);
    end
    
    if mod(i, 500) == 0
        fprintf('-> Processed %d / %d permutations...\n', i, NumPerm);
    end
end

elapsedTime = toc;
fprintf('\nMapping complete! Elapsed time: %.2f seconds (approx. %.2f minutes).\n', elapsedTime, elapsedTime/60);

%% 5. Save Results
[out_folder, ~, ~] = fileparts(OUTPUT_MAT);
if ~exist(out_folder, 'dir'), mkdir(out_folder); end

save(OUTPUT_MAT, 'Perm_biotype1_deviation');
fprintf('Spatial null distribution successfully saved to:\n%s\n', OUTPUT_MAT);
