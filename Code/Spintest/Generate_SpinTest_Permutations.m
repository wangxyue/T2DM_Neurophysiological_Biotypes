%% Parallelized Spin Test Partition Rotation 
% This script executes spherical rotations of surface parcellation labels
% (e.g., DK318) on the fsLR 32k template sphere. The generated null models
% preserve the original spatial autocorrelation structure of the brain maps.

clear; clc; close all;

SPHERE_SURF_LH = 'Q1-Q6_RelatedParcellation210.L.sphere.32k_fs_LR.surf.gii';
SPHERE_SURF_RH = 'Q1-Q6_RelatedParcellation210.R.sphere.32k_fs_LR.surf.gii';

LABEL_GII_LH   = 'lh.DK318.label.gii';
LABEL_GII_RH   = 'rh.DK318.label.gii';

OUTPUT_DIR     = 'path/to/Spintest/results/Perm_DK318_Label';

NUM_PERM       = 10000; 

%% 1. Environment Preparation & Parallel Pool Initialization
poolobj = gcp('nocreate');
if isempty(poolobj)
    fprintf('Starting local parallel pool for accelerated processing...\n');
    parpool('local'); 
end

if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
end

%% 2. Load Data
SphereSurf = cell(2,1);
SphereSurf{1} = SPHERE_SURF_LH;
SphereSurf{2} = SPHERE_SURF_RH;

LV = gifti(LABEL_GII_LH); 
RV = gifti(LABEL_GII_RH);

NumHemi = size(LV.cdata, 1);

ComLR = [LV.cdata; RV.cdata];
LabelLR = ComLR;

Label = cell(2,1);
Label{1} = LabelLR(1:NumHemi, 1);
Label{2} = LabelLR(NumHemi+1:end, 1);

%% 3. Execute Spatial Rotations
fprintf('\nLaunching %d parallel spatial permutation rotations...\n', NUM_PERM);
tic; 

parfor i = 1:NUM_PERM
    OutFile = fullfile(OUTPUT_DIR, ['Rand_' sprintf('%.5d', i)]);
    GetRotateLabel(SphereSurf, Label, OutFile);

    if mod(i, 500) == 0
        fprintf('-> Completed %d / %d permutations...\n', i, NUM_PERM);
    end
end
elapsedTime = toc; 

%% 4. Summary
fprintf('\nSuccess: All %d spatial permutations generated!\n', NUM_PERM);
fprintf('Total Processing Time: %.2f seconds (approx. %.2f minutes).\n', elapsedTime, elapsedTime/60);
fprintf('Permutation datasets saved to: %s\n', OUTPUT_DIR);
