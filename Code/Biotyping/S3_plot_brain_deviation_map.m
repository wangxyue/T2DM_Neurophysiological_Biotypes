%% Plot ROI-level deviation score
load('path/to/data/Biotype1_mean_deviation.mat')
inputdata = Biotype1_mean_deviation;

% mapping to fs_LR_32k
giiLpath = 'path/to/data/lh.DK318.label.gii';
giiRpath = 'path/to/data/rh.DK318.label.gii';

gii1 = gifti(giiLpath);
gii2 = gifti(giiRpath);

ParcelLabel = double([gii1.cdata;gii2.cdata]);

data_surf = zeros(length(ParcelLabel),1);

for i = 1:length(unique(ParcelLabel))-1
    data_surf(ParcelLabel==i) = inputdata(i);
end

outpath_txt = 'path/to/results/Mean_deviation/Biotype1.txt';
save(outpath_txt,'data_surf','-ascii');

% plot
surf = 'path/to/data/FSaverage_inflated_32K.nv';
BrainNet_MapCfg(surf,outpath_txt,'dev.mat');