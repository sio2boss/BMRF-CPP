%% BMRF DEMO
% Idea here is to demonstrate the robustness of the BMRF method,
% using simulated data prepared much like ER data

%% Select dataset
%dataset = 'sim';
dataset = 'loi';
%dataset = 'loi_big';

!date

%% Load PPI
disp 'Loading protein-protein-interaciton data..';
load([dataset '/ppi.mat'], 'ppiArray');

%% Load Genes
disp 'Loading simulated gene expression data...';
load([dataset '/genes.mat'], 'geneArray', 'geneIdArray', 'geneLabelArray');
geneNormalArray = log2(geneArray+4);

%% Load Seed Genes
disp 'Loading seed ids...';
load([dataset '/seed_gene_ids_only10.mat'], 'seedGeneIdArray');

%% Find subnetworks
disp 'Finding subnetworks with BMRF...'

%% bootstrapping MRF for subnetwork identification
DISTANCE = 2; TEMPERATURE = 1; BOOTSTRAPS = 100;

[bmrfNetworkIdArray, bmrfNetworkScore] = bmrf(...
    geneNormalArray, geneIdArray, geneLabelArray,...
    ppiArray, ...
    seedGeneIdArray, ...
    DISTANCE, TEMPERATURE, BOOTSTRAPS);


%% Save
disp('Saving results because that took a while...');
save([dataset '/results_matlab_' date '.mat'],'bmrfNetworkIdArray', 'bmrfNetworkScore');

!date
