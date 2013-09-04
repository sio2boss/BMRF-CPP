%% Validation

%% Load
disp 'Loading protein-protein-interaciton data..';
load('data/simulated_ppi.mat', 'ppiArray');

disp('Loading results to save time...');
load('data/results_matlab.mat','bmrfNetworkIdArray', 'bmrfNetworkScore');

disp('Loading truth genes...');
load('data/significant_genes.mat', 'sigGeneIdArray');

%% Visualize
vis

%% precision-recall of identified subnetwork compared with ground truth
% network
disp 'Recalling identified subnetworks...';
nMatchingGeneIds = length(intersect(sigGeneIdArray, bmrfNetworkIdArray{1}));

precision = nMatchingGeneIds/length(sigGeneIdArray);
recall = nMatchingGeneIds/length(bmrfNetworkIdArray{1});

disp '';
disp(['Results: precision=' num2str(precision) ' recall=' num2str(recall)]);
disp '';
