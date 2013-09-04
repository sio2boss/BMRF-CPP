
%% Visualize
% Call this bad boy after you have loaded results or run BMRF
subnetident = getppisubnet(ppiArray, bmrfNetworkIdArray{1,1}');
viewnet(savesif(subnetident, 'data/bmrf.sif'));

subnetident = getppisubnet(ppiArray, sigGeneIdArray');
viewnet(savesif(subnetident, 'data/sig.sif'));