%% Run GPU code with same imputs as matlab bmrf code

% Select dataset
dataset = 'sim';

!date

disp 'Finding subnetworks with BMRF...'

ppiFilename = [dataset '/ppi.mat'];
genesFilename = [dataset '/genes.mat'];
seedGeneFilename = [dataset '/seed_gene_ids.mat'];
outputFilename = 'results.mat';
DISTANCE = 2; TEMPERATURE = 1; BOOTSTRAPS = 100;

files = sprintf(' --ppi %s --genes %s --seeds %s --output %s', ...
    ppiFilename, genesFilename, seedGeneFilename, outputFilename);
args = sprintf(' --distance %i --temperature %i --bootstraps %i', ...
    DISTANCE, TEMPERATURE, BOOTSTRAPS);

tic
system(['LD_LIBRARY_PATH=/usr/lib64:/usr/local/cuda/lib64 ../build/bmrf-cpu ' files args],'-echo');
toc
