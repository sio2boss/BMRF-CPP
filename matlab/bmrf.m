% [sub_network_struct sub_network_score] = bmrf(data, label, ppi, seed_id, Distance, T)
%  Bootstrapping Markov Random Field (BMRF) subnetwork identification by
%  integrating protein-protein interaction network and gene expression
%  profiles
% Input:
%  data: microarray gene expression data, each column is a sample and each
%  row is a gene/feature
%  geneid: gene ID for each row in 'data'
%  label: binary (1,2) labels for data samples
%  ppi: protein-protein interaction networks
%  seed_id: seed gene ids
%  Distance: parameter used in 'search_network_MRF_SA2', default is 2
%  T: temperature parameter used in'search_network_MRF_SA2', default is 1
% Output: 
%  BMRF_network_ID: all subnetworks gene ids identified by BMRF
%  BMRF_network_score: all subnetwork score identified by BMRF

function [BMRF_network_ID BMRF_network_score] = bmrf(data, geneid, label, ppi, seed_id, Distance, T, nB)

    % use defaults if constant parameters are not provided
    if nargin < 8
        nB = 100;
    end
    if nargin < 7
        T = 1;
    end
    if nargin < 6
        Distance = 2;
    end

    tic
    fprintf(1, '  Normalize / Zscore Calculation...\n');
    c1 = find(label == 1); c2 = find(label == 2);
    nc1 = sum(label == 1); nc2 = sum(label == 2);
    ns = length(label);
    search_start = 1;
    search_end = length(seed_id);

    % normalize gene expression data
    data=(data-mean(data,2)*ones(1, size(data,2)))./(std(data, 0, 2)*ones(1, size(data,2)));
    Zscore = genescore(data(:,c1)', data(:,c2)');
    Zscore0 = Zscore;

    %for i=search_start: search_end,
    %    init_network = makenet(seed_id(i), 0);
    %    idx_network = i - search_start+1;
    %    [sub_network, netscore] = mrfsearchnet(init_network, ppi, geneid, Zscore, Distance, T);
    %    sub_network_struct{idx_network} = sub_network;
    %    sub_network_score{idx_network} = netscore;
    %end
    fprintf(1, 'normalize, %f', toc);


    % Bootstrapping for subnetwork
    tic
    fprintf(1, '  Bootstrap: (%i): ', nB);
    Boot_network = cell(nB, length(seed_id));
    Boot_network_score = zeros(nB, length(seed_id));
    for ib = 1:nB
        fprintf(1, '.');
        bc1 = randsample(c1, length(c1), true);   
        bc2 = randsample(c2, length(c2), true);
        Zscore = genescore(data(:,bc1)', data(:,bc2)');

        for i=search_start: search_end,
            init_network = makenet(seed_id(i), 0);
            idx_network = i - search_start+1;
            [sub_network, netscore] = mrfsearchnet(init_network, ppi, geneid, Zscore, Distance, T);
            Boot_network{ib, idx_network} = sub_network.id;
            Boot_network_score(ib,idx_network) = -netscore(end);
            %sub_network.id
            %-netscore(end)
        end
    end
    fprintf(1, '\n');
    fprintf(1, 'bootstrap, %f', toc);

    % generating baseline for bootstrapping
    tic
    Boot_network_baseline = cell(nB, length(seed_id));
    Boot_network_score_baseline = zeros(nB, length(seed_id));
    fprintf(1, '  Bootstrap: (%i): ', nB);
    for ib = 1:nB
        fprintf(1, '.');
        temp = randperm(ns);
        bc1 = temp(1:nc1);
        bc2 = temp((nc1+1):ns);
        Zscore = genescore(data(:,bc1)', data(:,bc2)');

        for i=search_start: search_end,
            init_network = makenet(seed_id(i), 0);
            idx_network = i - search_start+1;
            [sub_network, netscore] = mrfsearchnet(init_network, ppi, geneid, Zscore, Distance, T);
            Boot_network_baseline{ib, idx_network} = sub_network.id;
            Boot_network_score_baseline(ib,idx_network) = -netscore(end);
        end
    end
    fprintf(1, '\n');
    fprintf(1, 'baseline, %f', toc);

    
    % get the confident sub-networks using credibility analysis
    tic
    disp('  Credibility Analysis...');
    for i=search_start:search_end
         idx_network = i - search_start+1;
         dect_gid_all = [];
         dect_gid_all_baseline = [];
         for ib = 1:nB
             temp = Boot_network{ib,idx_network};
             dect_gid_all = union(dect_gid_all, temp);
             temp = Boot_network_baseline{ib,idx_network};
             dect_gid_all_baseline = union(dect_gid_all_baseline, temp);
         end
         freq_dect = zeros(length(dect_gid_all),1);
         freq_dect_baseline = zeros(length(dect_gid_all_baseline),1);
         for ib = 1:nB
             temp = Boot_network{ib,idx_network};
             [~,b] = intersect(dect_gid_all,temp);
             freq_dect(b) = freq_dect(b) + 1;
             temp = Boot_network_baseline{ib,idx_network};
             [~,b] = intersect(dect_gid_all_baseline, temp);
             freq_dect_baseline(b) = freq_dect_baseline(b) +1;
         end
         % calculate the CDF of condidence
         % discard the seed gene
         confidence_obs = freq_dect/nB;
         f_obs = 0.01:0.01:1;
         for kk = 1:length(f_obs)
            cdf_obs(kk) = sum(confidence_obs>=f_obs(kk));
         end
        cdf_obs = cdf_obs/length(confidence_obs);
        confidence_baseline = freq_dect_baseline/nB;
        f_baseline = 0.01:0.01:1;
        for kk = 1:length(f_baseline)
            cdf_baseline(kk) = sum(confidence_baseline>=f_baseline(kk));
        end
        cdf_baseline = cdf_baseline/length(confidence_baseline);
        % calculate FDR
        fdr = cdf_baseline./cdf_obs;
        tt = find(fdr<=0.05);
        if ~isempty(tt)
            thrd = tt(1);
            thrd_freq = max(thrd,15);
        else 
            thrd_freq = 15;
        end
        BMRF_network_ID{idx_network} = dect_gid_all(freq_dect>=thrd_freq);
        BMRF_network_score(idx_network) = -mrfnetscore(geneid, BMRF_network_ID{idx_network}, Zscore0, ppi);
    end
    fprintf(1, 'credibility, %f', toc);
    
end

