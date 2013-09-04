% [sub_network netscore netscore2 netsize] = mrfSearchNetwork(prev_network, ppi, geneid, Zscore, upper_distance)
% simulated annealing search for optimal subnetwork given seed gene using
% MAP-MRF to score the posterior of a subnetwork
% Input:
%  pre_network: previous network structure
%    .id is the gene id in the network
%    .distance is the distance to the seed gene
%  ppi: protein-protein interaction networks
%  geneid: gene id vector
%  Zscore: z-score of the genes in the 'geneid'
%  upper_distance: upper limite of searching jump
%  T: temperature for simulated annealing
% Output:
%  sub_network: current subnetwork structure
%  netscore: negative of network scores for all distinctive subnetworks during the search
%  netscore2: negative of network scores for all subnetworks during the
%  search
%  netsize: network size for all distinctive subnetworks during the search
%  

function [sub_network netscore] = mrfsearchnet(prev_network, ppi, geneid, Zscore, upper_distance, T)

    if nargin < 6
        T = 1;
    end
    
    % inital network score, the z-score of seed node
    sub_network = prev_network;
    [~,ind_gene] = intersect(geneid, prev_network.id);
    k = 1;
    netscore(k) = -Zscore(ind_gene);
    
    % Maximumal number of iteration 
    iter = 1;
    while iter<1000
        
        % Generate the candidate nodes in the network with one step further
        cand_net = netcand(ppi, geneid, prev_network, upper_distance);
        
        % We might get genes, continue if not.
        if(isempty(intersect(geneid, cand_net.id))==true)
            iter = iter + 1;
            T = 0.9*T;
            continue;
        end
        
        % Simulated Annealing method to add/delete one more gene
        % assume prior probability of selecting genes are equal
        % randomly select one for candidate genes
        [~,ind_gene] = intersect(geneid, cand_net.id);
        aidx_gene = randperm(length(ind_gene));
        idx_gene = aidx_gene(1);

        % given the candidate gene, calculate potential energy of new sub-network
        new_id = cand_net.id(idx_gene);
        new_distance = cand_net.distance(idx_gene);
        if isempty(intersect(prev_network.id, new_id))
            cand_network_id = [prev_network.id, new_id];
        else
            cand_network_id = setdiff(prev_network.id, new_id);
        end
        if length(cand_network_id) > 1
            netscore_cand = mrfnetscore(geneid, cand_network_id, Zscore, ppi);
        else 
            netscore_cand = netscore(1);
        end
        
        % if new_energy < old_energy then ADD new node
        % else ADD new node with probability of exp(-(Difference of energy)/T);
        netscore_diff = netscore_cand - netscore(k);
        if netscore_diff < 0 
            % add/delete new node into currence sub-netowrk
            if isempty(intersect(prev_network.id, new_id))
                sub_network.id = [sub_network.id, new_id];
                sub_network.distance = [sub_network.distance, new_distance];
            else
                [~,b] = intersect(sub_network.id, new_id);
                sub_network.id(b) = [];
                sub_network.distance(b) = [];
            end
            k = k+1;
            netscore(k) = netscore_cand;
        else  % ADD/DELETE new node with probability of exp(-(Difference of energy)/T);
            sprob = exp(-netscore_diff/T);
            randidx = randsample([0 1], 1, true, [1-sprob sprob]);
            
			%disp(['netscore_diff ' num2str(netscore_diff) ', sprob ' num2str(sprob) ', temp ' num2str(T) ', fix ' num2str(randidx)]);
            if randidx == 1
                if isempty(intersect(prev_network.id, new_id))
                    sub_network.id = [sub_network.id, new_id];
                    sub_network.distance = [sub_network.distance, new_distance];
                else
                    [~,b] = intersect(sub_network.id, new_id);
                    sub_network.id(b) = [];
                    sub_network.distance(b) = [];
                end
                k = k+1;
                netscore(k) = netscore_cand;
            end
        end
        
        % Get ready for next iteration
        iter = iter + 1;
        T = 0.9*T;
        prev_network = sub_network;
        
    end
    
end

