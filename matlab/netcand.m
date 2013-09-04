function [ cand_id ] = netcand( ppi, geneid, prev_network, upper_distance )
%NETCAND Generates network based on previous network
%   The previous network is modified based on random variable to add or
%   drop a random node.

    % Generate the candidate nodes in the network with one step further
    cand_id = makenet();
    sign = randsample([0 1],1); % to determine add or delete one node
    if (sign==1) && (length(prev_network.id)<20)
        for i=1: length(prev_network.id),
            seed_id = makenet(prev_network.id(i), prev_network.distance(i));
            cand_id_temp.id = union(ppi(ppi(:,1)==seed_id.id, 2)', ppi(ppi(:,2)==seed_id.id, 1)');
            cand_id_temp.id = intersect(geneid,cand_id_temp.id);

            %%Since we don't consider topology change at this time, all previously
            %%included nodes are removed
            cand_id_temp.id = setdiff(cand_id_temp.id, prev_network.id);
            cand_id_temp.distance(1: length(cand_id_temp.id)) = seed_id.distance+1;
            cand_id_expand.id = union(cand_id.id, cand_id_temp.id);
            for j =1: length(cand_id_expand.id)
                dist_temp1 = inf;
                dist_temp2 = inf;
                idx =find(cand_id.id==cand_id_expand.id(j));
                if(isempty(idx)==false),
                    dist_temp1 = cand_id.distance(idx);
                end
                idx = find(cand_id_temp.id==cand_id_expand.id(j));
                if(isempty(idx)==false),
                    dist_temp2 = cand_id_temp.distance(idx);
                end
                cand_id_expand.distance(j) = min(dist_temp1, dist_temp2);
            end
            cand_id=cand_id_expand;
        end
    else
        % get the outer bound of genes into candidate gene list
        if (max(prev_network.distance) > 0)
%             if(length(prev_network.id) == 2)
%                 cand_id.id = prev_network.id(2);
%             else
%                 % Determine if network is connected
%                 tempid = prev_network.id(2:end);
%                 tempdist = prev_network.distance(2:end);
%                 for count = 1 : length(tempid)
%                     sig_gid = setdiff(prev_network.id,tempid(count));
%                     sig_ppi = getppisubnet(ppi,sig_gid);
%                     S = zeros(length(sig_gid));
%                     for i = 1 : length(sig_gid)
%                         gconn = union(sig_ppi(find(sig_ppi(:,1)==sig_gid(i)),2), sig_ppi(find(sig_ppi(:,2)==sig_gid(i)),1));
%                         [~,igconn] = intersect(sig_gid, gconn);
%                         S(i,igconn) = 1;
%                         S(igconn,i) = 1;
%                     end
%                     if(isnetconnected(S))
%                         cand_id.id = [cand_id.id, tempid(count)];
%                         cand_id.distance = [cand_id.distance tempdist(count)];
%                     end
%                 end
                maxdistance = max(prev_network.distance);
                t = find(prev_network.distance == maxdistance);
                cand_id.id = [cand_id.id prev_network.id(t)];
                cand_id.distance = [cand_id.distance prev_network.distance(t)];
            end
    end
    
    % Prune all the nodes whose distance is larger than upper_distance
    if(isempty(cand_id.id)==false)
        cand_id.id = cand_id.id(cand_id.distance <= upper_distance);
        cand_id.distance = cand_id.distance(cand_id.distance <= upper_distance);
    end
    
end

