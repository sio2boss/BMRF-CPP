% netscore = networkscore(geneid, cand_network_id, Zscore, ppi,
% lambda, gamma)
% calculate network energy given genes in ppi network
% Input: 
%  geneid: gene id vector
%  cand_network_id: gene ids in the interested network
%  Zscore: z-score of gene ids in geneid
%  ppi: protein-portein interaction network
%  lambda, gamma: two non-negative parameters
% Ouput:
%  netscore: network potentials, negative of netework score
function netscore = mrfnetscore(geneid, cand_network_id, Zscore, ppi, lambda, gamma)

    if nargin < 6
        gamma = 1;
    end
    if nargin < 5
        lambda = 1;
    end

    sppi = getppisubnet(ppi, cand_network_id);
    sgeneid = unique(sppi(:));
    [~,b] = intersect(geneid, sgeneid);
    sZscore = Zscore(b)';
    
    ng = length(sgeneid);
    ne = length(sppi);
    L = sparse(ng,ng);
    for i = 1:ng
        gconn = union(sppi(find(sppi(:,1)==sgeneid(i)),2), sppi(find(sppi(:,2)==sgeneid(i)),1));
        [a,igconn] = intersect(sgeneid, gconn);
        L(i,igconn) = -1;
        L(igconn,i) = -1;
        L(i,i) = length(a);
    end
    LL = L;
    for i=1:ng
        L(i,:) = L(i,:)/sqrt(LL(i,i));
    end
    for i=1:ng
        L(:,i) = L(:,i)/sqrt(LL(i,i));
    end

    x = pinv(2*lambda*L/ne+gamma*eye(size(L))/ng)*(1+sZscore*gamma)/ng;
    netscore = -mean(x) + lambda*x'*L*x/ne + gamma*(x-sZscore)'*(x-sZscore)/(2*ng);

end
