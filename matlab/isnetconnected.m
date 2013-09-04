function [ isc ] = isnetconnected(S)
%ISNETCONNECTED Determines if network is connected with eigenvalue decomp.
%  INPUT: S - connectivity matrix

    isc = false;
    
    M = length(S);
    delta = zeros(M,M);
    for i=1:M
        delta(i,i) = length(find(S(i,:)==1));
    end

    % Use eigenvalues to determine if network is connected
    eig_value=sort(eig(delta-S), 'ascend');
    if eig_value(2)>0.001
        isc = true;
    end

end