% given gid, find the sppi in ppi

function sppi = getppisubnet(ppi, gid)

    index1 = [];
    index2 = [];
    for i=1:length(gid)
        index1 = [index1 find(ppi(:,1) == gid(i))'];
        index2 = [index2 find(ppi(:,2) == gid(i))'];
    end

    index = intersect(index1, index2);

    sppi = ppi(index,:);

end

        