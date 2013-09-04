function [ filename ] = savesif(data, filename )
%SIF Summary of this function goes here
%   Detailed explanation goes here
% args:
%   data - N x 2 matrix of interaction
%   filename - output filename

    fid = fopen(filename, 'w');

    dim = size(data);
    for id=1:dim(1)
        fprintf(fid,'%i pp %i\n', data(id,1), data(id,2));
    end
    fclose(fid);

end

