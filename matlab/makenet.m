function [ net ] = makenet( ids, dist )
%MAKENET Summary of this function goes here
%   Detailed explanation goes here
    net = struct('id', [], 'distance', []);
    
    if nargin > 1
        net.id = ids;
    end
    
    if nargin == 2
        net.distance = dist;
    end
end

