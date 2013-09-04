function [ c ] = converttocell( filename )
%CONVERTTOCELL Summary of this function goes here
%   Detailed explanation goes here

load(filename);

score = bmrfNetworkScoreArray;
len = bmrfNetworkLengthArray;
ids = bmrfNetworkIdArray;

for seed=1:length(score)
    for x=1:len(seed)
        c{seed,x} = ids(seed,x);
    end
end

end

