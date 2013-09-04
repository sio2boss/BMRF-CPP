function [] = viewnet(filename)

%cyto_path='/Applications/Cytoscape_v2.8.3/cytoscape.sh';
cyto_path='/home/sio2/Applications/Cytoscape_v2.8.3/cytoscape.sh';
%cyto_path='/home/sio2/Applications/cytoscape_v3.0.0-beta1/cytoscape.sh';
system([cyto_path ' --network ' filename ' &']);

end