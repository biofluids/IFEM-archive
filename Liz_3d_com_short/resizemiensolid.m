% extract the connectivity matrix from the mien_solid.in
clear all
close all
load('mien_solid.in')
fid_out = fopen('connectsolid.in','wt');
[ne,nen]=size(mien_solid);

for i=1:ne
    for j=1:nen-1
    fprintf(fid_out,'%8d',mien_solid(i,j));
    end
  fprintf(fid_out,'\n');
end