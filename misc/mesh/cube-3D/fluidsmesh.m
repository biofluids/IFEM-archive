% this file is to read *.msh file that is generated from gmsh
% it will produce 3 files, 
% mxyz - coordinates
% mien - connectivity
% mrng - boundary face
clear;
fid=fopen('fluid.msh','r');
tline=fgetl(fid)
%fscanf(fid, '%s')
coord = fscanf(fid,'%f');
nn = coord(1);
for i=1:nn
    index=(i-1)*4+1;
    nod_num=coord(index+1);
    x(i)=coord(index+2);
    y(i)=coord(index+3);
    z(i)=coord(index+4);
end
tline=fgetl(fid);
tline=fgetl(fid);
ne=fgetl(fid);
for i=1:ne
    
%for i = 2:length(coord)
fclose(fid);
