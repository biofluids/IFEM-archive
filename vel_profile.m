fid=fopen('fem.vel00010','r');
l=fgetl(fid);
v(2010)=0.0;

%for i=1:2010
c=fscanf(fid, '%e');
%v(i)=c;
%vel(i)=c;
%end
fclose(fid);
dim=3;
nn=length(c)/dim;
for i=1:nn
    vel_x(i)=c((i-1)*dim+1);
    vel_y(i)=c((i-1)*dim+2);
    vel_z(i)=c((i-1)*dim+3);
end
nn_x=101;
nn_y=21;

j=15;dx=0.1; dy=0.1;
for k=1:nn_y
    nn_j(k)=(k-1)*nn_x+j;
    y(k)=dy*(k-1);
end

plot(vel_x(nn_j(:)), y)
