clear;
fidl=fopen('lift.dat-100','r');
fidd=fopen('drag.dat-100','r');
n=151;
%n=31
t(1:n)=0;
drag(1:n)=0;
lift(1:n)=0;
%Fl=fscanf(fidl,'%f %f\n');
%Fd=fscanf(fidd,'%f %f\n');
Fl=fscanf(fidl,'%f %f %f\n');
Fd=fscanf(fidd,'%f %f %f\n');
for i=1:n
 t(i)=Fl(3*i-2);
 drag(i)=Fd(3*i-1);
 lift(i)=Fl(3*i-1);
 coefd(i)=Fd(3*i);
 coefl(i)=Fl(3*i);
end
figure(1);
plot(t(1:n),coefd(1:n))
xlabel('Time(s)','FontSize',30);
ylabel('Drag Coefficient','FontSize',30);
title('Drag Coefficient','FontSize',30)
figure(2);
plot(t(1:n),coefl(1:n)*20)
xlabel('Time(s)','FontSize',30);
ylabel('Lift Coefficient','FontSize',30);
title('Lift Coefficient','FontSize',30)
