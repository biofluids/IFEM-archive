clear;
x=[10 20 50 100 500 1000];
y=[7.437 3.7185 1.4872 .7437 .14874 0.07437];
y2=[3.3 2.28 1.57 1.41 1.36 1.298];
semilogx(x,y,x,y2);
axis([10 1e3 0 5]);
grid;
xlabel('Reynolds Number','FontSize',16);
ylabel('Drag Coefficient','FontSize',16);
title('Drag Coefficient of the Flow Past a Cylinder','FontSize',16)