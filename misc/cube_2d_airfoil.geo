// Flat 3D domain for "2D" simulation

ax = 10;
ay = 0.2;
az = 5.5;

char_length = 0.4;

phys_vol   = 1000;

p1 = newp; Point(p1) = {0 ,-ay/2,-az/2,char_length-0.3};
p2 = newp; Point(p2) = {ax,-ay/2,-az/2,char_length};
p3 = newp; Point(p3) = {ax,-ay/2, az/2,char_length};
p4 = newp; Point(p4) = {0 ,-ay/2, az/2,char_length};
p5 = newp; Point(p5) = {0 , ay/2,-az/2,char_length-0.3};
p6 = newp; Point(p6) = {ax, ay/2,-az/2,char_length};
p7 = newp; Point(p7) = {ax, ay/2, az/2,char_length};
p8 = newp; Point(p8) = {0 , ay/2, az/2,char_length};

l1  = 50; Line(l1 ) = {p1,p2};
l2  = newreg; Line(l2 ) = {p2,p3};
l3  = newreg; Line(l3 ) = {p3,p4};
l4  = newreg; Line(l4 ) = {p4,p1};
l5  = newreg; Line(l5 ) = {p5,p6};
l6  = newreg; Line(l6 ) = {p6,p7};
l7  = newreg; Line(l7 ) = {p7,p8};
l8  = newreg; Line(l8 ) = {p8,p5};
l9  = newreg; Line(l9 ) = {p1,p5};
l10 = newreg; Line(l10) = {p2,p6};
l11 = newreg; Line(l11) = {p3,p7};
l12 = newreg; Line(l12) = {p4,p8};

ll1 = newreg; Line Loop(ll1) = { l1 , l2 , l3 , l4 };
ll2 = newreg; Line Loop(ll2) = {-l3 , l11, l7 ,-l12};
ll3 = newreg; Line Loop(ll3) = { l5 , l6 , l7 , l8 };
ll4 = newreg; Line Loop(ll4) = { l1 , l10,-l5 ,-l9 };
ll5 = newreg; Line Loop(ll5) = { l10, l6 ,-l11,-l2 };
ll6 = newreg; Line Loop(ll6) = { l9 ,-l8 ,-l12, l4 };

plane1 = newreg; Plane Surface(plane1) = {ll1};
plane2 = newreg; Plane Surface(plane2) = {ll2};
plane3 = newreg; Plane Surface(plane3) = {ll3};
plane4 = newreg; Plane Surface(plane4) = {ll4};
plane5 = newreg; Plane Surface(plane5) = {ll5};
plane6 = newreg; Plane Surface(plane6) = {ll6};

sl1 = newreg; Surface Loop(sl1) = {-plane1,plane2,plane3,-plane4,plane5,-plane6};
vol1 = newreg; Volume(vol1) = sl1;


Physical Volume(phys_vol) = {vol1};
Physical Surface(1) = {plane1};
Physical Surface(2) = {plane2};
Physical Surface(3) = {plane3};
Physical Surface(4) = {plane4};
Physical Surface(5) = {plane5};
Physical Surface(6) = {plane6};

