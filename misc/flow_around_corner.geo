// Flat 3D domain for "2D" simulation

x1 = 1;
x2 = 0.5;
x3 = 1;

y1 = 0.7;
y2 = 0.3;
y3 = 1;

z1 = 0.1;

char_length = 0.035;

phys_vol   = 1000;

//p1  = newp; Point(p1)  = {-x1-x2,0     ,0,char_length * 4};
//p2  = newp; Point(p2)  = {   -x2,0     ,0,char_length * 4};
//p3  = newp; Point(p3)  = {0     ,0     ,0,char_length * 0.8};
//p4  = newp; Point(p4)  = {0     ,-y2   ,0,char_length * 4};
//p5  = newp; Point(p5)  = {0     ,-y2-y1,0,char_length * 4};
//p6  = newp; Point(p6)  = {    x3,-y2-y1,0,char_length * 4};
//p7  = newp; Point(p7)  = {    x3,-y2   ,0,char_length * 4};
//p8  = newp; Point(p8)  = {    x3,0     ,0,char_length * 4};
//p9  = newp; Point(p9)  = {    x3, y3   ,0,char_length * 4};
//p10 = newp; Point(p10) = {0     , y3   ,0,char_length * 4};
//p11 = newp; Point(p11) = {   -x2, y3   ,0,char_length * 4};
//p12 = newp; Point(p12) = {-x1-x2, y3   ,0,char_length * 4};

p1  = newp; Point(p1)  = {-x1-x2,0     ,0,char_length * 0.7};
p2  = newp; Point(p2)  = {   -x2,0     ,0,char_length * 0.7};
p3  = newp; Point(p3)  = {0     ,0     ,0,char_length * 0.3};
p4  = newp; Point(p4)  = {0     ,-y2   ,0,char_length * 1};
p5  = newp; Point(p5)  = {0     ,-y2-y1,0,char_length * 2};
p6  = newp; Point(p6)  = {    x3,-y2-y1,0,char_length * 6};
p7  = newp; Point(p7)  = {    x3,-y2   ,0,char_length * 6};
p8  = newp; Point(p8)  = {    x3,0     ,0,char_length * 6};
p9  = newp; Point(p9)  = {    x3, y3   ,0,char_length * 6};
p10 = newp; Point(p10) = {0     , y3   ,0,char_length * 5};
p11 = newp; Point(p11) = {   -x2, y3   ,0,char_length * 5};
p12 = newp; Point(p12) = {-x1-x2, y3   ,0,char_length * 6};


l1  = 50; Line(l1 ) = {p1,p2};
l2  = newreg; Line(l2 ) = {p2,p3};
l3  = newreg; Line(l3 ) = {p3,p4};
l4  = newreg; Line(l4 ) = {p4,p5};
l5  = newreg; Line(l5 ) = {p5,p6};
l6  = newreg; Line(l6 ) = {p6,p7};
l7  = newreg; Line(l7 ) = {p7,p8};
l8  = newreg; Line(l8 ) = {p8,p9};
l9  = newreg; Line(l9 ) = {p9,p10};
l10 = newreg; Line(l10) = {p10,p11};
l11 = newreg; Line(l11) = {p11,p12};
l12 = newreg; Line(l12) = {p12,p1};

ll1 = newreg; Line Loop(ll1) = { l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12};
//ll2 = newreg; Line Loop(ll2) = {-l3 , l11, l7 ,-l12};
//ll3 = newreg; Line Loop(ll3) = { l5 , l6 , l7 , l8 };
//ll4 = newreg; Line Loop(ll4) = { l1 , l10,-l5 ,-l9 };
//ll5 = newreg; Line Loop(ll5) = { l10, l6 ,-l11,-l2 };
//ll6 = newreg; Line Loop(ll6) = { l9 ,-l8 ,-l12, l4 };

plane1 = newreg; Plane Surface(plane1) = {ll1};
//plane2 = newreg; Plane Surface(plane2) = {ll2};
//plane3 = newreg; Plane Surface(plane3) = {ll3};
//plane4 = newreg; Plane Surface(plane4) = {ll4};
//plane5 = newreg; Plane Surface(plane5) = {ll5};
//plane6 = newreg; Plane Surface(plane6) = {ll6};

//sl1 = newreg; Surface Loop(sl1) = {-plane1,plane2,plane3,-plane4,plane5,-plane6};
//vol1 = newreg; Volume(vol1) = sl1;

Extrude Surface {63, {0.0,0.0,z1}};
sl1 = newreg; Surface Loop(sl1) = {125,-63,80,84,88,92,96,100,104,108,112,116,120,124};
vol1 = newreg; Volume(vol1) = sl1;

Physical Volume(phys_vol) = {vol1};
Physical Surface(1) = {125};
Physical Surface(2) = {63};
Physical Surface(3) = {124};
Physical Surface(4) = {88,92};
Physical Surface(5) = {100,104,108};
Physical Surface(6) = {112,116,120};
Physical Surface(7) = {80,84};
Physical Surface(8) = {96};


