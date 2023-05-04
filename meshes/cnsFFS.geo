r = DefineNumber[ 0.075];
xmin = DefineNumber[ 0.0];
xmax = DefineNumber[ 3.0];
ymin = DefineNumber[ 0.0];
ymax = DefineNumber[ 1.0];

xstp = DefineNumber[ 0.6];
ystp = DefineNumber[ 0.2];

Point(1) = {xmin, ymin, 0.0, r};
Point(2) = {xstp, ymin, 0.0, 0.1*r};
Point(3) = {xstp, ystp, 0.0, 0.1*r};
Point(4) = {xmax, ystp, 0.0, r};
Point(5) = {xmax, ymax, 0.0, r};
Point(6) = {xmin, ymax, 0.0, r};

//+
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};
//+
Line Loop(5) = {1, 2, 3, 4, 5, 6};
//+
Plane Surface(6) = {5};
Physical Surface("Domain", 9) = {6};
//+
Physical Line("Wall", 1) = {1,2,3,5};
Physical Line("Inflow", 2) = {6};
Physical Line("Outflow", 3) = {4};
//+