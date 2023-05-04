// Fine
// xa  = DefineNumber[0.02];
// xdi = DefineNumber[0.25];
// xdo = DefineNumber[1.5];

xa  = DefineNumber[0.01];
xdi = DefineNumber[0.075];
xdo = DefineNumber[0.075];

// Angle of attack in radians
angle = DefineNumber[-8*Pi/180];
trs   = DefineNumber[-3];

Point(1) = {0.0, 0.5, 0, xa}; 
Point(2) = {0.0, -0.5, 0, xa}; 
Point(3) = {-1.0, 0.0, 0, xa}; 

Line(1) = {1,2}; 
Line(2) = {2,3};
Line(3) = {3,1};

//+
Translate {trs, 0, 0} {Line{1}; Line{2}; Line{3}; }
// Rotate {{0, 0, 1}, {0, 0, 0}, angle} {Curve{1}; Curve{2}; }

Point(800) = { 0,  0.0, 0.0, xdo}; 
Point(801) = { 0.0, -10,  0.0, xdo}; 
Point(802) = { 5.0, -10,  0.0, xdo}; 
Point(803) = { 5.0,  10,  0.0, xdo}; 
Point(804) = { 0.0,  10,  0.0, xdo}; 

Point(805) = { trs, -2,  0.0, xdi}; 
Point(806) = { 5.0, -3,  0.0, xdi}; 
Point(807) = { 5.0,  3,  0.0, xdi}; 
Point(808) = { trs,  2,  0.0, xdi}; 

Point(809) = { -10.0,  0,  0.0, xdo}; 
Point(810) = { trs,  0.0,  0.0, xdi}; 

Line(4)  = {801,802}; 
Line(5)  = {802,806}; 
Line(6)  = {806,807}; 
Line(7)  = {807,803}; 
Line(8)  = {803,804}; 
Line(9)  = {805,806}; 
Line(10) = {807,808}; 

Circle(11) = {804, 800, 809};
Circle(12) = {809, 800, 801};
Circle(13) = {808, 810, 805};


Curve Loop(1) = {1, 2, 3};
Curve Loop(2) = {9, 6, 10, 13};
Curve Loop(3) = {4, 5, 6, 7, 8, 11, 12};

Plane Surface(1) = {2, 3};
Plane Surface(2) = {1, 2};
Coherence;

Physical Line("Wall",1)     = {1,2, 3};
Physical Line("Inflow",2)   = {4, 8, 11,12};
Physical Line("Outflow",3)  = {5, 6, 7};
Physical Surface("Domain",9) = {1,2};

Coherence;


