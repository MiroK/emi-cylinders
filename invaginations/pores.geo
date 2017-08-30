SetFactory("OpenCASCADE");

// Main cell
Cylinder(1) = {0, 0, 0, 1, 0, 0, 0.5, 2*Pi};

// Pores
Cylinder(2) = {0.75, 0.4, 0.0, 0, 1, 0, 0.10, 2*Pi};
Cylinder(3) = {0.45, 0.4, 0.0, 0, 1, 0, 0.10, 2*Pi};
Cylinder(4) = {0.15, 0.4, 0.0, 0, 1, 0, 0.10, 2*Pi};
Cylinder(5) = {0.75, 0.0, 0.4, 0, 0, 1, 0.10, 2*Pi};
Cylinder(6) = {0.45, 0.0, 0.4, 0, 0, 1, 0.10, 2*Pi};
Cylinder(7) = {0.15, 0.0, 0.4, 0, 0, 1, 0.10, 2*Pi};

// Main \ Pores
For i In {2:7}
	BooleanDifference{ Volume{1}; Delete; }{ Volume{i}; Delete; }
EndFor

// OuterBox
Box(2) = {-0.3, -0.8, -0.7, 1.5, 1.5, 1.5};

BooleanDifference{ Volume{2}; Delete; }{ Volume{1}; }

// Surface of inner; where multipler lives
Physical Surface(1) = {33, 28, 24, 21, 31, 23, 26, 32, 37, 35, 30, 29, 34, 22, 27, 25};

// Surface of the box; for bcs
Physical Surface(2) = {39, 41, 36, 40, 38, 37};

Physical Volume(1) = {1, 2};

// Mesh as ./gmsh -3 -clscale 0.125 pores.geo
// gmsh 3.0.0 and higher is needed
