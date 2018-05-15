// 2D
// This is a simple case of 2 connected cells where we want to stody domain
// decomposition. To this end domain is broken into half and the shared 
// surfaces/curves marked as well

// NOTE: needs special dolfin convert to have points?

DefineConstant[
length = {100, Name "length of cell (body)"}
length_x = {40, Name "length of connection in x direction"}
length_y = {20, Name "length of connection in y direction"}
padx = {10, Name "bounding box padding in x direction"}
pady = {10, Name "bounding box padding in y direction"}
size_cell = {4, Name "mesh size on cell surface"}
size_box = {4, Name "mesh size on bounding box surface"}
];
// The length units here micro meters
SetFactory("OpenCASCADE");

nx = 2;
ny = 1;

//////////////////////////////////////////////////////////////////////

Point(1) = {-length/2-length_x, length_y/2, 0, size_cell};
Point(2) = {-length/2, length_y/2, 0, size_cell};
Point(3) = {-length/2, length/2, 0, size_cell};
Point(4) = {-length_x/2, length/2, 0, size_cell};
Point(5) = {-length_x/2, length/2+length_y, 0, size_cell};
Point(6) = {length_x/2, length/2+length_y, 0, size_cell};
Point(7) = {length_x/2, length/2, 0, size_cell};
Point(8) = {length/2, length/2, 0, size_cell};
Point(9) = {length/2, length_y/2, 0, size_cell};
Point(10) = {length/2+length_x, length_y/2, 0, size_cell};
//
Point(20) = {-length/2-length_x, -length_y/2, 0, size_cell};
Point(19) = {-length/2, -length_y/2, 0, size_cell};
Point(18) = {-length/2, -length/2, 0, size_cell};
Point(17) = {-length_x/2, -length/2, 0, size_cell};
Point(16) = {-length_x/2, -length/2-length_y, 0, size_cell};
Point(15) = {length_x/2,-length/2-length_y, 0, size_cell};
Point(14) = {length_x/2, -length/2, 0, size_cell};
Point(13) = {length/2, -length/2, 0, size_cell};
Point(12) = {length/2, -length_y/2, 0, size_cell};
Point(11) = {length/2+length_x, -length_y/2, 0, size_cell};

l = newl;
For i In {1:19}
  Line(l) = {i, i+1};
  l += 1;
EndFor
Line(newl) = {20, 1};

// Bounding
min_x = -length/2 - length_x - padx;
min_y = -length/2 - length_y - pady;

max_x = length/2 + length_x;
max_y = length/2 + length_y+pady;

Point(21) = {min_x, min_y, 0, size_box};
Point(22) = {max_x, min_y, 0, size_box};
Point(23) = {max_x, max_y, 0, size_box};
Point(24) = {min_x, max_y, 0, size_box};

//+
Line(21) = {22, 21};
//+
Line(22) = {21, 24};
//+
Line(23) = {24, 23};
//+
Line(24) = {23, 10};
//+
Line(25) = {11, 22};
//+
Curve Loop(1) = {4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {23, 24, -9, -8, -7, -6, -5, -4, -3, -2, -1, -20, -19, -18, -17, -16, -15, -14, -13, -12, -11, 25, 21, 22};
//+
Plane Surface(2) = {2};


cylinder_flip() = Symmetry {-1, 0, 0, max_x} { Duplicata { Surface{1}; } };
box_flip() = Symmetry {-1, 0, 0, max_x} { Duplicata { Surface{2}; } };

Physical Surface(1) = {1};
Physical Surface(10) = {2};
Physical Surface(2) = {cylinder_flip[]};
Physical Surface(20) = {box_flip[]};

// Bdry of 1
Physical Curve(1) = {7, 6, 5, 4, 3, 2, 1, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 9, 8};
// Bdry of 2
Physical Curve(2) = {27, 26, 45, 44, 42, 43, 41, 40, 39, 38, 37, 36, 34, 35, 33, 31, 30, 29, 28};
// Outer bdry of 10
Physical Curve(10) = {23, 22, 21};
// Outer bdry of 20
Physical Curve(20) = {46, 69, 68};
// Shared of cylinders
Physical Curve(12) = {10};
// Shared of boxes
Physical Curve(1020) = {24, 25};
// Tripple shared
Physical Point(3) = {10, 11};
