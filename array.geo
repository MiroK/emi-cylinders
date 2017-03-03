n_cylinders = 4;
// The first cylinder has center at x0, y0, z0
x0 = 0;
y0 = 0;
z0 = 0;
// The height of each cylinder is H with char mesh size size_c
R = 0.5;
H = 1.5;
size_R = 0.2;
// The height of the joint regions is h with char mesh size size_r
r = 0.25;
h = 0.5;
size_r = 0.2;

// The cylinders are enclosed in a bbox which leaves dx, dy gaps around the
// square that bounds the cylinder crossection. In z direction the gap is dz.
// For bbox the char size is size_b
dx = 0.2;
dy = 0.2;
dz = 0.2;
size_b = 0.3;
// Each cylinder volume gets number 1, ..., n_cylinders. These markers are also
// inherited by surfaces that bound the cylinder except the surface that is
// shared by 2 cylinders. The label for the cylinder is (volume label +
// n_cylinders). The outer volume, i.e. bounding box - cylinder union is tagged
// as 0
//! Cut here

// We draw the first guy by hand
// Joint
p = newp;
Point(p) = {x0, y0, z0-H/2-h/2, size_r};
Point(p+1) = {x0+r, y0, z0-H/2-h/2, size_r};
Point(p+2) = {x0, y0+r, z0-H/2-h/2, size_r};
Point(p+3) = {x0-r, y0, z0-H/2-h/2, size_r};
Point(p+4) = {x0, y0-r, z0-H/2-h/2, size_r};
// Base outer
p = newp;
Point(p) = {x0, y0, z0-H/2, size_R};
Point(p+1) = {x0+R, y0, z0-H/2, size_R};
Point(p+2) = {x0, y0+R, z0-H/2, size_R};
Point(p+3) = {x0-R, y0, z0-H/2, size_R};
Point(p+4) = {x0, y0-R, z0-H/2, size_R};
// Base inner
p = newp;
Point(p) = {x0, y0, z0-H/2, size_r};
Point(p+1) = {x0+r, y0, z0-H/2, size_r};
Point(p+2) = {x0, y0+r, z0-H/2, size_r};
Point(p+3) = {x0-r, y0, z0-H/2, size_r};
Point(p+4) = {x0, y0-r, z0-H/2, size_r};
// Surface
Circle(1) = {4, 1, 3};
Circle(2) = {3, 1, 2};
Circle(3) = {2, 1, 5};
Circle(4) = {5, 1, 4};
Circle(5) = {14, 6, 13};
Circle(6) = {13, 6, 12};
Circle(7) = {12, 6, 15};
Circle(8) = {15, 6, 14};
Circle(9) = {9, 6, 8};
Circle(10) = {8, 6, 7};
Circle(11) = {7, 6, 10};
Circle(12) = {10, 6, 9};
Line Loop(13) = {4, 1, 2, 3};
Plane Surface(14) = {13};
Line Loop(15) = {12, 9, 10, 11};
Line Loop(16) = {8, 5, 6, 7};
Plane Surface(17) = {15, 16};
Line(18) = {15, 5};
Line(19) = {12, 2};
Line(20) = {3, 13};
Line(21) = {4, 14};
Line Loop(22) = {18, -3, -19, 7};
Ruled Surface(23) = {22};
Line Loop(24) = {2, -19, -6, -20};
Ruled Surface(25) = {24};
Line Loop(26) = {5, -20, -1, 21};
Ruled Surface(27) = {26};
Line Loop(28) = {8, -21, -4, -18};
Ruled Surface(29) = {28};
///////////////////////////////////////////////////////////////////////////////
// Mirror on the top
///////////////////////////////////////////////////////////////////////////////
// Joint
p = newp;
Point(p) = {x0, y0, z0+H/2+h/2, size_r};
Point(p+1) = {x0+r, y0, z0+H/2+h/2, size_r};
Point(p+2) = {x0, y0+r, z0+H/2+h/2, size_r};
Point(p+3) = {x0-r, y0, z0+H/2+h/2, size_r};
Point(p+4) = {x0, y0-r, z0+H/2+h/2, size_r};
// Base outer
p = newp;
Point(p) = {x0, y0, z0+H/2, size_R};
Point(p+1) = {x0+R, y0, z0+H/2, size_R};
Point(p+2) = {x0, y0+R, z0+H/2, size_R};
Point(p+3) = {x0-R, y0, z0+H/2, size_R};
Point(p+4) = {x0, y0-R, z0+H/2, size_R};
// Base inner
p = newp;
Point(p) = {x0, y0, z0+H/2, size_r};
Point(p+1) = {x0+r, y0, z0+H/2, size_r};
Point(p+2) = {x0, y0+r, z0+H/2, size_r};
Point(p+3) = {x0-r, y0, z0+H/2, size_r};
Point(p+4) = {x0, y0-r, z0+H/2, size_r};
Circle(30) = {19, 16, 18};
Circle(31) = {18, 16, 17};
Circle(32) = {17, 16, 20};
Circle(33) = {20, 16, 19};
Circle(34) = {29, 21, 28};
Circle(35) = {28, 21, 27};
Circle(36) = {27, 21, 30};
Circle(37) = {30, 21, 29};
Circle(38) = {24, 21, 23};
Circle(39) = {23, 21, 22};
Circle(40) = {22, 21, 25};
Circle(41) = {25, 21, 24};
Line Loop(42) = {30, 31, 32, 33};
Plane Surface(43) = {42};
Line Loop(44) = {41, 38, 39, 40};
Line Loop(45) = {37, 34, 35, 36};
Plane Surface(46) = {44, 45};
Line(47) = {30, 20};
Line(48) = {17, 27};
Line(49) = {18, 28};
Line(50) = {19, 29};
Line Loop(51) = {48, 36, 47, -32};
Ruled Surface(52) = {51};
Line Loop(53) = {35, -48, -31, 49};
Ruled Surface(54) = {53};
Line Loop(55) = {30, 49, -34, -50};
Ruled Surface(56) = {55};
Line Loop(57) = {33, 50, -37, 47};
Ruled Surface(58) = {57};
// Finally the sides
Line(59) = {9, 24};
Line(60) = {10, 25};
Line(61) = {7, 22};
Line(62) = {8, 23};
Line Loop(63) = {12, 59, -41, -60};
Ruled Surface(64) = {63};
Line Loop(65) = {60, -40, -61, 11};
Ruled Surface(66) = {65};
Line Loop(67) = {61, -39, -62, 10};
Ruled Surface(68) = {67};
Line Loop(69) = {62, -38, -59, 9};
Ruled Surface(70) = {69};
// And the volume
Surface Loop(71) = {14, 56, 54, 46, 64, 17, 70, 68, 66, 23, 29, 27, 25, 52, 58, 43};
Volume(72) = {71};

// Mark the volume of the cylinder here. Surfaces will follow
Physical Volume(1) = {72};
Physical Surface(1) = {56, 54, 46, 64, 17, 70, 68, 66, 23, 29, 27, 25, 52, 58};
// All but the first tail 
cylinder_bbox[] = {14, 56, 54, 46, 64, 17, 70, 68, 66, 23, 29, 27, 25, 52, 58};

last_volume = 72;
tail = 43;
Physical Surface(n_cylinders+1) = {tail};
For i In {2:n_cylinders}
    new_volume[] = Translate {0, 0, H+h} { Duplicata { Volume{last_volume}; } };
    last_volume = new_volume[0];

    b() = Boundary{ Volume{last_volume}; } ;
    // The 'middle' surfaces are marked 
    cylinder_shell[] = {};
    For j In {1:14}
        cylinder_shell[] += {b[j]};
    EndFor
    tail = b[15];

    cylinder_bbox[] += cylinder_shell[];

    If(i < n_cylinders)
        Physical Surface(i+n_cylinders) = {tail};
    EndIf    

    Physical Surface(i) = {cylinder_shell[]};
    Physical Volume(i) = {last_volume};
EndFor

Physical Surface(1) += {14};
Physical Surface(n_cylinders) += {tail};

// Finally collect the bounding surfaces of the cylinder
cylinder_bbox[] += {tail}; 

bbox_surfaces[] = {};
///////////////////////////////////////////////////////////////////////////////
// Bounding box
///////////////////////////////////////////////////////////////////////////////
// Bottom corner
p = newp;
dx = dx + R;
dy = dy + R;
Point(p) =   {x0+dx, y0+dy, z0-H/2-h/2-dz, size_b};
Point(p+1) = {x0+dx, y0-dy, z0-H/2-h/2-dz, size_b};
Point(p+2) = {x0-dx, y0-dy, z0-H/2-h/2-dz, size_b};
Point(p+3) = {x0-dx, y0+dy, z0-H/2-h/2-dz, size_b};
// Bottom plane
l = newl;
Line(l) = {p, p+1};
Line(l+1) = {p+1, p+2};
Line(l+2) = {p+2, p+3};
Line(l+3) = {p+3, p};
s = news;
Line Loop(s) = {l, l+1, l+2, l+3};
Plane Surface(s+1) = {s};
bbox_surfaces[] += {s+1};

// Top corner
P = newp;
up = n_cylinders*(H+h) + dz;
Point(P) =   {x0+dx, y0+dy, z0-h/2-H/2+up, size_b};
Point(P+1) = {x0+dx, y0-dy, z0-h/2-H/2+up, size_b};
Point(P+2) = {x0-dx, y0-dy, z0-h/2-H/2+up, size_b};
Point(P+3) = {x0-dx, y0+dy, z0-h/2-H/2+up, size_b};
// Top plane;
L = newl;
Line(L) = {P, P+1};
Line(L+1) = {P+1, P+2};
Line(L+2) = {P+2, P+3};
Line(L+3) = {P+3, P};
s = news;
Line Loop(s) = {L, L+1, L+2, L+3};
Plane Surface(s+1) = {s};
bbox_surfaces[] += {s+1};

// Side surfaces
ls = newl;
Line(ls) = {p, P};
Line(ls+1) = {p+1, P+1};
Line(ls+2) = {p+2, P+2};
Line(ls+3) = {p+3, P+3};

s = news;
Line Loop(s) = {L, -(ls+1), -l, ls};
Plane Surface(s+1) = {s};
bbox_surfaces[] += {s+1};

s = news;
Line Loop(s) = {L+1, -(ls+2), -(l+1), (ls+1)};
Plane Surface(s+1) = {s};
bbox_surfaces[] += {s+1};

s = news;
Line Loop(s) = {L+2, -(ls+3), -(l+2), (ls+2)};
Plane Surface(s+1) = {s};
bbox_surfaces[] += {s+1};

s = news;
Line Loop(s) = {L+3, -(ls), -(l+3), (ls+3)};
Plane Surface(s+1) = {s};
bbox_surfaces[] += {s+1};

// Outer = bbox - cylinders
v = newv;
Surface Loop(v) = {bbox_surfaces[]};
Surface Loop(v+1) = {cylinder_bbox[]};
Volume(v+2) = {v, v+1};
Physical Volume(0) = {v+2};
