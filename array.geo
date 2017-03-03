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
Surface Loop(71) = {43, 56, 54, 46, 64, 17, 70, 68, 66, 23, 29, 27, 25, 14, 52, 58};
Volume(72) = {71};

// We can mark stuff now
Physical Volume(1) = {72};
Physical Surface(1) = {56, 54, 46, 64, 17, 70, 68, 66, 23, 29, 27, 25, 14, 52, 58};
Physical Surface(1+n_cylinders) = {43};

last_cylinder[] = {14, 56, 54, 46, 64, 17, 70, 68, 66, 23, 29, 27, 25, 52, 58, 43};
tip = 14;
tail = 43;
For i In {1:n_cylinders}
    out[] = Translate {0, 0, H+h} { Duplicata { Surface{last_cylinder[]}; } };

    tip = tail;
    tail = out[14];
    out_filter[] = {tip};
    For j In {1:14}
        out_filter[] += {out[j]};
    EndFor

    sl = newl;
    Surface Loop(sl) = out_filter[];
    v = newv;
    Volume(v) = {sl};
    Physical Volume(i+1) = {v};

    last_cylinder[] = out_filter[];
EndFor
