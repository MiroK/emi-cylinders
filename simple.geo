n_cylinders = 4;
// The first cylinder has center at x0, y0, z0
x0 = 0;
y0 = 0;
z0 = 0;
// The hight of each cylinder is h with char mesh size size_c
r = 0.5;
h = 1.5;
size_c = 0.2;
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
///////////////////////////////////////////////////////////////////////////////
// Cells
///////////////////////////////////////////////////////////////////////////////
// We grow cylinders from base up
Point(1) = {x0, y0, z0-h/2, size_c};
Point(2) = {x0+r, y0, z0-h/2, size_c};
Point(3) = {x0-r, y0, z0-h/2, size_c};
Point(4) = {x0, y0+r, z0-h/2, size_c};
Point(5) = {x0, y0-r, z0-h/2, size_c};

// Base
Circle(1) = {4, 1, 2};
Circle(2) = {2, 1, 5};
Circle(3) = {5, 1, 3};
Circle(4) = {3, 1, 4};
Line Loop(5) = {4, 1, 2, 3};

Plane Surface(6) = {5};

cylinder_volumes[] = {};
cylinder_facets[] = {6};

last = 6;
For i In {1:n_cylinders}
    Out = Extrude {0, 0, h} { Surface{last}; };
	
    cylinder_volumes += {Out[1]};
    cylinder_facets += {Out[2], Out[3], Out[4], Out[5]};
    // Collect physical facets and volumes
    Physical Volume(i) = {Out[1]}; 
    Physical Surface(i) = {Out[2], Out[3], Out[4], Out[5]};
    // Membrane between cylinders 
    If(i < n_cylinders)
        Physical Surface(i+n_cylinders) = {Out[0]};
    EndIf    

    last = Out[0];
EndFor
// First and last surfaces need top and bottom
cylinder_facets += {last};
Physical Surface(1) += {6};
Physical Surface(n_cylinders) += {last};

bbox_surfaces[] = {};
///////////////////////////////////////////////////////////////////////////////
// Bounding box
///////////////////////////////////////////////////////////////////////////////
// Bottom corner
p = newp;
dx = dx + r;
dy = dy + r;
Point(p) = {x0+dx, y0+dy, z0-h/2-dz, size_b};
Point(p+1) = {x0+dx, y0-dy, z0-h/2-dz, size_b};
Point(p+2) = {x0-dx, y0-dy, z0-h/2-dz, size_b};
Point(p+3) = {x0-dx, y0+dy, z0-h/2-dz, size_b};
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
H = n_cylinders*h + dz;
Point(P) = {x0+dx, y0+dy, z0-h/2+H, size_b};
Point(P+1) = {x0+dx, y0-dy, z0-h/2+H, size_b};
Point(P+2) = {x0-dx, y0-dy, z0-h/2+H, size_b};
Point(P+3) = {x0-dx, y0+dy, z0-h/2+H, size_b};
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
Surface Loop(v+1) = {cylinder_facets[]};
Volume(v+2) = {v, v+1};
Physical Volume(0) = {v+2};
