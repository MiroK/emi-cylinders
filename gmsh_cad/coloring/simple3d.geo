// This is a simple case of 2 connected cells where we want to stody domain
// decomposition. To this end domain is broken into half and the shared 
// surfaces/curves marked as well

// NOTE: needs special dolfin convert to have lines

DefineConstant[
radius = {10, Name "radius of cell"}
radius_x = {8, Name "radius of connection in x direction"}
radius_y = {6, Name "radius of connection in y direction"}
length = {100, Name "length of cell (body)"}
length_x = {4, Name "length of connection in x direction"}
length_y = {10, Name "length of connection in y direction"}
padx = {10, Name "bounding box padding in x direction"}
pady = {10, Name "bounding box padding in y direction"}
padz = {10, Name "bounding box padding in z direction"}
size_cell = {4, Name "mesh size on cell surface"}
size_box = {4, Name "mesh size on bounding box surface"}
];
// The length units here micro meters
SetFactory("OpenCASCADE");

nx = 2;
ny = 1;

//////////////////////////////////////////////////////////////////////
depth = Sqrt(radius*radius-radius_y*radius_y);
x_shift = length + 2*length_x;
y_shift = 2*depth + 2*length_y;

v = newv;
// Main
Cylinder(v) = {0-length/2, 0, 0, length, 0, 0, radius};
// X connection
Cylinder(v+1) = {length/2, 0, 0, length_x, 0, 0, radius_x};
Cylinder(v+2) = {0-length/2-length_x, 0, 0, length_x, 0, 0, radius_x};
// Y connection
Cylinder(v+3) = {0, depth, 0, 0, length_y, 0, radius_y};
Cylinder(v+4) = {0, -depth, 0, 0, -length_y, 0, radius_y};

cylinder() = BooleanUnion{ Volume{v}; Delete;}{ Volume{v+1, v+2, v+3, v+4}; Delete;};

// Define bounding for the piece
min_x = -length/2 - length_x - padx;
min_y = -depth - length_y - pady;
min_z = -radius - padz;

max_x = -length/2 - length_x + 1*x_shift ;
max_y = -depth - length_y + ny*y_shift + pady;
max_z = radius + padz;

// Exterior
box = newv;
Box(box) = {min_x, min_y, min_z, max_x-min_x, max_y-min_y, max_z-min_z};

// Done with first half
// Flip for second half
cylinder_flip() = Symmetry {-1, 0, 0, max_x} { Duplicata { Volume{cylinder}; } };
box_flip() = Symmetry {-1, 0, 0, max_x} { Duplicata { Volume{box}; } };
// Intersections
foo[] = BooleanFragments {Volume{box, box_flip}; Delete; }{Volume{cylinder_flip, cylinder}; Delete; };

cylider = foo[0];
box = foo[1];
Physical Volume(1) = {cylinder};
Physical Volume(10) = {box};

cylider_flip = foo[2];
box_flip = foo[3];
Physical Volume(2) = {cylinder_flip};
Physical Volume(20) = {box_flip};

cs = Unique(Abs(Boundary{ Volume{cylinder}; }));
fcs = Unique(Abs(Boundary{ Volume{cylinder_flip}; }));
bs = Unique(Abs(Boundary{ Volume{box}; }));
fbs = Unique(Abs(Boundary{ Volume{box_flip}; }));

//Cylinder vs exterior
cs_e[] = {};
fcs_e[] = {};
bs_fbs[] = {};
// Cyl vs exterior
ccount = #cs[];
bcount = #bs[];
For i In {0:(ccount-1)}
  ai = cs[i];
  For j In {0:(bcount-1)}
    bj = bs[j];
    If (ai == bj)
      cs_e[] += {bj};
    EndIf
  EndFor
EndFor
Physical Surface(1) = {cs_e[]};

// Flipped cyl vs exterior
ccount = #fcs[];
bcount = #fbs[];
For i In {0:(ccount-1)}
  ai = fcs[i];
  For j In {0:(bcount-1)}
    bj = fbs[j];
    If (ai == bj)
      fcs_e[] += {bj};
    EndIf
  EndFor
EndFor
Physical Surface(2) = {fcs_e[]};

// Box shared
Physical Surface(1020) = {6};
// Cylinder shared
Physical Surface(12) = {32};
// Bdry of 10 box
Physical Surface(10) = {2, 1, 4, 3, 5};
// Bdry of 20 box
Physical Surface(20) = {18, 19, 21, 20, 17};
// Tripply shared curve
Physical Curve(1021) = {3};

Characteristic Length{cs_e[]} = size_cell;
Characteristic Length{fcs_e[]} = size_cell;
Characteristic Length{6, 32} = size_cell;

Characteristic Length{2, 1, 4, 3, 5, 18, 19, 21, 20, 17} = size_box;