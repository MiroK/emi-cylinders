Print.X3dPrecision = 1E-15;
Geometry.Tolerance = 1E-12;

DefineConstant[
radius = {10, Name "radius of cell"}
radius_x = {8, Name "radius of connection in x direction"}
radius_y = {6, Name "radius of connection in y direction"}
length = {100, Name "length of cell (body)"}
length_x = {4, Name "length of connection in x direction"}
length_y = {10, Name "length of connection in y direction"}
padz = {20, Name "bounding box padding in z direction"}
];

// The length units here micro meters
SetFactory("OpenCASCADE");

//////////////////////////////////////////////////////////////////////
depth = Sqrt(radius*radius-radius_y*radius_y);
x_shift = length + 2*length_x;
y_shift = 2*depth + 2*length_y;

min_x = -length/2 - length_x;
min_y = -depth - length_y;
min_z = -radius - padz;

max_x = length/2 + length_x;
max_y = depth + length_y;
max_z = radius + padz;

v = newv;
// The cylinder
Cylinder(v) = {0-length/2-length_x, 0, 0, length+2*length_x, 0, 0, radius};
Cylinder(v+1) = {0, -depth-length_y, 0, 0, y_shift, 0, radius_y};
cylinder() = BooleanUnion{ Volume{v}; Delete;}{ Volume{v, v+1}; Delete;};

// Bounding box
box = newv;
Box(box) = {min_x, min_y, min_z, max_x-min_x, max_y-min_y, max_z-min_z}; 

cell_bdry() = Unique(Abs(Boundary{ Volume{cylinder}; }));
box_bdry() = Unique(Abs(Boundary{ Volume{box}; }));

v() = BooleanFragments {Volume{box}; Delete; }{Volume{cylinder}; Delete; };

cylinder = v[0];
box = v[1];

// Periodicity maps
surfMaster = 10;
surfSlave = 11;
boundMaster[] = {5};
boundSlave[] = {17};
Periodic Surface surfSlave { boundSlave[] } = surfMaster { boundMaster[] };

surfMaster = 13;
surfSlave = 12;
boundMaster[] = {9};
boundSlave[] = {15};
Periodic Surface surfSlave { boundSlave[] } = surfMaster { boundMaster[] };

surfMaster = 1;
surfSlave = 7;
boundMaster[] = {1, 2, 4, 3, 5};
boundSlave[] = {7, 11, 14, 13, 17};
Periodic Surface surfSlave { boundSlave[] } = surfMaster { boundMaster[] };

surfMaster = 5;
surfSlave = 2;
boundMaster[] = {4, 10, 14, 12, 15};
boundSlave[] = {1, 8, 7, 6, 9};
Periodic Surface surfSlave { boundSlave[] } = surfMaster { boundMaster[] };

// Physical volumes and surfaces
Physical Volume(1) = {cylinder};
Physical Volume(2) = {box};

interfaces[] = Unique(Abs(Boundary{ Volume{cylinder}; }));  
boundary[] = Unique(Abs(Boundary{ Volume{box}; }));
boundary[] -= {interfaces[]};

Physical Surface(2) = {boundary[]};
Physical Surface(1) = {interfaces[]};