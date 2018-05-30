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
Cylinder(v) = {0-length/2, 0, 0, length, 0, 0, radius};
Cylinder(v+1) = {0-length/2-length_x, 0, 0, length+2*length_x, 0, 0, radius_x};
Cylinder(v+2) = {0, -depth-length_y, 0, 0, y_shift, 0, radius_y};

rest[] = {v+1, v+2};
cylinder = BooleanUnion{ Volume{v}; Delete; }{ Volume{rest[]}; Delete; };
// Bounding box
box = newv;
Box(box) = {min_x, min_y, min_z, max_x-min_x, max_y-min_y, max_z-min_z}; 

v() = BooleanFragments {Volume{box}; Delete; }{Volume{cylinder}; Delete; };

cylinder = v[0];
box = v[1];

//Periodicity maps
surfMaster = 17;
surfSlave = 14;
boundMaster[] = {5};
boundSlave[] = {18};
Periodic Surface surfSlave { boundSlave[] } = surfMaster { boundMaster[] };

surfMaster = 1;
surfSlave = 7;
boundMaster[] = {2, 4, 3, 1, 5};
boundSlave[] = {11, 14, 13, 7, 18};
Periodic Surface surfSlave { boundSlave[] } = surfMaster { boundMaster[] };

surfMaster = 15;
surfSlave = 16;
boundMaster[] = {15};
boundSlave[] = {9};
Periodic Surface surfSlave { boundSlave[] } = surfMaster { boundMaster[] };

surfMaster = 2;
surfSlave = 5;
boundMaster[] = {6, 7, 8, 1, 9};
boundSlave[] = {12, 14, 10, 4, 15};
Periodic Surface surfSlave { boundSlave[] } = surfMaster { boundMaster[] };

// // Physical volumes and surfaces
Physical Volume(1) = {cylinder};
//Physical Volume(2) = {box};

interfaces[] = Unique(Abs(Boundary{ Volume{cylinder}; }));  
boundary[] = Unique(Abs(Boundary{ Volume{box}; }));
boundary[] -= {interfaces[]};

//Physical Surface(2) = {boundary[]};
Physical Surface(1) = {interfaces[]};