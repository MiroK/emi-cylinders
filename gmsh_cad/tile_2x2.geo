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

v = newv;
// The cylinder
Cylinder(v) = {0-length/2-length_x, 0, 0, length+2*length_x, 0, 0, radius};
Cylinder(v+1) = {0, -depth-length_y, 0, 0, y_shift, 0, radius_y};
cylinder() = BooleanUnion{ Volume{v}; Delete;}{ Volume{v, v+1}; Delete;};


c_right() = Translate {x_shift, 0, 0}{ Duplicata{ Volume{cylinder[]}; } };
c_up() = Translate {0, y_shift, 0}{ Duplicata{ Volume{cylinder[]}; } };
c_up_right() = Translate {x_shift, y_shift, 0}{ Duplicata{ Volume{cylinder[]}; } };
cells() = BooleanFragments {Volume{cylinder[]}; Delete; }{Volume{c_right[], c_up[], c_up_right[]}; Delete; };

min_x = -length/2 - length_x;
min_y = -depth - length_y;
min_z = -radius - padz;

max_x = min_x + 2*x_shift;
max_y = min_y + 2*y_shift;
max_z = radius + padz;

box = newv;
Box(box) = {min_x, min_y, min_z, max_x-min_x, max_y-min_y, max_z-min_z}; 

v() = BooleanFragments {Volume{box}; Delete; }{Volume{cells()}; Delete; };

cylinders[] = {v[0], v[1], v[2], v[3]};
exteriors[] = {v[4]};

// X periodicity
surfMaster = 35;
surfSlave = 39;
boundMaster[] = {25};
boundSlave[] = {45};
Periodic Surface surfSlave { boundSlave[] } = surfMaster { boundMaster[] };

surfMaster = 37;
surfSlave = 41;
boundMaster[] = {26};
boundSlave[] = {46};
Periodic Surface surfSlave { boundSlave[] } = surfMaster { boundMaster[] };

surfMaster = 17;
surfSlave = 24;
boundMaster[] = {21, 23, 24, 22, 25, 26};
boundSlave[] = {28, 35, 36, 33, 45, 46};
Periodic Surface surfSlave { boundSlave[] } = surfMaster { boundMaster[] };

// Y periodicity
surfMaster = 38;
surfSlave = 36;
boundMaster[] = {37};
boundSlave[] = {30};
Periodic Surface surfSlave { boundSlave[] } = surfMaster { boundMaster[] };

surfMaster = 42;
surfSlave = 40;
boundMaster[] = {38};
boundSlave[] = {31};
Periodic Surface surfSlave { boundSlave[] } = surfMaster { boundMaster[] };

surfMaster = 21;
surfSlave = 18;
boundMaster[] = {24, 34, 36, 32, 37, 38};
boundSlave[] = {21, 27, 28, 29, 30, 31};
Periodic Surface surfSlave { boundSlave[] } = surfMaster { boundMaster[] };

Physical Volume(1) = {cylinders[]};
Physical Volume(2) = {exteriors[]};

interfaces[] = Unique(Abs(Boundary{ Volume{cylinders[]}; }));  
boundary[] = Unique(Abs(Boundary{ Volume{exteriors[]}; }));
boundary[] -= {interfaces[]};

Physical Surface(2) = {boundary[]};
Physical Surface(1) = {interfaces[]};