DefineConstant[
radius = {10, Name "radius of cell"}
radius_x = {8, Name "radius of connection in x direction"}
radius_y = {6, Name "radius of connection in y direction"}
length = {100, Name "length of cell (body)"}
length_x = {4, Name "length of connection in x direction"}
length_y = {10, Name "length of connection in y direction"}
nx = {2, Name "number of cells in x direction"}
ny = {2, Name "number o cells in y direction"}
padx = {50, Name "bounding box padding in x direction"}
pady = {50, Name "bounding box padding in y direction"}
padz = {50, Name "bounding box padding in z direction"}
size_cell = {4, Name "mesh size on cell surface"}
size_box = {4, Name "mesh size on bounding box surface"}
];
// The length units here micro meters
SetFactory("OpenCASCADE");

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

// Generate the x arrays
row[] = {cylinder};

For i In {1:(nx-1)}
  next = Translate {x_shift, 0, 0}{ 
    Duplicata{ Volume{row[i-1]}; }
  };
  row[] += {next};
EndFor

volumes[] = {row[]};  // interior
// Shift rows to get the tailing in y
For j In {1:(ny-1)}
  row[] = Translate {0, y_shift, 0}{ 
    Duplicata{ Volume{row[]}; }
  };
  volumes[] += {row[]};
EndFor

// Define bounding box
min_x = -length/2 - length_x - padx;
min_y = -depth - length_y;
min_z = -radius - padz;

max_x = -length/2 - length_x + nx*x_shift;
max_y = -depth - length_y + ny*y_shift;
max_z = radius + padz;

// Exterior
v = newv;
Box(v) = {min_x, min_y, min_z, max_x-min_x, max_y-min_y, max_z-min_z};

cell_bdry() = Unique(Abs(Boundary{ Volume{volumes[]}; }));
box_bdry() = Unique(Abs(Boundary{ Volume{v}; }));
// Make sure that the mesh conforms to circles
BooleanFragments {Surface{box_bdry[]}; Delete; }{Surface{cell_bdry[]}; Delete; }

BooleanFragments {Volume{v}; Delete; }{Volume{volumes[]}; Delete; }

// Mark interior and exterior domain differently to set material parameters
// 1 is for inside
Physical Volume(1) = {volumes[]};
// 2 is out
box = nx*ny + 1;
Physical Volume(2) = {box};

exterior_surface[] = Unique(Abs(Boundary{ Volume{box}; }));
next = #exterior_surface[];

cc_interface[] = {};
ce_interface[] = {};

For v In {1:(nx*ny)}
  volume_surface[] = Unique(Abs(Boundary{ Volume{v}; }));
  a = #volume_surface[];
    
  this_ce[] = {};
  // Each volume surface is either an cell-exterior type
  For j In {0:(next-1)}
    bj = exterior_surface[j];
    For i In {0:(a-1)}
      ai = volume_surface[i];
      If (ai == bj)
        this_ce[] += {bj};
      EndIf
    EndFor
  EndFor
  ce_interface[] += {this_ce[]};
  volume_surface[] -= {this_ce[]};
    
  // The remaining surfaces are 3
  cc_interface[] += {volume_surface[]};
EndFor

// Interfaces between cells and exterior
Physical Surface(1) = {ce_interface[]};
// Interfaces between cells 
Physical Surface(3) = {cc_interface[]};

box_bdry[] -= {ce_interface[]};

Physical Surface(2) = {box_bdry[]};

Characteristic Length{ce_interface[]} = size_cell;
Characteristic Length{cc_interface[]} = size_cell;
Characteristic Length{box_bdry[]} = size_box;