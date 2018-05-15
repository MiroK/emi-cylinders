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
min_y = -depth - length_y - pady;
min_z = -radius - padz;

max_x = -length/2 - length_x + nx*x_shift + padx;
max_y = -depth - length_y + ny*y_shift + pady;
max_z = radius + padz;

// Exterior
v = newv;
Box(v) = {min_x, min_y, min_z, max_x-min_x, max_y-min_y, max_z-min_z};
BooleanFragments {Volume{v}; Delete; }{Volume{volumes[]}; Delete; }

// COlro each cylinder individually
For row In {1:ny}
  For col In {1:nx}
    v = (row-1)*nx + col;
    Physical Volume(v) = {v};
  EndFor
EndFor
// External is largest
ncells = nx*ny;
box = ncells + 1;
Physical Volume(box) = {box};

ce_interfaces[] = {};  // Cell and exterior
cc_interfaces[] = {};
// Mark surfaces of volumes
For row In {1:ny}
  For col In {1:nx}
    v = (row-1)*nx + col;
    
    v_surfaces() = Unique(Abs(Boundary{ Volume{v}; }));
    vcount = #v_surfaces[];
    shared[] = {};
    
    vleft = v+1;
    vright = v-1;
    vup = v + nx;
    vdown = v -nx;
    
    neighbors[] = {vleft, vright, vup, vdown};
    For n In {0:3}
      neighbor = neighbors[n];
      If(1 <= neighbor && neighbor <= ncells)
        neighbor_surfaces() = Unique(Abs(Boundary{ Volume{neighbor}; }));
        // Check for intersection
      
        ncount = #neighbor_surfaces[];
        For i In {0:(ncount-1)}
          ai = neighbor_surfaces[i];
          For j In {0:(vcount-1)}
            bj = v_surfaces[j];
            If (ai == bj)
              shared[] += {bj};
              cc_interfaces[] += {bj};
              
              If (v > neighbor)
                x = Round(Log10(neighbor)) + 1;
                y = v;
                tens = Round(Exp(x*Log(10)));
                tag = tens*neighbor + v;
                Physical Surface(tag) = {ai};
              EndIf
            EndIf
          EndFor
        EndFor
      EndIf
    EndFor
    v_surfaces[] -= {shared[]};
    ce_interfaces[] += {v_surfaces[]};
  EndFor
EndFor

Physical Surface(1) = {ce_interfaces[]};

// Let's also mark the bounding box; it is a difference of the bounding
// surfaces of cells
all[] = Unique(Abs(Boundary{ Volume{box}; }));
all[] -= {cc_interfaces[]};
all[] -= {ce_interfaces[]};
// Mark it as 2
Physical Surface(2) = {all[]};

// Characteristic Length{interfaces[]} = size_cell;
// Characteristic Length{cc_interfaces[]} = size_cell;
// Characteristic Length{all[]} = size_box;