// A rectangular domain is tiles by cells which are of the form
//       ___
//   ___|   |____
//  |___     ____|
//      |___|
//  

DefineConstant[
length = {100, Name "length of cell (body)"}
length_x = {40, Name "length of connection in x direction"}
length_y = {40, Name "length of connection in y direction"}
padx = {40, Name "bounding box padding in x direction"}
pady = {40, Name "bounding box padding in y direction"}
nx = {2, Name "number of cells in x direction"}
ny = {2, Name "number o cells in y direction"}
size_cell = {60, Name "mesh size on cell surface"}
size_box = {60, Name "mesh size on bounding box surface"}
];
// The length units here micro meters

//////////////////////////////////////////////////////////////////////

x_shift = length + 2*length_x;
y_shift = length + 2*length_y;

Macro Tile
  cp[] = {};
  p = newp;
  For i In {0:19}
      cp[i] = p;
      p += 1;
  EndFor
  
  Point(cp[0]) = {-length/2-length_x, length_y/2, 0, size_cell};
  Point(cp[1]) = {-length/2, length_y/2, 0, size_cell};
  Point(cp[2]) = {-length/2, length/2, 0, size_cell};
  Point(cp[3]) = {-length_x/2, length/2, 0, size_cell};
  Point(cp[4]) = {-length_x/2, length/2+length_y, 0, size_cell};
  Point(cp[5]) = {length_x/2, length/2+length_y, 0, size_cell};
  Point(cp[6]) = {length_x/2, length/2, 0, size_cell};
  Point(cp[7]) = {length/2, length/2, 0, size_cell};
  Point(cp[8]) = {length/2, length_y/2, 0, size_cell};
  Point(cp[9]) = {length/2+length_x, length_y/2, 0, size_cell};
  //
  Point(cp[19]) = {-length/2-length_x, -length_y/2, 0, size_cell};
  Point(cp[18]) = {-length/2, -length_y/2, 0, size_cell};
  Point(cp[17]) = {-length/2, -length/2, 0, size_cell};
  Point(cp[16]) = {-length_x/2, -length/2, 0, size_cell};
  Point(cp[15]) = {-length_x/2, -length/2-length_y, 0, size_cell};
  Point(cp[14]) = {length_x/2,-length/2-length_y, 0, size_cell};
  Point(cp[13]) = {length_x/2, -length/2, 0, size_cell};
  Point(cp[12]) = {length/2, -length/2, 0, size_cell};
  Point(cp[11]) = {length/2, -length_y/2, 0, size_cell};
  Point(cp[10]) = {length/2+length_x, -length_y/2, 0, size_cell};
  
  // Trace the boundary
  l = newl;
  cl[] = {};
  For i In {0:18}
    Line(l) = {cp[i], cp[i+1]};
    
    cl[] += {l};
    l += 1;
  EndFor
  Line(l) = {cp[19], cp[0]};
  cl[] += l;
  
  cell_loop = newll;
  Line Loop(cell_loop) = {
    cl[0], cl[1], cl[2], cl[3], cl[4],
    cl[5], cl[6], cl[7], cl[8], cl[9],
    cl[10], cl[11], cl[12], cl[13], cl[14],
    cl[15], cl[16], cl[17], cl[18], cl[19]};
  
  cell_surface = news;
  Plane Surface(cell_surface) = {cell_loop};
  
  // Bounding
  min_x = -length/2 - length_x;
  min_y = -length/2 - length_y;

  max_x = length/2 + length_x;
  max_y = length/2 + length_y;

  p = newp;
  Point(p) = {min_x, min_y, 0, size_box};
  Point(p+1) = {max_x, min_y, 0, size_box};
  Point(p+2) = {max_x, max_y, 0, size_box};
  Point(p+3) = {min_x, max_y, 0, size_box};

  l = newl;
  Line(l) = {cl[19], p};
  Line(l+1) = {p, cl[15]};
  Line(l+2) = {cl[14], p+1};
  Line(l+3) = {p+1, cl[10]};
  Line(l+4) = {cl[9], p+2};
  Line(l+5) = {p+2, cl[5]};
  Line(l+6) = {cl[4], p+3};
  Line(l+7) = {p+3, cl[0]};
  
  bl[] = {l, l+1, l+2, l+3, l+4, l+5, l+6, l+7};
  // We have 4 exterior loops
  ext_cl1 = newll;
  Line Loop(ext_cl1) = {bl[2], bl[3], cl[10], cl[11], cl[12], cl[13]};
  ext_s1 = news;
  Plane Surface(ext_s1) = {ext_cl1};
  
  ext_cl2 = newll;
  Line Loop(ext_cl2) = {bl[4], bl[5], cl[5], cl[6], cl[7], cl[8]};
  ext_s2 = news;
  Plane Surface(ext_s2) = {ext_cl2};

  ext_cl3 = newll;
  Line Loop(ext_cl3) = {bl[6], bl[7], cl[0], cl[1], cl[2], cl[3]};
  ext_s3 = news;
  Plane Surface(ext_s3) = {ext_cl3};

  ext_cl4 = newll;
  Line Loop(ext_cl4) = {bl[0], bl[1], cl[15], cl[16], cl[17], cl[18]};
  ext_s4 = news;
  Plane Surface(ext_s4) = {ext_cl4};
  
  exterior_surfaces[] = {ext_s1, ext_s2, ext_s3, ext_s4};
Return 

// Tile it
Call Tile;
// First cells
// Generate the x arrays
crow[] = {cell_surface};
For i In {1:(nx-1)}
  next = Translate {x_shift, 0, 0}{ 
    Duplicata{ Surface{crow[i-1]}; }
  };
  crow[] += {next};
EndFor

cells[] = {crow[]};  // interior
// Shift rows to get the tailing in y
For j In {1:(ny-1)}
  crow[] = Translate {0, y_shift, 0}{ 
    Duplicata{ Surface{crow[]}; }
  };
  cells[] += {crow[]};
EndFor

// Now exteriors
// Generate the x arrays
crow[] = {exterior_surfaces[]};
next[] = {exterior_surfaces[]};
For i In {1:(nx-1)}
  next[] = Translate {x_shift, 0, 0}{ 
    Duplicata{ Surface{next[]}; }
  };
  crow[] += {next[]};
EndFor

exteriors[] = {crow[]};  // interior
// Shift rows to get the tailing in y
For j In {1:(ny-1)}
  crow[] = Translate {0, y_shift, 0}{ 
    Duplicata{ Surface{crow[]}; }
  };
  exteriors[] += {crow[]};
EndFor

domains[] = {cells[]};
domains[] += {exteriors[]};
// Doing all the Abs messes up the orientation!!!!
// I know the type of edges that are wrong
shared[] = {};
n = #domains[];
For d In {0:(n-1)}
  domain_0 = domains[d];
  surfaces_0 = Unique(Abs(Boundary{ Surface{domain_0}; }));
  ns0 = #surfaces_0[];  
  For D In {(d+1):(n-1)}
    domain_1 = domains[D];      
    surfaces_1 = Unique(Abs(Boundary{ Surface{domain_1}; }));
    ns1 = #surfaces_1[];
    
    For i In {0:(ns0-1)}
      bdry = surfaces_0[i];
      For j In {0:(ns1-1)}
        match = surfaces_1[j];
        If(bdry == match)
          shared[] += {bdry};
        EndIf
      EndFor
    EndFor
    
  EndFor
EndFor

interfaces[] = Unique(Boundary{ Surface{domains[]}; });
shared[] = Unique(shared[]);

ni = #interfaces[];
ns = #shared[];

not_shared[] = {};
For i In {0:(ni-1)}
  ei = interfaces[i];
  flag = 0;
  For s In {0:(ns-1)}
    es = shared[s];
   
    If(es == ei || es == -ei)
      flag += 1;
    EndIf
  EndFor
  
  If(flag == 0)
    not_shared[] += {ei};
  EndIf
EndFor
 
// The frame
min_x = -length/2 - length_x - padx;
min_y = -length/2 - length_y - pady;

max_x = min_x + nx*x_shift + 2*padx;
max_y = min_y + ny*y_shift + 2*pady;

p = newp;
Point(p) = {min_x, min_y, 0, size_box};
Point(p+1) = {max_x, min_y, 0, size_box};
Point(p+2) = {max_x, max_y, 0, size_box};
Point(p+3) = {min_x, max_y, 0, size_box};

l = newl;
Line(l) = {p, p+1};
Line(l+1) = {p+1, p+2};
Line(l+2) = {p+2, p+3};
Line(l+3) = {p+3, p};

frame[] = {l, l+1, l+2, l+3};
frame_ll = newl;
Line Loop(frame_ll) = {frame[]};

// Go every boundary volume 
// then its every surface that is in shared should be flipped
x = #not_shared[];
For i In {0:(x-1)}
  edge = not_shared[i];
 
  found = 0;
  For cell In {0:(nx*ny-1)}
    bdry() = Unique(Boundary{ Surface{cells[cell]}; });
    y = #bdry[];
    
    For j In {0:(y-1)}
      match = bdry[j];
      If(match == edge || match == -edge)
        found += 1;
      EndIf
    EndFor
  EndFor
  
  If(found > 0)
    not_shared[i] = -edge;
  EndIf
EndFor

inner_ll = newl;
Line Loop(inner_ll) = {not_shared[]};

s = news;
Plane Surface(s) = {inner_ll, frame_ll};

// Mark volumes
//
Physical Surface(1) = {cells[]};  // cell are 1
Physical Surface(2) = {exteriors[]};  // rim and exterior are 2
Physical Surface(2) += {s};

cell_interfaces[] = Unique(Abs(Boundary{ Surface{cells[]}; }));
Physical Curve(1) = {cell_interfaces[]};

Physical Curve(2) = {frame[]};
