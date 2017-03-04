z_cylinders = 4;
y_cylinders = 3;
// The height of each cylinder is H with char mesh size size_c
R = 0.8;
H = 1.5;
size_R = 0.2;
// The height of the joint regions in z-direction is h with char mesh size size_r
r = 0.5;
h = 0.5;
size_r = 0.2;
// The height of the joint regions in y-direction is l with char mesh size size_l
t = 0.25;
l = 0.5;
size_l = 0.2;
// Note that to make the transition of the joint to the cylindrical surface
// there is a weld which results in effective radius of q > t. Below q is the
// scaling factor so that q[radius] is q*t
q = 1.1;
// The cylinders are enclosed in a bbox which leaves dx, dy gaps around the
// square that bounds the cylinder crossection. In z direction the gap is dz.
// For bbox the char size is size_b
dx = 0.2;
dy = 0.2;
dz = 0.2;
size_b = 0.3;
// The cylinders are numbered (the volume) starting from on first growing in 
// y-direction then in z-direction. The volume id determines also id of all the
// surfaces that are not shared. For the shared surfaces the y-shared surface is
// volume_id + (z_cylinders*y_cylinders). In z direction this is 
// volume_id + 2*(z_cylinders*y_cylinders). Note that in this policy it is the first
// cylinder which hits the surface that gets to name it.
// The outer volume, i.e. bounding box - cylinder union is tagged as 0
//! Cut here

// The first cylinder has center at x0, y0, z0
// NOTE: all of this is written assuming x0, y0, z0 = origin so that's why it is
// not a parameter
x0 = 0;
y0 = 0;
z0 = 0;
// Mesh quality/size parameters
Mesh.Smoothing = 4;
Mesh.SmoothNormals = 4;

// Let's compute stuff for the y-joint
q = q*t;            // Control the curvature of the joint part
sA = q/R;
cA = Sqrt(R^2-q^2)/R;
tA = sA/cA;
dp = tA*(q-t);
d = (q-t)/cA;
A = Asin(q/R);
X = d*Tan((Pi/2+A)/2);
Cx = cA*(R+X);
Cy = sA*(R+X);

// The first cylinder is drawn by hand using a lot of symmetries
Point(1) = {0, 0, H/2+h/2, size_r};
Point(2) = {r, 0, H/2+h/2, size_r};
Point(3) = {q, Sqrt(r^2-q^2), H/2+h/2, size_r};
Point(4) = {0, r, H/2+h/2, size_r};

Point(5) = {0, 0, H/2, size_R};
Point(6) = {r, 0, H/2, size_R};
Point(7) = {q, Sqrt(r^2-q^2), H/2, size_R};
Point(8) = {0, r, H/2, size_R};

Point(9) = {R, 0, H/2, size_R};
Point(10) = {q, Sqrt(R^2-q^2), H/2, size_R};
Point(11) = {0, R, H/2, size_R};

Point(12) = {R, 0, 0, size_R};
Point(13) = {q, Sqrt(R^2-q^2), 0, size_R};

Point(14) = {0, R+l/2, 0, size_l};
Point(15) = {0, R+l/2, t, size_l};
Point(16) = {t, R+l/2, 0, size_l};
Point(17) = {0, R+(q-t), t, size_l};

Point(18) = {0, 0, 0, size_R};
Point(19) = {0, R, t, size_R};
Point(20) = {0, R, q, size_R};
Point(21) = {0, R, q, size_R};
Point(22) = {0, R+(q-t), q, size_R};

Point(23) = {t, Sqrt(R^2-t^2), 0, size_R};
Point(25) = {t, Sqrt(R^2-q^2)+dp+d, 0, size_R};
Point(26) = {Cy, Cx, 0, size_R};
Point(27) = {0, R+(q-t), 0, size_R};
Point(28) = {0, R, 0, size_R};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {6, 5, 7};
Circle(4) = {7, 5, 8};
Line(5) = {4, 8};
Line(6) = {3, 7};
Line(7) = {2, 6};
Circle(8) = {9, 5, 10};
Circle(9) = {10, 5, 11};
Circle(10) = {12, 18, 13};
Line(11) = {6, 9};
Line(12) = {7, 10};
Line(13) = {8, 11};
Line(14) = {9, 12};
Line(15) = {10, 13};
Circle(16) = {15, 14, 16};
Circle(17) = {20, 22, 17};
Line(18) = {11, 20};
Line(19) = {17, 15};
Circle(20) = {13, 26, 25};
Line(21) = {16, 25};

Ellipse(22) = {17, 27, 17, 25};
Ellipse(23) = {20, 28, 20, 13};
Line Loop(24) = {19, 16, 21, -22};
Ruled Surface(25) = {24};
Line Loop(26) = {23, 20, -22, -17};
Ruled Surface(27) = {26};
Line Loop(28) = {18, 23, -15, 9};
Ruled Surface(29) = {28};
Line Loop(30) = {8, 15, -10, -14};
Ruled Surface(31) = {30};
Line Loop(32) = {1, 6, -3, -7};
Ruled Surface(33) = {32};
Line Loop(34) = {2, 5, -4, -6};
Ruled Surface(35) = {34};
Line Loop(36) = {11, 8, -12, -3};
Plane Surface(37) = {36};
Line Loop(38) = {4, 13, -9, -12};
Plane Surface(39) = {38};
// Now this is 1/8 of the shape

eight[] = {25, 27, 29, 31, 37, 39, 33, 35};
mirror_z[] = Symmetry {0, 0, 1, 0} { Duplicata { Surface{eight[]}; } };

quarter[] = eight[];
quarter[] += mirror_z[];
mirror_y[] = Symmetry {0, 1, 0, 0} { Duplicata { Surface{quarter[]}; } };

half[] = quarter[];
half[] += mirror_y[];
mirror_x[] = Symmetry {1, 0, 0, 0} { Duplicata { Surface{half[]}; } };

root[] = half[];
root[] += mirror_x[];
// Let's get the closing surfaces for y
Line Loop(313) = {80, -120, 278, -238};
Plane Surface(314) = {313};
Line Loop(315) = {158, -198, 42, -16};
Plane Surface(316) = {315};
// Z are a bit more involved

Line(317) = {1015, 144};
Line(318) = {144, 430};
Line(319) = {144, 715};
Line(320) = {144, 145};
Line Loop(321) = {317, 319, -227, 307};
Plane Surface(322) = {321};
Line Loop(323) = {317, 318, 154, -312};
Plane Surface(324) = {323};
Line Loop(325) = {318, -149, 71, -320};
Plane Surface(326) = {325};
Line Loop(327) = {320, 76, -232, -319};
Plane Surface(328) = {327};

Line(329) = {280, 1};
Line(330) = {1, 865};
Line(331) = {1, 565};
Line(332) = {1, 3};
Line Loop(333) = {114, -272, -330, -329};
Plane Surface(334) = {333};
Line Loop(335) = {329, 332, -1, 109};
Plane Surface(336) = {335};
Line Loop(337) = {332, 2, -192, -331};
Plane Surface(338) = {337};
Line Loop(339) = {330, -267, 187, -331};
Plane Surface(340) = {339};

cylinder[] = {314, 328, 326, 324, 322};
cylinder[] += {316, 340, 338, 336, 334};
cylinder[] += root[];

sl = newsl;
Surface Loop(sl) = cylinder[];
v = newv;
Volume(v) = {sl};
pv = 1;
Physical Volume(pv) = {v};
Physical Surface(pv) = root[]; // Don't forget to add caps
Physical Surface(pv) += {328, 326, 324, 322}; // bottom
Physical Surface(pv) += {314}; // left
// Finally claim the shared y
n_cylinders = y_cylinders*z_cylinders;
Physical Surface(pv+n_cylinders) = {316};
// And z shared
Physical Surface(pv+2*n_cylinders) = {340, 338, 336, 334};

// Collect cylinders that we have defined
all_cylinders[] = {v};

parent = v;
For i In {1:z_cylinders}
  // Mirror in y
  last_volume = parent;
  For j In {2:y_cylinders}
    new_volume[] = Translate {0, 2*(R+l/2), 0} { Duplicata { Volume{last_volume}; } };

    last_volume = new_volume[0];
    all_cylinders[] += {last_volume};
    pv += 1;
    Physical Volume(pv) = {last_volume};

    // Extract shell
    b() = Boundary{ Volume{last_volume}; } ;
    cylinder_shell[] = {};
    For k In {10:73}
        cylinder_shell[] += {b[k]};
    EndFor
    Physical Surface(pv) = cylinder_shell[];
    // Bottom claims
    If(i==1)
      Physical Surface(pv) += {b[1], b[2], b[3], b[4]};
    EndIf
    // Top claims
    If(i == z_cylinders)
       Physical Surface(pv) += {b[6], b[7], b[8], b[9]};
    EndIf
    // Right claims
    If(j==y_cylinders)
        Physical Surface(pv) += {b[5]};
    EndIf
    // Claim y-shared
    If(j < y_cylinders)
        Physical Surface(pv+n_cylinders) = {b[5]};
    EndIf
    // Claim z-shared
    If(i < z_cylinders)
        Physical Surface(pv+2*n_cylinders) = {b[6], b[7], b[8], b[9]};
    EndIf
  EndFor 

  If(i < z_cylinders)
    new_volume[] = Translate {0, 0, H+h} { Duplicata { Volume{parent}; } };
    parent = new_volume[0];
    pv += 1;
    Physical Volume(pv) = {parent};
    all_cylinders[] += {parent};

    // Extract shell
    b() = Boundary{ Volume{parent}; } ;
    cylinder_shell[] = {};
    For k In {10:73}
        cylinder_shell[] += {b[k]};
    EndFor
    Physical Surface(pv) = cylinder_shell[];
    // Left claims
    Physical Surface(pv) += {b[0]};
    // Claim y-shared
    Physical Surface(pv+n_cylinders) = {b[5]};
    // Top claims
    If((i+1)==z_cylinders)
       Physical Surface(pv) += {b[6], b[7], b[8], b[9]};
    EndIf
    If((i+1)<z_cylinders)
        Physical Surface(pv+2*n_cylinders) = {b[6], b[7], b[8], b[9]};
    EndIf
  EndIf

EndFor

// Get the bounding surfaces for cylinder
bbox_cylinders = CombinedBoundary{ Volume{all_cylinders[]}; } ;

bbox_surfaces[] = {};
///////////////////////////////////////////////////////////////////////////////
// Bounding box
///////////////////////////////////////////////////////////////////////////////
// Bottom corner
p = newp;
dx = dx + R+l/2;
ddy = dy + R+l/2;
ddY = dy + 2*(R+l/2)*(y_cylinders-1) + R+l/2;
Point(p) =   {x0+dx, y0+ddY, z0-H/2-h/2-dz, size_b};
Point(p+1) = {x0+dx, y0-ddy, z0-H/2-h/2-dz, size_b};
Point(p+2) = {x0-dx, y0-ddy, z0-H/2-h/2-dz, size_b};
Point(p+3) = {x0-dx, y0+ddY, z0-H/2-h/2-dz, size_b};
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

dz = dz + (H+h)*(z_cylinders-1) + H/2 + h/2;
// Top corner
P = newp;
Point(P) =   {x0+dx, y0+ddY, z0+dz, size_b};
Point(P+1) = {x0+dx, y0-ddy, z0+dz, size_b};
Point(P+2) = {x0-dx, y0-ddy, z0+dz, size_b};
Point(P+3) = {x0-dx, y0+ddY, z0+dz, size_b};
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
//
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
//v = newv;
Surface Loop(v) = {bbox_surfaces[]};
Surface Loop(v+1) = {bbox_cylinders[]};
Volume(v+2) = {v, v+1};
Physical Volume(0) = {v+2};
