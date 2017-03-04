x_cylinders = 3;
y_cylinders = 4;
// The height of each cylinder is H with char mesh size size_R
R = 0.8;
H = 1.5;
size_R = 0.2;
// The length of the joint regions in x/y-direction is h/2 with char mesh size size_r
r = 0.2;
h = 0.5;
size_r = 0.2;
// Note that to make the transition of the joint to the cylindrical surface
// there is a weld which results in effective radius of q > t. Below q is the
// scaling factor so that q[radius] is q*t
q = 1.1;
// The cylinders are enclosed in a bbox which leaves dx, gaps around the
// square that bounds the cylinder crossection. In z direction the gap is dz.
// For bbox the char size is size_b
dx = 0.2;
dz = 0.2;
size_b = 0.3;
// The cylinders are numbered (the volume) starting from on first growing in 
// x-direction then in y-direction. The volume id determines also id of all the
// surfaces that are not shared. For the shared surfaces the x-shared surface is
// volume_id + (y_cylinders*x_cylinders). In y direction this is 
// volume_id + 2*(x_cylinders*y_cylinders). Note that in this policy it is the first
// cylinder which hits the surface that gets to name it.
// The outer volume, i.e. bounding box - cylinder union is tagged as 0
//! Cut here

// Mesh quality/size parameters
Mesh.Smoothing = 4;
Mesh.SmoothNormals = 4;

// Let's compute stuff for the joints
q = q*r;            // Control the curvature of the joint part
sA = q/R;
cA = Sqrt(R^2-q^2)/R;
tA = sA/cA;
dp = tA*(q-r);
d = (q-r)/cA;
A = Asin(q/R);
X = d*Tan((Pi/2+A)/2);
Cx = cA*(R+X);
Cy = sA*(R+X);

// The first cylinder is drawn by hand using a lot of symmetries
fact = Sqrt(2)/2;
Point(1) = {R*fact, R*fact, H/2, size_R};
Point(2) = {R*fact, R*fact, 0, size_R};
Point(3) = {R+h/2, 0, r, size_r};
Point(4) = {R+h/2, 0, 0, size_r};
Point(5) = {R+h/2, r, 0, size_r};
Point(6) = {R, 0, H/2, size_R};
Point(7) = {R, 0, q, size_r};
Point(8) = {R+(q-r), 0, q, size_r};
Point(9) = {R, 0, r, size_R};
Point(10) = {R+(q-r), 0, r, size_r};
Point(11) = {R, 0, 0, size_R};
Point(12) = {Sqrt(R^2-q^2), q, 0, size_r};
Point(13) = {Sqrt(R^2-q^2), q, H/2, size_r};
Point(14) = {Sqrt(R^2-q^2)+dp+d, r, 0, size_r};
Point(15) = {Cx, Cy, 0, size_R};
Point(16) = {Sqrt(R^2-r^2), r, 0, size_r};
Point(17) = {0, 0, H/2, size_R};
Point(18) = {0, 0, 0, size_R};
Circle(1) = {13, 17, 1};
Circle(2) = {12, 18, 2};
Line(3) = {1, 2};
Line(4) = {12, 13};
Circle(5) = {6, 17, 13};
Circle(6) = {3, 4, 5};
Line(7) = {3, 10};
Line(8) = {6, 7};
Circle(9) = {7, 8, 10};
Line(10) = {5, 14};
Circle(11) = {14, 15, 12};
Ellipse(12) = {7, 11, 7, 12};
Ellipse(13) = {10, 11, 10, 14};
Line Loop(14) = {7, 13, -10, -6};
Ruled Surface(15) = {14};
Line Loop(16) = {12, -11, -13, -9};
Ruled Surface(17) = {16};
Line Loop(18) = {8, 12, 4, -5};
Ruled Surface(19) = {18};
Line Loop(20) = {4, 1, 3, -2};
Ruled Surface(21) = {20};

// The rest is mirroring
sixteenth[] = {15, 17, 19, 21};
mirror_z[] = Symmetry {0, 0, 1, 0} { Duplicata { Surface{sixteenth[]}; } };

eigth[] = sixteenth[];
eigth[] += mirror_z[];
mirror_y[] = Symmetry {0, 1, 0, 0} { Duplicata { Surface{eigth[]}; } };

x_pos[] = eigth[];
x_pos[] += mirror_y[];

y_neg[] = Symmetry {1, 1, 0, 0} { Duplicata { Surface{x_pos[]}; } };
y_pos[] = Symmetry {1, -1, 0, 0} { Duplicata { Surface{x_pos[]}; } };

x_neg[] = Symmetry {-1, -1, 0, 0} { Duplicata { Surface{y_pos[]}; } };

shell[] = x_pos[];
shell[] += y_neg[];
shell[] += y_pos[];
shell[] += x_neg[];

// Let see it from the top and then bottom
Line(315) = {297, 17};
Line(316) = {871, 17};
Line(317) = {17, 860};
Line(318) = {17, 1027};
Line(319) = {740, 17};
Line(320) = {17, 726};
Line(321) = {17, 559};
Line(322) = {570, 17};
Line(323) = {17, 1};
Line(324) = {13, 17};
Line(325) = {17, 6};
Line(326) = {138, 17};
Line(327) = {152, 17};
Line(328) = {17, 439};
Line(329) = {272, 17};
Line(330) = {17, 283};
Line Loop(331) = {316, -315, -254};
Plane Surface(332) = {331};
Line Loop(333) = {315, 330, 97};
Plane Surface(334) = {333};
Line Loop(335) = {330, 94, 329};
Plane Surface(336) = {335};
Line Loop(337) = {329, 328, 134};
Plane Surface(338) = {337};
Line Loop(339) = {328, 137, 327};
Plane Surface(340) = {339};
Line Loop(341) = {327, -326, 58};
Plane Surface(342) = {341};
Line Loop(343) = {326, 325, -55};
Plane Surface(344) = {343};
Line Loop(345) = {325, 5, 324};
Plane Surface(346) = {345};
Line Loop(347) = {324, 323, -1};
Plane Surface(348) = {347};
Line Loop(349) = {323, -175, 322};
Plane Surface(350) = {349};
Line Loop(351) = {322, 321, -172};
Plane Surface(352) = {351};
Line Loop(353) = {321, -212, -320};
Plane Surface(354) = {353};
Line Loop(355) = {320, 215, 319};
Plane Surface(356) = {355};
Line Loop(357) = {319, 318, 294};
Plane Surface(358) = {357};
Line Loop(359) = {318, 291, -317};
Plane Surface(360) = {359};
Line Loop(361) = {317, -251, 316};
Plane Surface(362) = {361};

top_cap[] = {332, 334, 336, 338, 340, 342, 344, 346, 348, 350, 352, 354, 356, 358, 360, 362};

Line(363) = {949, 75};
Line(364) = {375, 75};
Line(365) = {75, 361};
Line(366) = {75, 350};
Line(367) = {75, 517};
Line(368) = {230, 75};
Line(369) = {75, 216};
Line(370) = {60, 75};
Line(371) = {75, 71};
Line(372) = {85, 75};
Line(373) = {75, 648};
Line(374) = {637, 75};
Line(375) = {75, 804};
Line(376) = {75, 818};
Line(377) = {1105, 75};
Line(378) = {75, 938};
Line Loop(379) = {364, -363, 274};
Plane Surface(380) = {379};
Line Loop(381) = {363, 378, -271};
Plane Surface(382) = {381};
Line Loop(383) = {378, -311, 377};
Plane Surface(384) = {383};
Line Loop(385) = {377, 376, -314};
Plane Surface(386) = {385};
Line Loop(387) = {376, -235, -375};
Plane Surface(388) = {387};
Line Loop(389) = {375, 232, 374};
Plane Surface(390) = {389};
Line Loop(391) = {374, 373, 192};
Plane Surface(392) = {391};
Line Loop(393) = {373, 195, 372};
Plane Surface(394) = {393};
Line Loop(395) = {372, 371, 39};
Plane Surface(396) = {395};
Line Loop(397) = {371, 36, 370};
Plane Surface(398) = {397};
Line Loop(399) = {370, 369, 75};
Plane Surface(400) = {399};
Line Loop(401) = {369, 78, 368};
Plane Surface(402) = {401};
Line Loop(403) = {368, 367, 157};
Plane Surface(404) = {403};
Line Loop(405) = {367, 154, -366};
Plane Surface(406) = {405};
Line Loop(407) = {364, 365, 117};
Plane Surface(408) = {407};
Line Loop(409) = {365, 114, -366};
Plane Surface(410) = {409};


bot_cap[] = {380, 382, 384, 386, 388, 390, 392, 394, 396, 398, 400, 402, 404, 406, 408, 410};

// Finally surfaces on joint regions
Line Loop(411) = {281, -241, 261, -301};
Plane Surface(412) = {411};
x_neg = 412;

Line Loop(413) = {124, -84, 104, -144};
Plane Surface(414) = {413};
y_neg = 414;

Line Loop(415) = {26, -65, 45, 6};
Plane Surface(416) = {415};
x_pos = 416;

Line Loop(417) = {222, -182, 162, -202};
Plane Surface(418) = {417};
y_pos = 418;

// The whole cylinder is now
shell[] += top_cap[];
shell[] += bot_cap[];
cylinder[] = {y_neg, x_neg, y_pos, x_pos};
cylinder[] += shell[];

sl = newsl;
Surface Loop(sl) = {cylinder[]};
v = newv;
Volume(v) = {sl};
pv = 1;
Physical Volume(pv) = {v};
Physical Surface(pv) = shell[]; // Don't forget to add caps
Physical Surface(pv) += {x_neg}; // bottom == x_neg
Physical Surface(pv) += {y_neg};                // left == y_neg

n_cylinders = x_cylinders*y_cylinders;
Physical Surface(pv+n_cylinders) = {y_pos};
// And x shared
Physical Surface(pv+2*n_cylinders) = {x_pos};

///////////////////////////////////////////////////////////////////////////////
// Grid generation
///////////////////////////////////////////////////////////////////////////////
all_cylinders[] = {v};
parent = v;
For i In {1:x_cylinders}
  // Mirror in y
  last_volume = parent;
  For j In {2:y_cylinders}
    new_volume[] = Translate {0, 2*(R+h/2), 0} { Duplicata { Volume{last_volume}; } };
    last_volume = new_volume[0];
    all_cylinders[] += {last_volume};

    pv += 1;
    Physical Volume(pv) = {last_volume};
    
    // The surfaces that I own
    b() = Boundary{ Volume{last_volume}; } ;
    not_shared[] = {};
    For k In {4:99}
        not_shared[] += {b[k]};
    EndFor
    Physical Surface(pv) = not_shared[];
    // If this is last I claim right y_pos
    If(j == y_cylinders)
      Physical Surface(pv) += {b[2]};
    EndIf
    // If this is x-first we take bottom = x_neg
    If(i == 1)
      Physical Surface(pv) += {b[1]};
    EndIf
    // If this is x_last we take top = x_pos
    If(i == x_cylinders)
      Physical Surface(pv) += {b[3]};
    EndIf
    // Not the shared surfaces	
    // Claim y-shared, that is right
    If(j < y_cylinders)
        Physical Surface(pv+n_cylinders) = {b[2]};
    EndIf
    // Claim x-shared x_pos
    If(i < x_cylinders)
        Physical Surface(pv+2*n_cylinders) = {b[3]};
    EndIf

  EndFor 

  If(i < x_cylinders)
    new_volume[] = Translate {2*(R+h/2), 0, 0} { Duplicata { Volume{parent}; } };
    parent = new_volume[0];
    all_cylinders[] += {parent};

    pv += 1;
    Physical Volume(pv) = {parent};

    // The surfaces that I own
    b() = Boundary{ Volume{parent}; } ;
    not_shared[] = {};
    For k In {4:99}
        not_shared[] += {b[k]};
    EndFor
    Physical Surface(pv) = not_shared[];

    // Always claim left y_neg
    Physical Surface(pv) += {b[0]};
    // The right is shared
    Physical Surface(pv+n_cylinders) = {b[2]};

    If((i+1) == x_cylinders)
      Physical Surface(pv) += {b[3]};
    EndIf
    // Claim x-shared x_pos
    If((i + 1)< x_cylinders)
        Physical Surface(pv+2*n_cylinders) = {b[3]};
    EndIf

  EndIf

EndFor
// Get the bounding surfaces for cylinder
bbox_cylinders = CombinedBoundary{ Volume{all_cylinders[]}; } ;

bbox_surfaces[] = {};
///////////////////////////////////////////////////////////////////////////////
// Bounding box
///////////////////////////////////////////////////////////////////////////////
dx = R+h/2+dx;
dz = H/2+dz;
// Bottom corner
p = newp;
Point(p) =     {-dx,                        -dx, -dz, size_b};
Point(p+1) =   {-dx,                        (2*R+h)*(y_cylinders-1)+dx, -dz, size_b};
Point(p+2) =   {(2*R+h)*(x_cylinders-1)+dx, (2*R+h)*(y_cylinders-1)+dx, -dz, size_b};
Point(p+3) =   {(2*R+h)*(x_cylinders-1)+dx, -dx, -dz, size_b};
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
Point(P) =     {-dx,                        -dx, dz, size_b};
Point(P+1) =   {-dx,                        (2*R+h)*(y_cylinders-1)+dx, dz, size_b};
Point(P+2) =   {(2*R+h)*(x_cylinders-1)+dx, (2*R+h)*(y_cylinders-1)+dx, dz, size_b};
Point(P+3) =   {(2*R+h)*(x_cylinders-1)+dx, -dx, dz, size_b};
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
