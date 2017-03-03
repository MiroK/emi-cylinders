Point(1) = {0, 0, 0, 1};
Point(2) = {0, 1, 0, 1};
Point(3) = {1, 1, 0, 1};
Point(4) = {1, 0, 0, 1};
Point(5) = {0, 0, 1, 1};
Point(6) = {0, 1, 1, 1};
Point(7) = {1, 1, 1, 1};
Point(8) = {1, 0, 1, 1};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {1, 5};
Line(6) = {5, 8};
Line(7) = {8, 4};
Line(8) = {3, 7};
Line(9) = {7, 8};
Line(10) = {5, 6};
Line(11) = {6, 7};
Line(12) = {6, 2};

Line Loop(13) = {6, 7, 4, 5};
Plane Surface(14) = {13};
Line Loop(15) = {7, -3, 8, 9};
Plane Surface(16) = {15};
Line Loop(17) = {8, -11, 12, 2};
Plane Surface(18) = {17};
Line Loop(19) = {12, -1, 5, 10};
Plane Surface(20) = {19};
Line Loop(21) = {10, 11, 9, -6};
Plane Surface(22) = {21};
Line Loop(23) = {2, 3, 4, 1};
Plane Surface(24) = {23};
Surface Loop(25) = {22, 20, 18, 16, 14, 24};
Volume(26) = {25};

Out[] = Translate {0, 0, 2} {
  Duplicata { Surface{22, 20, 18, 16, 14, 24};}
};

Printf("%g", Out[0]);
Printf("%g", Out[1]);
