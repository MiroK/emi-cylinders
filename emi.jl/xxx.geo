SIZE=0.1;

Point(1) = {0, 0, 0., SIZE};
Point(2) = {1, 0, 0., SIZE};
Point(3) = {1, 1, 0., SIZE};
Point(4) = {0, 1, 0., SIZE};
Line(5) = {2, 1};
Line(6) = {3, 2};
Line(7) = {4, 3};
Line(8) = {1, 4};
Line Loop(9) = {5, 6, 7, 8};
Physical Surface(2) = {9};
Physical Line(2) = {5, 6, 7, 8};
