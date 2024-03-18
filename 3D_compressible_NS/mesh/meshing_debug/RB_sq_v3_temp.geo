cl__1=1;
dim1 = 0.001;
dim2 = 0.0005;
nElem = 42;

Point(1) = {-dim1, -dim1, -dim1, cl__1};
Point(2) = {dim1, -dim1, -dim1, cl__1};
Point(3) = {dim1, dim1, -dim1, cl__1};
Point(4) = {-dim1, dim1, -dim1, cl__1};
Point(5) = {-dim1, -dim1, dim1, cl__1};
Point(6) = {dim1, -dim1, dim1, cl__1};
Point(7) = {dim1, dim1, dim1, cl__1};
Point(8) = {-dim1, dim1, dim1, cl__1};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

// back surface
Line Loop(13) = {1, 2, 3, 4};
Plane Surface(14) = {13};

// front surface
Line Loop(15) = {5, 6, 7, 8};
Plane Surface(16) = {15};

// right surface
Line Loop(17) = {-10, 2, 11, -6};
Plane Surface(18) = {17};

// left surface
Line Loop(19) = {-9, -4, 12, 8};
Plane Surface(20) = {19};

// top surface
Line Loop(21) = {-7, -11, 3, 12};
Plane Surface(22) = {21};

// bottom surface
Line Loop(23) = {-1, 9, 5, -10};
Plane Surface(24) = {23};

Surface Loop(25) = {18, 24, 14, 22, 16, 20};
Volume(26) = {25};

Transfinite Line {-1, -2, 3, 4, -5, -6, 7, 8, -9, -10, -11, -12} = nElem Using Progression 1;
Transfinite Surface {14, 16, 18, 20, 22, 24};
Recombine Surface {14, 16, 18, 20, 22, 24};

Physical Surface("LeftBoundary") = {20};
Physical Surface("RightBoundary") = {18};
Physical Surface("TopBoundary") = {22};
Physical Surface("BottomBoundary") = {24};
Physical Surface("FrontBoundary") = {16};
Physical Surface("BackBoundary") = {14};
Physical Volume("fluid") = {26};


