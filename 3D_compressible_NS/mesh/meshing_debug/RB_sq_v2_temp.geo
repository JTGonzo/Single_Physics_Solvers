cl__1=1;
dim1 = 0.015;
dim2 = 0.0005;
nElem = 7;
Point(1) = {-dim1, -dim1, -dim1, cl__1};
Point(2) = {dim1, -dim1, -dim1, cl__1};
Point(3) = {dim1, dim1, -dim1, cl__1};
Point(4) = {-dim1, dim1, -dim1, cl__1};
Point(5) = {0.0, 0.0, -dim1, cl__1};
Point(6) = {-dim2, -dim2, -dim1, cl__1};
Point(7) = {dim2, -dim2, -dim1, cl__1};
Point(8) = {dim2, dim2, -dim1, cl__1};
Point(9) = {-dim2, dim2, -dim1, cl__1};
Point(10) = {dim2, dim1, -dim1, cl__1};
Point(11) = {-dim2, dim1, -dim1, cl__1};
Point(12) = {-dim1, dim2, -dim1, cl__1};
Point(13) = {-dim1, -dim2, -dim1, cl__1};
Point(14) = {-dim2, -dim1, -dim1, cl__1};
Point(15) = {dim2, -dim1, -dim1, cl__1};
Point(16) = {dim1, -dim2, -dim1, cl__1};
Point(17) = {dim1, dim2, -dim1, cl__1};

Line(23) = {3, 10};
Line(24) = {10, 11};
Line(25) = {11, 4};
Line(26) = {4, 12};
Line(27) = {12, 13};
Line(28) = {13, 1};
Line(29) = {1, 14};
Line(30) = {14, 15};
Line(31) = {15, 2};
Line(32) = {2, 16};
Line(33) = {16, 17};
Line(34) = {17, 3};
Line(35) = {10, 8};
Line(36) = {8, 17};
Line(37) = {16, 7};
Line(38) = {7, 15};
Line(39) = {14, 6};
Line(40) = {6, 13};
Line(41) = {12, 9};
Line(42) = {9, 11};
Line(43) = {8, 9};
Line(44) = {9, 6};
Line(45) = {6, 7};
Line(46) = {7, 8};

Line Loop(47) = {23, 35, 36, 34};
Plane Surface(48) = {47};
Line Loop(49) = {25, 26, 41, 42};
Plane Surface(50) = {49};
Line Loop(51) = {40, 28, 29, 39};
Plane Surface(52) = {51};
Line Loop(53) = {31, 32, 37, 38};
Plane Surface(54) = {53};
Line Loop(55) = {24, -42, -43, -35};
Plane Surface(56) = {55};
Line Loop(57) = {27, -40, -44, -41};
Plane Surface(58) = {57};
Line Loop(59) = {30, -38, -45, -39};
Plane Surface(60) = {59};
Line Loop(61) = {33, -36, -46, -37};
Plane Surface(62) = {61};
Line Loop(63) = {43, 44, 45, 46};
Plane Surface(64) = {63};

Transfinite Line {-23, 25, -26, 28, -29, 31, -32, 36, -35, 42, -41, 40, -39, 38, -37, 34} = nElem Using Progression 1;
Transfinite Line {24, 27, 30, 33, 43, 44, 45, 46} = 7 Using Progression 1;
Transfinite Surface {48, 50, 52, 54, 56, 58, 60, 62, 64};
Recombine Surface {48, 50, 52, 54, 56, 58, 60, 62, 64};

Extrude {0, 0, dim1-dim2} {
  Surface{48, 50, 52, 54, 56, 58, 60, 62, 64};Layers{6};Recombine;
}

Extrude {0, 0, 2*dim2} {
  Surface{130, 86, 152, 196, 218, 240, 174, 262, 108};Layers{6};Recombine;
}

Extrude {0, 0, dim1-dim2} {
  Surface{306, 460, 438, 416, 328, 372, 350, 394, 284};Layers{6};Recombine;
}

Physical Surface("LeftBoundary") = {183,275,99,121,495,601,649,337,451};
Physical Surface("RightBoundary") = {227,305,319,85,143,561,623,381,481};
Physical Surface("TopBoundary") = {161,293,403,73,95,535,447,469,491};
Physical Surface("BottomBoundary") = {205,279,315,125,139,557,579,653,359};
Physical Surface("FrontBoundary") = {504, 526, 548, 570, 592, 614, 636, 658, 482};
Physical Surface("BackBoundary") = {48, 50, 52, 54, 56, 58, 60, 62, 64};
Physical Volume("fluid") = {1:27};
