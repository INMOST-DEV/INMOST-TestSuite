// Gmsh project
// Domain to emulate ventricular appendage with internal muscles
SetFactory("OpenCASCADE");
x=-0.2; //torus center
y=-0;   //torus center
r1 = 0.9; // torus radius
r2 = 0.16; // torus section radius
alpha = Pi*0.35; //torus angle
rc = 0.01; //radius of cylinders
tc = 0.6; //offset cylinders by section radius
rc1 = r1 - tc*r2; // radius to place even cylinders
rc2 = r1 + tc*r2; // radius to place odd cylinders
h = 0.16; //cylinders height
N = 6; //number of cylinders
//+
Torus(1) = {x, y, 0, r1, r2, alpha};
list() = {}; //empty list
ic = 2;
For k In {1:N}
	beta = alpha*(k-0.5)/N; //angle for cylinder positions
	xc = x + Cos(beta)*rc1;
	yc = y + Sin(beta)*rc1;
	Cylinder(ic) = {xc,yc,-h,0,0,2*h,rc,2*Pi};
	list() += ic;
	Printf("Cylinder %g (k %g) at %g %g angle %g",ic,k,xc,yc,beta);
	ic++;
	
	xc = x + Cos(beta)*rc2;
	yc = y + Sin(beta)*rc2;
	Cylinder(ic) = {xc,yc,-h,0,0,2*h,rc,2*Pi};
	list() += ic;
	Printf("Cylinder %g (k %g) at %g %g angle %g",ic,k,xc,yc,beta);
	ic++;
EndFor

//+
BooleanDifference{ Volume{1}; Delete; }{ Volume{list()}; Delete; }
//+
Physical Surface(1) = {3}; //inlet surface
Physical Surface(2) = {2}; //outlet surface

Field[1] = Box;
Field[1].VIn = 0.015;
Field[1].VOut = 0.1;
Field[1].XMax = 1;
Field[1].YMax = 1;
Field[1].ZMax = 1;
Field[1].ZMin = -1;
Background Field = 1;

