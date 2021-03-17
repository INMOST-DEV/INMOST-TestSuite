// Gmsh
//SetFactory("OpenCASCADE");
L=0.137; //crevice size (mm)
C=1.8+L;
h=0.1;
Vin = 0.002;
Vout = 0.02;
Vnodes=0.002;
Point(1) = {0,0,0};
Point(2) = {1.5,0,0,Vnodes};
Point(3) = {1.5,-0.125,0,Vnodes};
Point(4) = {1.5+L,-0.125,0,Vnodes};
Point(5) = {1.5+L,0,0,Vnodes};
Point(6) = {1.8+L,0,0};
Point(7) = {1.8+L,1.5,0};
Point(8) = {0,1.5,0};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,1};



Curve Loop(1) = {1,2,3,4,5,6,7,8};
Plane Surface(1) = {1};

Recombine Surface {1};



//Transfinite Line {2,5} = 20 Using Progression 1;


Extrude {0,0,h} {
  Surface{1}; Layers{1}; Recombine;
}



//refine in crevice
Field[1] = Box;
Field[1].VIn = Vin;
Field[1].VOut = Vout;
Field[1].XMin = 1.5-0.065;
Field[1].XMax = 1.5+L+0.065;
Field[1].YMin = -0.125;
Field[1].YMax = 0.075;
Field[1].ZMin = -0.1;
Field[1].ZMax = 0.2;

//+
Background Field = 1;

Mesh.CharacteristicLengthFromPoints = 1;
Mesh.CharacteristicLengthFromCurvature = 1;
Mesh.CharacteristicLengthExtendFromBoundary = 1;

//Mesh.Algorithm = 8;
//Mesh 3;
