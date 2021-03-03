#include "test3.h"

using namespace INMOST;

//shortcut for types
typedef Storage::real real;
typedef Storage::enumerator enumerator;
typedef Storage::var_array var_array;
typedef Storage::real_array real_array;
typedef Storage::reference_array ref_array;

double Test3::Gravity() const {return 0;}
double Test3::FluidDensity() const {return 1;}
double Test3::FluidViscosity() const {return 1;}
double Test3::SolidDensity() const {return 2;}
double Test3::Porosity(double _x, double _y, double _z) const {return 0.2;}
double Test3::InverseBiotModulus(double _x, double _y, double _z) const {return 1;}



Test3::Test3()// : C0(GenMatrix1(6)), K0(GenMatrix1(3)), B0(GenMatrix1(3))
{
	
	real vE0[36] =
	{
		93,           46,           22,           13,           72,           35,
		46,           95,           41,           62,           56,           24,
		22,           41,           89,           25,           33,           21,
		13,           62,           25,           87,           13,           25,
		72,           56,           33,           13,           99,           57,
		35,           24,           21,           25,           57,           78
	};
	real vE1[36] =
	{
		81,           33,            5,           14,           37,           32,
		33,          100,           26,           58,           58,           73,
		5,           26,           90,           23,           60,           28,
		14,           58,           23,           77,           29,            3,
		37,           58,           60,           29,           84,           60,
		32,           73,           28,            3,           60,           99

	};
	real vK0[9] =
	{
		25,            2,           39,
		2,           42,            7,
		39,            7,          100
	};
	real vK1[9] =
	{
		1,            2,            2,
		2,           68,           42,
		2,           42,           49
	};
	real vB0[9] =
	{
		1,            6,            5,
		6,           67,           27,
		5,           27,           76
	};
	real vB1[9] =
	{
		1,            1,            1,
		1,           68,           13,
		1,           13,           73
	};
	E0 = rMatrix(vE0,6,6);
	K0 = rMatrix(vK0,3,3);
	B0 = rMatrix(vB0,3,3);
	E1 = rMatrix(vE1,6,6);
	K1 = rMatrix(vK1,3,3);
	B1 = rMatrix(vB1,3,3);
	/*
	int seed = 0;
	seed = GenTensor(E0,6,seed);
	seed = GenTensor(E1,6,seed+1);
	seed = 0;
	seed = GenTensor(K0,3,seed);
	seed = GenTensor(B0,3,seed+1);
	seed = GenTensor(K1,3,seed+1);
	seed = GenTensor(B1,3,seed+1);
	 */
	std::cout << "E0: " << std::endl;
	E0.Print();
	rCheckEigen(E0,true);
	std::cout << "K0: " << std::endl;
	K0.Print();
	rCheckEigen(K0,true);
	std::cout << "B0: " << std::endl;
	B0.Print();
	rCheckEigen(B0,true);
	std::cout << "E1: " << std::endl;
	E1.Print();
	rCheckEigen(E1,true);
	std::cout << "K1: " << std::endl;
	K1.Print();
	rCheckEigen(K1,true);
	std::cout << "B1: " << std::endl;
	B1.Print();
	rCheckEigen(B1,true);
	real rho_f = FluidDensity();
	real g = Gravity();
	rMatrix T0(3,3), T1(3,3); //3 by 3 projections in materials A and B
	rMatrix M0(3,9), M1(3,9); //transversal directions in materials A and B
	rMatrix n(3,1), x(3,1), z(3,1);
	n(0,0) = 1;
	n(1,0) = 0;
	n(2,0) = 0;
	x(0,0) = 0.5;
	x(1,0) = 0.5;
	x(2,0) = 0.5;
	z(0,0) = 0;
	z(1,0) = 0;
	z(2,0) = 1;
	
	SplitConormals(real_array(E0.data(),36),n,T0,M0,false);
	SplitConormals(real_array(E1.data(),36),n,T1,M1,false);
	real mu = FluidViscosity();
	
	real l0 = n.DotProduct(K0*n);
	real l1 = n.DotProduct(K1*n);
	real d = n.DotProduct(x);
	
	real ps;
	//rMatrix gp1(3,1), gp2(3,1);
	gp1.Resize(3,1);
	gp2.Resize(3,1);
	p1 = 1.0;
	gp1(0,0) = 2;
	gp1(1,0) = 0;
	gp1(2,0) = 0;
	ps = p1 + gp1.DotProduct(x);
	gp2 = (rMatrix::Unit(3) + 1.0/l1*(K0-K1)*n*n.Transpose()).Transpose()*gp1;
	gp2-= rho_f*g/l1*((K0-K1)*n*n.Transpose()).Transpose()*z;
	p2 = p1 - d/l1*((K0-K1)*n).DotProduct(gp1);
	p2+= rho_f*g*d/l1*((K0-K1)*n).DotProduct(z);
	
	
	//rMatrix gu1(9,1), gu2(9,1);
	//rMatrix u1(3,1), u2(3,1);
	u1.Resize(3,1);
	u2.Resize(3,1);
	gu1.Resize(9,1);
	gu2.Resize(9,1);
	
	real vu1[3] = {1,14,76};
	real vgu1[9] = {46, 54, 22, 5, 68, 68, 94, 39, 52};
	
	u1 = rMatrix(vu1,3,1);
	gu1 = rMatrix(vgu1,9,1);
	
	//for(int k = 0; k < 3; ++k)
	//	u1(k,0) = ceil((rand()/(1.0*RAND_MAX))*100);
	
	//for(int k = 0; k < 9; ++k)
	//	gu1(k,0) = ceil((rand()/(1.0*RAND_MAX))*100);
	
	gu2 = (rMatrix::Unit(9) + rMatrix::Unit(3).Kronecker(n)*T1.Invert()*(M0-M1))*gu1 + rMatrix::Unit(3).Kronecker(n)*T1.Invert()*(B1-B0)*n*ps;
	u2 = u1 - d*T1.Invert()*(M0-M1)*gu1 -  d*T1.Invert()*(B1-B0)*n*ps;
	
	std::cout << "p1 " << p1 << " gp1: ";
	gp1.Transpose().Print();
	std::cout << "p2 " << p2 << " gp2: ";
	gp2.Transpose().Print();
	std::cout << "u1 " << u1(0,0) << " " << u1(1,0) << " " << u1(2,0) << " gu1: ";
	gu1.Transpose().Print();
	std::cout << "u2 " << u2(0,0) << " " << u2(1,0) << " " << u2(2,0) << " gu2: ";
	gu2.Transpose().Print();
	std::cout << "check argument continuity at " << x(0,0) << " " << x(1,0) << " " << x(2,0) <<  std::endl;
	std::cout << "p: " << p1+x.DotProduct(gp1) << " " << p2+x.DotProduct(gp2) << std::endl;
	std::cout << "u1: ";
	(u1+rMatrix::Unit(3).Kronecker(x.Transpose())*gu1).Transpose().Print();
	std::cout << "u2: ";
	(u2+rMatrix::Unit(3).Kronecker(x.Transpose())*gu2).Transpose().Print();
	std::cout << "check flux continuity at " << x(0,0) << " " << x(1,0) << " " << x(2,0) <<  std::endl;
	std::cout << "p: " << (K0*n).DotProduct(gp1) << " " << (K1*n).DotProduct(gp2) << std::endl;
	std::cout << "u1: ";
	(M0*gu1 - B0*n*(p1+x.DotProduct(gp1))).Transpose().Print();
	std::cout << "u2: ";
	(M1*gu2 - B1*n*(p2+x.DotProduct(gp2))).Transpose().Print();
	x(1,0) = 0.1;
	x(2,0) = 0.3;
	std::cout << "check argument continuity at " << x(0,0) << " " << x(1,0) << " " << x(2,0) <<  std::endl;
	std::cout << "p: " << p1+x.DotProduct(gp1) << " " << p2+x.DotProduct(gp2) << std::endl;
	std::cout << "u1: ";
	(u1+rMatrix::Unit(3).Kronecker(x.Transpose())*gu1).Transpose().Print();
	std::cout << "u2: ";
	(u2+rMatrix::Unit(3).Kronecker(x.Transpose())*gu2).Transpose().Print();
	std::cout << "check flux continuity at " << x(0,0) << " " << x(1,0) << " " << x(2,0) <<  std::endl;
	std::cout << "p: " << (K0*n).DotProduct(gp1) << " " << (K1*n).DotProduct(gp2) << std::endl;
	std::cout << "u1: ";
	(M0*gu1 - B0*n*(p1+x.DotProduct(gp1))).Transpose().Print();
	std::cout << "u2: ";
	(M1*gu2 - B1*n*(p2+x.DotProduct(gp2))).Transpose().Print();
}

vMatrix Test3::ElasticTensor(double _x, double _y, double _z) const {return (_x < 0.5 ? E0 : E1);}
vMatrix Test3::BiotTensor(double _x, double _y, double _z) const { return (_x < 0.5 ? B0 : B1); }
vMatrix Test3::PermTensor(double _x, double _y, double _z) const
{
	return (_x < 0.5 ? K0 : K1);
}

hMatrix Test3::Displacement(double _x, double _y, double _z, double _t) const
{
	unknown x(_x,0), y(_y,1), z(_z,2), t(_t,3);
	vMatrix xyz(3,1);
	xyz(0,0)= x;
	xyz(1,0)= y;
	xyz(2,0)= z;
	hMatrix sol(3,1);
	if( _x < 0.5 ) return u1 + rMatrix::Unit(3).Kronecker(xyz.Transpose())*gu1;
	else return u2 + rMatrix::Unit(3).Kronecker(xyz.Transpose())*gu2;
		
}

hessian_variable Test3::Pressure(double _x, double _y, double _z, double _t) const
{
	const real pi = 3.1415926536;
	unknown x(_x,0), y(_y,1), z(_z,2), t(_t,3);
	vMatrix xyz(3,1);
	xyz(0,0)= x;
	xyz(1,0)= y;
	xyz(2,0)= z;
	//return 16*x+17*y+18*z+19*t+20;
	//return 9*x+10*y+11*z+1;
	//return (x*x*x*y*y*z + x*sin(2*pi*x*z)*sin(2*pi*x*y)*sin(2*pi*z))*exp(-t);
	if( _x < 0.5 ) return p1 + xyz.DotProduct(gp1);
	else return p2 + xyz.DotProduct(gp2);
}


rMatrix Test3::BCFlow_coef(double _x, double _y, double _z, double _t) const
{
	rMatrix ret(2,1);
	ret(0,0) = 1; //dirichlet
	ret(1,0) = 0;
	return ret;
}

rMatrix Test3::BCMech_coef(double _x, double _y, double _z, double _t) const
{
	rMatrix ret(4,1);
	ret(0,0) = 1; //dirichlet perp
	ret(1,0) = 0;
	ret(2,0) = 1; //dirichlet parallel
	ret(3,0) = 0;
	return ret;
}










