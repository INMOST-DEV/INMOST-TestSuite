#include "test7.h"

using namespace INMOST;

//shortcut for types
typedef Storage::real real;
typedef Storage::enumerator enumerator;
typedef Storage::var_array var_array;
typedef Storage::real_array real_array;
typedef Storage::reference_array ref_array;

double Test7::Gravity() const {return 1;}
double Test7::FluidDensity() const {return 1;}
double Test7::FluidViscosity() const {return 1;}
double Test7::SolidDensity() const {return 2;}
double Test7::Porosity(double _x, double _y, double _z) const {return 0.25;}
double Test7::InverseBiotModulus(double _x, double _y, double _z) const {return 1;}



Test7::Test7()
{
	double ratioE = 50;
	ratio = 100;
	real E1, E2, E3, nu12, nu21, nu13, nu31, nu23, nu32, G23, G13, G12;
	rMatrix C(6,6);
	C.Zero();
	
	E1 = E2 = 1;
	E3 = 1.0*ratioE;
	nu12 = nu21 = 1.0/ratioE;
	nu31 = 2.0;///ratio;
	nu13 = nu31*E1/E3;
	nu32 = nu31;
	nu23 = nu13;
	G23 = G13 = 0.5;
	G12 = 0.5*ratioE;
	C(0,0) =   1.0/E1;
	C(0,1) = -nu21/E2;
	C(0,2) = -nu31/E3;
	C(1,0) = -nu12/E1;
	C(1,1) =   1.0/E2;
	C(1,2) = -nu32/E3;
	C(2,0) = -nu13/E1;
	C(2,1) = -nu23/E2;
	C(2,2) =   1.0/E3;
	C(3,3) = 1.0/(2.0*G23);
	C(4,4) = 1.0/(2.0*G13);
	C(5,5) = 1.0/(2.0*G12);
	std::cout << "C:" <<std::endl;
	C.Print();
	rCheckEigen(C);
	std::cout << "E:" <<std::endl;
	E0 = C.Invert();
	E0.Print();
	rCheckEigen(E0);
	K0.Resize(3,3);
	B0.Resize(3,3);
	K0.Zero(); B0.Zero();
	B0(0,0) = B0(1,1) = K0(0,0) = K0(1,1) = 1;
	K0(2,2) = ratio;
	B0(2,2) = ratio;
}

vMatrix Test7::ElasticTensor(double _x, double _y, double _z) const
{
	return E0;
}

vMatrix Test7::BiotTensor(double _x, double _y, double _z) const
{
	return B0;
}

vMatrix Test7::PermTensor(double _x, double _y, double _z) const
{
	return K0;
}

hMatrix Test7::Displacement(double _x, double _y, double _z, double _t) const
{
	unknown x(_x,0), y(_y,1), z(_z,2), t(_t,3);
	hMatrix sol(3,1);
	const real pi = 3.1415926536;
	sol(0,0) = sin(2*pi*x)*cos(pi*t/8.0); //u
	sol(1,0) = sin(2*pi*y)*cos(pi*t/8.0); //v
	sol(2,0) = sin(2*pi*z)*cos(pi*t/8.0); //w
	
	return sol;
}

hessian_variable Test7::Pressure(double _x, double _y, double _z, double _t) const
{
	const real pi = 3.1415926536;
	unknown x(_x,0), y(_y,1), z(_z,2), t(_t,3);
	//return sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z)*cos(pi*t/6.0);
	return exp(-2*pi*sqrt(1.0/ratio)*x)*exp(-2*pi*sqrt(1.0/ratio)*y)*sin(2*pi*z)*cos(pi*t/6.0);
}

rMatrix Test7::BCFlow_coef(double _x, double _y, double _z, double _t) const
{
	rMatrix ret(2,1);
	if( 1 - _x - _y - _z > 0 )
	{
		ret(0,0) = 0; //dirichlet parallel
		ret(1,0) = 1; //neumann parallel
	}
	else
	{
		ret(0,0) = 1; //dirichlet
		ret(1,0) = 0;
	}
	return ret;
}

rMatrix Test7::BCMech_coef(double _x, double _y, double _z, double _t) const
{
	rMatrix ret(4,1);
	if( 1 - _x - _y - _z > 0 )
	{
		ret(0,0) = 0; //dirichlet perp
		ret(1,0) = 1; //neumann perp
		ret(2,0) = 0; //dirichlet parallel
		ret(3,0) = 1; //neumann parallel
	}
	else
	{
		ret(0,0) = 1; //dirichlet perp
		ret(1,0) = 0; //neumann perp
		ret(2,0) = 1; //dirichlet parallel
		ret(3,0) = 0; //neumann parallel
	}
	return ret;
}










