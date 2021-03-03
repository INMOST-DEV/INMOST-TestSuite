#include "test2.h"

using namespace INMOST;

//shortcut for types
typedef Storage::real real;
typedef Storage::enumerator enumerator;
typedef Storage::var_array var_array;
typedef Storage::real_array real_array;
typedef Storage::reference_array ref_array;

double Test2::Gravity() const {return -9.8;}
double Test2::FluidDensity() const {return 978;}
double Test2::FluidViscosity() const {return 0.2;}
double Test2::SolidDensity() const {return 2500;}
double Test2::Porosity(double _x, double _y, double _z) const {return 0.25;}
double Test2::InverseBiotModulus(double _x, double _y, double _z) const {return 0.1;}

Test2::Test2()  {}

vMatrix Test2::ElasticTensor(double _x, double _y, double _z) const
{
	unknown x(_x,0), y(_y,1), z(_z,2);
	real a = 0.01, b = 0.1, c = 1, d = 0.002;
	vMatrix C(6,6);
	double Cv[21] =
	{
		1.323 , 0.0726 , 0.263 , 0.108 , -0.08 , -0.239,
		1.276 , -0.318 , 0.383 , 0.108 , 0.501 ,
		0.943 , -0.183 , 0.146 , 0.182 ,
		1.517 , -0.0127 , -0.304 ,
		1.209 , -0.326 ,
		1.373
	};
	C = rMatrix::FromTensor(Cv,21,6);
	//vCheckEigen(C);
	return C;
}

vMatrix Test2::BiotTensor(double _x, double _y, double _z) const
{
	unknown x(_x,0), y(_y,1), z(_z,2);
	vMatrix B(3,3);
	B(0,0) = B(1,1) = B(2,2) = 1.5;
	B(1,0) = B(0,1) = 0.1;
	B(2,0) = B(0,2) = 0.5;
	B(1,2) = B(2,1) = 0.15;
	//vCheckEigen(B);
	return B;
}

vMatrix Test2::PermTensor(double _x, double _y, double _z) const
{
	unknown x(_x,0), y(_y,1), z(_z,2);
	vMatrix K(3,3);
	K(0,0) = K(1,1) = K(2,2) = 1.5;
	K(1,0) = K(0,1) = 0.5;
	K(2,0) = K(0,2) = 0.35;
	K(1,2) = K(2,1) = 0.45;
	//vCheckEigen(K);
	return K;
}

hMatrix Test2::Displacement(double _x, double _y, double _z, double _t) const
{
	unknown x(_x,0), y(_y,1), z(_z,2), t(_t,3);
	hMatrix sol(3,1);
	
	sol(0,0) = ((x-0.5)*(x-0.5)-y-z)*(1+t*t); //u
	sol(1,0) = ((y-0.5)*(y-0.5)-x-z)*(1+t*t); //v
	sol(2,0) = ((z-0.5)*(z-0.5)-x-y)*(1+t*t); //w
	
	return sol;
}

hessian_variable Test2::Pressure(double _x, double _y, double _z, double _t) const
{
	const real pi = 3.1415926536;
	unknown x(_x,0), y(_y,1), z(_z,2), t(_t,3);
	return sin((1-x)*(1-y)*(1-z))/(2.0*sin(1.0))+0.5*(1-x)*(1-x)*(1-x)*(1-y)*(1-y)*(1-z)*(1+t*t);
}

rMatrix Test2::BCFlow_coef(double _x, double _y, double _z, double _t) const
{
	rMatrix ret(2,1);
	ret(0,0) = 1; //dirichlet
	ret(1,0) = 0;
	return ret;
}

rMatrix Test2::BCMech_coef(double _x, double _y, double _z, double _t) const
{
	rMatrix ret(4,1);
	ret(0,0) = 1; //dirichlet perp
	ret(1,0) = 0;
	ret(2,0) = 1; //dirichlet parallel
	ret(3,0) = 0;
	return ret;
}










