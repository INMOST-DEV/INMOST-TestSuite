#include "test4.h"

using namespace INMOST;

//shortcut for types
typedef Storage::real real;
typedef Storage::enumerator enumerator;
typedef Storage::var_array var_array;
typedef Storage::real_array real_array;
typedef Storage::reference_array ref_array;

double Test4::Gravity() const {return -9.8;}
double Test4::FluidDensity() const {return 1;}
double Test4::FluidViscosity() const {return 1;}
double Test4::SolidDensity() const {return 2;}
double Test4::Porosity(double _x, double _y, double _z) const {return 0.2 + 0.2 * sin(_x) * _y;}
double Test4::InverseBiotModulus(double _x, double _y, double _z) const {return 1 + sin(_z) * _x;}

Test4::Test4()  {}

vMatrix Test4::ElasticTensor(double _x, double _y, double _z) const
{
	unknown x(_x,0), y(_y,1), z(_z,2);
	real a = 0.01, b = 0.1, c = 1, d = 0.002;
	vMatrix C(6,6);
	C(0,0) = a*x*x + b*y*y + c*z*z+d;
	C(1,1) = c*x*x + a*y*y + b*z*z+d;
	C(2,2) = b*x*x + c*y*y + a*z*z+d;
	C(3,3) = c*x*x + b*y*y + a*z*z+d;
	C(4,4) = a*x*x + c*y*y + b*z*z+d;
	C(5,5) = b*x*x + a*y*y + c*z*z+d;
	C(0,1) = C(1,0) = -a*x*y;
	C(0,2) = C(2,0) = -b*x*z;
	C(0,3) = C(3,0) = -a*y*z;
	C(0,4) = C(4,0) = -a*(x*y+y*z);
	C(0,5) = C(5,0) = -a*(x*z+y*z);
	C(1,2) = C(2,1) = -b*(x*y+x*z);
	C(1,3) = C(3,1) = -a*(x*z);
	C(1,4) = C(4,1) = -a*(y*z);
	C(1,5) = C(5,1) = -b*(x*y);
	C(2,3) = C(3,2) = -a*(x*y+x*z);
	C(2,4) = C(4,2) = -a*(y*z+x*y);
	C(2,5) = C(5,2) = -b*(x*z+y*z);
	C(3,4) = C(4,3) = -a*(x*y*z);
	C(3,5) = C(5,3) = -a*(x*z-y*z);
	C(4,5) = C(5,4) = -a*(x*y*z);
	//vCheckEigen(C);
	return C;
}

vMatrix Test4::BiotTensor(double _x, double _y, double _z) const
{
	unknown x(_x,0), y(_y,1), z(_z,2);
	vMatrix B(3,3);
	// B as a function of C?
	real a = 0.0;
	variable div = 1; //((x-1)*(x-1)+(y-1)*(y-1)+(z-1)*(z-1)+1.0e-5);
	B(0,0) = (a*x*x + z*z + y*y + 1)/div;
	B(0,1) = B(1,0) = (a-1)*x*y/div;
	B(0,2) = B(2,0) = (a-1)*x*z/div;
	B(1,1) = (a*y*y + x*x + z*z + 1)/div;
	B(1,2) = B(2,1) = (a-1)*y*z/div;
	B(2,2) = (a*z*z + x*x + y*y + 1)/div;
	//vCheckEigen(B);
	return B;
}

vMatrix Test4::PermTensor(double _x, double _y, double _z) const
{
	unknown x(_x,0), y(_y,1), z(_z,2);
	vMatrix K(3,3);
	real a = 0.5;
	variable div = 1;//(x*x+y*y+z*z+1.0e-5);
	K(0,0) = (a*x*x + z*z + y*y + 1)/div;
	K(0,1) = K(1,0) = (a-1)*x*y/div;
	K(0,2) = K(2,0) = (a-1)*x*z/div;
	K(1,1) = (a*y*y + x*x + z*z + 1)/div;
	K(1,2) = K(2,1) = (a-1)*y*z/div;
	K(2,2) = (a*z*z + x*x + y*y + 1)/div;
	//vCheckEigen(K);
	return K;
}

hMatrix Test4::Displacement(double _x, double _y, double _z, double _t) const
{
	unknown x(_x,0), y(_y,1), z(_z,2), t(_t,3);
	hMatrix sol(3,1);
	const real a = 2*3.1415926536;
	sol(0,0) = 0.05*(sin(a*(x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) + (z-0.5)*(z-0.5))+1)*exp(-t); //u
	sol(1,0) = 0.05*(sin((x-0.5)*(x-0.5) + a*(y-0.5)*(y-0.5) + (z-0.5)*(z-0.5))+1)*exp(-t); //v
	sol(2,0) = 0.05*(sin((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) + a*(z-0.5)*(z-0.5))+1)*exp(-t); //w
	
	return sol;
}

hessian_variable Test4::Pressure(double _x, double _y, double _z, double _t) const
{
	const real pi = 3.1415926536;
	unknown x(_x,0), y(_y,1), z(_z,2), t(_t,3);
	return (x*x*x*y*y*z + x*sin(2*pi*x*z)*sin(2*pi*x*y)*sin(2*pi*z))*exp(-t);
}

rMatrix Test4::BCFlow_coef(double _x, double _y, double _z, double _t) const
{
	rMatrix ret(2,1);
	ret(0,0) = 1; //dirichlet
	ret(1,0) = 0;
	return ret;
}

rMatrix Test4::BCMech_coef(double _x, double _y, double _z, double _t) const
{
	rMatrix ret(4,1);
	ret(0,0) = 1; //dirichlet perp
	ret(1,0) = 0;
	ret(2,0) = 1; //dirichlet parallel
	ret(3,0) = 0;
	return ret;
}










