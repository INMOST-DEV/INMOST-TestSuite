#include "test1.h"

using namespace INMOST;

//shortcut for types
typedef Storage::real real;
typedef Storage::enumerator enumerator;
typedef Storage::var_array var_array;
typedef Storage::real_array real_array;
typedef Storage::reference_array ref_array;

double Test1::Gravity() const {return 1.0;}
double Test1::FluidDensity() const {return 1000;}
double Test1::FluidViscosity() const {return 0.1;}
double Test1::SolidDensity() const {return 2000;}
double Test1::Porosity(double _x, double _y, double _z) const {return 0.2;}
double Test1::InverseBiotModulus(double _x, double _y, double _z) const {return 10;}

Test1::Test1()// : C0(GenMatrix1(6)), K0(GenMatrix1(3)), B0(GenMatrix1(3))
{
	real vC0[36] =
	{
		81,           33,           5,            14,           37,           32,
		33,           100,          26,           58,           58,           73,
		5,            26,           90,           23,           60,           28,
		14,           58,           23,           77,           29,           3,
		37,           58,           60,           29,           84,           60,
		32,           73,           28,            3,           60,           99
		
	};
	real vC1[36] =
	{
		93,           46,           22,           13,           72,           35,
		46,           95,           41,           62,           56,           24,
		22,           41,           89,           25,           33,           21,
		13,           62,           25,           87,           13,           25,
		72,           56,           33,           13,           99,           57,
		35,           24,           21,           25,           57,           78
	};
	real vK0[9] =
	{
		1,            2,            2,
		2,           68,           42,
		2,           42,           49
	};
	real vK1[9] =
	{
		25,            2,           39,
		2,           42,            7,
		39,            7,          100
	};
	real vB0[9] =
	{
		1,            1,            1,
		1,           68,           13,
		1,           13,           73
	};
	real vB1[9] =
	{
		1,            6,            5,
		6,           67,           27,
		5,           27,           76
	};
	C0 = rMatrix(vC1,6,6);
	K0 = rMatrix(vK1,3,3)/1000.0;
	B0 = rMatrix(vB1,3,3)/100.0;
	//int seed = 0;
	//GenTensor(C0,6,seed);
	//seed = GenTensor(K0,3,seed);
	//seed = GenTensor(B0,3,seed+1);
	std::cout << "C0: " << std::endl;
	C0.Print();
	rCheckEigen(C0,true);
	std::cout << "K0: " << std::endl;
	K0.Print();
	rCheckEigen(K0,true);
	std::cout << "B0: " << std::endl;
	B0.Print();
	rCheckEigen(B0,true);
}

vMatrix Test1::ElasticTensor(double _x, double _y, double _z) const {return C0;}
vMatrix Test1::BiotTensor(double _x, double _y, double _z) const { return B0; }
vMatrix Test1::PermTensor(double _x, double _y, double _z) const
{
	return K0;
	unknown x(_x,0), y(_y,1), z(_z,2);
	vMatrix K(3,3);
	/*
	real a = 0.001;
	variable div = 1;//(x*x+y*y+z*z+1.0e-5);
	K(0,0) = (a*x*x + z*z + y*y + 1)/div;
	K(0,1) = K(1,0) = (a-1)*x*y/div;
	K(0,2) = K(2,0) = (a-1)*x*z/div;
	K(1,1) = (a*y*y + x*x + z*z + 1)/div;
	K(1,2) = K(2,1) = (a-1)*y*z/div;
	K(2,2) = (a*z*z + x*x + y*y + 1)/div;
	CheckEigen(K);
	 */
	K = rMatrix::Unit(3);
	return K;
}

hMatrix Test1::Displacement(double _x, double _y, double _z, double _t) const
{
	unknown x(_x,0), y(_y,1), z(_z,2), t(_t,3);
	hMatrix sol(3,1);
	/*
	sol(0,0) = (sin(((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) + (z-0.5)*(z-0.5))*8)+1)*exp(-t); //u
	sol(1,0) = (sin(((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) + (z-0.5)*(z-0.5))*8)+1)*exp(-t); //v
	sol(2,0) = (sin(((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) + (z-0.5)*(z-0.5))*8)+1)*exp(-t); //w
	 */
	
	
	sol(0,0) = 1*x+2*y+3*z+4*t+5;
	sol(1,0) = 6*x+7*y+8*z+9*t+10;
	sol(2,0) = 11*x+12*y+13*z+14*t+15;
	 
	
	/*
	sol(0,0) = 1*x+2*y+3*z+5;
	sol(1,0) = 6*x+7*y+8*z+10;
	sol(2,0) = 11*x+12*y+13*z+15;
	*/
	
	/*
	sol(0,0) = 1*x+2*y+3*z+1;
	sol(1,0) = 4*x+3*y+2*z+1;
	sol(2,0) = 5*x+6*y+7*z+1;
	*/
	
	return sol;
}

hessian_variable Test1::Pressure(double _x, double _y, double _z, double _t) const
{
	const real pi = 3.1415926536;
	unknown x(_x,0), y(_y,1), z(_z,2), t(_t,3);
	//return 19*t+20;
	return 16*x+17*y+18*z+19*t+20;
	//return 16*x+17*y+18*z+20;
	//return 9*x+10*y+11*z+1;
	//return (x*x*x*y*y*z + x*sin(2*pi*x*z)*sin(2*pi*x*y)*sin(2*pi*z))*exp(-t);
}

rMatrix Test1::BCFlow_coef(double _x, double _y, double _z, double _t) const
{
	rMatrix ret(2,1);
	if( 1 - _x - _y - _z > 0 )
	{
		ret(0,0) = 1; //dirichlet parallel
		ret(1,0) = 0; //neumann parallel
		//ret(0,0) = 1; //dirichlet parallel
		//ret(1,0) = 0; //neumann parallel
	}
	else
	{
		ret(0,0) = 0; //dirichlet
		ret(1,0) = 1;
	}
	return ret;
}

rMatrix Test1::BCMech_coef(double _x, double _y, double _z, double _t) const
{
	rMatrix ret(4,1);
	if( 1 - _x - _y - _z > 0 )
	{
		ret(0,0) = 0; //dirichlet perp
		ret(1,0) = 1; //neumann perp
		ret(2,0) = 0; //dirichlet parallel
		ret(3,0) = 1; //neumann parallel
		//ret(0,0) = 1; //dirichlet perp
		//ret(1,0) = 0; //neumann perp
		//ret(2,0) = 1; //dirichlet parallel
		//ret(3,0) = 0; //neumann parallel
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










