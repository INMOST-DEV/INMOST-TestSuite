#ifndef _TEST3_
#define _TEST3_

#include "inmost.h"
#include "abstract_test.h"


class Test3 : public AbstractTest
{
	INMOST::rMatrix E0, K0, B0, G0, d0;
	INMOST::rMatrix E1, K1, B1, G1, d1;
	double p1, p2;
	INMOST::rMatrix gp1,gp2, u1, u2, gu1,gu2, gu2t;
public:
	Test3();
	//media characteristics, tensors
	INMOST::vMatrix ElasticTensor(double x, double y, double z) const; //matrix 6 x 6
	INMOST::vMatrix BiotTensor(double x, double y, double z) const; //matrix 3 x 3
	INMOST::vMatrix PermTensor(double x, double y, double z) const; //matrix 3 x 3
	//media characteristics, scalars
	double Porosity(double x, double y, double z) const;
	double InverseBiotModulus(double x, double y, double z) const;
	double FluidDensity() const;
	double FluidViscosity() const;
	double SolidDensity() const;
	double Gravity() const;
	//solution
	INMOST::hMatrix Displacement(double x, double y, double z, double t) const; //vector 3 x 1: u,v,w
	INMOST::hessian_variable Pressure(double x, double y, double z, double t) const; // p
	//boundary condition coefficients
	INMOST::rMatrix BCMech_coef(double x, double y, double z, double t) const;
	INMOST::rMatrix BCFlow_coef(double x, double y, double z, double t) const;
	INMOST::MarkerType DarcyDomain() const {return 0;}
	bool Analytic() const {return true;}
	int Number() const {return 3;}
	std::string ProblemName() const {return "Linear analytic heterogeneous";}
};



#endif //_TEST1_
