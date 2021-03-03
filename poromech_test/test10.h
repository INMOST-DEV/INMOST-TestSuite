#ifndef _TEST10_
#define _TEST10_

#include "inmost.h"
#include "abstract_test.h"


class Test10 : public AbstractTest
{
	double WI[3], pbhp[3];
	INMOST::Cell cc[3];
public:
	Test10();
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
	bool Analytic() const {return false;}
	int Number() const {return 10;}
	std::string ProblemName() const {return "Norne field";}
	
	void Init(INMOST::Mesh & m);
	void SetBC(INMOST::Mesh & m, double T, INMOST::MarkerType boundary) const;
	void SetForce(INMOST::Mesh & m, const INMOST::dynamic_variable & p,  double T) const;
	void SetProperty(INMOST::Mesh & m) const;
};



#endif //_TEST10_
