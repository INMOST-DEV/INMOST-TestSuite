#ifndef _TEST8_
#define _TEST8_

#include "inmost.h"
#include "abstract_test.h"


class Test8 : public AbstractTest
{
	const double a; //Size of the box in x direction
	const double b; //Size of the box in y direction
	const double E; //Young modulus
	const double nu; //Poisson ratio
	const double K; //Permeability coefficient
	const double alpha; //Biot coefficient
	const double M; //specific storage coefficient
	const double mu; //fluid viscosity
	const double F; // force
	const double lambda, G; // Lame coefficients
	const double Ku; //undrained bulk modulus
	const double B; //skempton's coefficient
	const double nuu; // undrained Poisson ratio
	const double c; // consolidation coefficient
	
	INMOST::rMatrix C0, K0, B0;
public:
	Test8();
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
	int Number() const {return 6;}
	std::string ProblemName() const {return "Mandel";}
	
	void SetBC(INMOST::Mesh & m, double T, INMOST::MarkerType boundary) const;
	
	void SetForce(INMOST::Mesh& m, const INMOST::dynamic_variable& p, double T) const;
};



#endif //_TEST8_
