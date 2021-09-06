#ifndef _TEST12_
#define _TEST12_

#include "inmost.h"
#include "abstract_test.h"


class Test12 : public AbstractTest
{
	const double a; //Size of the box in x direction
	const double b; //Size of the box in y direction
	const double c; //Size of the box in z direction
	const double k; //fluid permeability
	const double mu; //fluid viscosity
	const double Ks, Kf; // solid and fluid bulk modulus
	const double phi; // porosity
	const double F; // force
	double M; //specific storage coefficient
	double A11, A12, A13, A14, A15;
	double B11, B12, B13;
	double p0; //initial pressure
	double c1, c2, c3; //consolidation coefficients
	double Dmcoef; //coefficient in series
	double lambda0;
	const bool forcebc; // set top force or displacement
	INMOST::rMatrix C, K, B;
public:
	Test12();
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
	std::string ProblemName() const {return "Anisotropic Mandel";}
	
	void SetBC(INMOST::Mesh & m, double T, INMOST::MarkerType boundary) const;
	void Init(INMOST::Mesh& m);
	void SetForce(INMOST::Mesh& m, const INMOST::dynamic_variable& p, double T) const;
};



#endif //_TEST12_
