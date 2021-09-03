#ifndef _TEST9_
#define _TEST9_

#include "inmost.h"
#include "abstract_test.h"


class Test9 : public AbstractTest
{
	const double a1; //Size of the top layer in x direction
	const double a2; //Size of the bottom layer in x direction
	const double b; //Size of the box in y direction
	const double E1,E2; //Young modulus
	const double nu1,nu2; //Poisson ratio
	const double k1,k2; //Permeability coefficient
	const double alpha1,alpha2; //Biot coefficient
	const double M1; //specific storage coefficient
	const double mu; //fluid viscosity
	const double F; // force
	//dependent parameters
	const double S1;
	const double lambda1, G1, lambda2, G2; // Lame coefficients
	const double m1, m2;
	const double M2,S2;
	//const double Ku1, Ku2; //undrained bulk modulus
	const double B1,B2; //skempton's coefficient
	//const double nuu1,nuu2; // undrained Poisson ratio
	const double c1,c2; // consolidation coefficient
	const double theta;
	const double beta;
	INMOST::rMatrix C1, K1, bB1, C2, K2, bB2;
	std::vector<double> wi; //trigonometric equation roots
public:
	Test9();
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
	std::string ProblemName() const {return "Terzaghi double-layered";}
	
	void SetBC(INMOST::Mesh & m, double T, INMOST::MarkerType boundary) const;
	void Init(INMOST::Mesh& m);
	void SetForce(INMOST::Mesh& m, const INMOST::dynamic_variable& p, double T) const;

	double Qfunc(double wi, double t) const;
	double P1func(double wi, double x, double t) const;
	double intP1func(double wi, double x, double t) const;
	double P2func(double wi, double x, double t) const;
	double intP2func(double wi, double x, double t) const;
};



#endif //_TEST9_
