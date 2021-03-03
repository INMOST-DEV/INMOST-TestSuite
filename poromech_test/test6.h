#ifndef _TEST6_
#define _TEST6_

#include "inmost.h"
#include "abstract_test.h"


class Test6 : public AbstractTest
{
	double a; //Size of the box in x direction
	double b; //Size of the box in y direction
	double E; //Young modulus
	double nu; //Poisson ratio
	double alpha; //Biot coefficient
	double K; //Permeability coefficient
	double x0, y0; //source position
	double c0; //specific storage coefficient
	double lambda, mu; //lame coefficients
	double beta;
	INMOST::rMatrix C0, K0, B0;
	INMOST::Cell source;
public:
	Test6();
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
	std::string ProblemName() const {return "Barry-Mercer";}
	
	void Init(INMOST::Mesh & m);
	void SetBC(INMOST::Mesh & m, double T, INMOST::MarkerType boundary) const;
	void SetForce(INMOST::Mesh & m, const INMOST::dynamic_variable & p,  double T) const;
	
	double p_tilde(int n, int q, double T) const;
	double u_tilde(int n, int q, double T) const;
	double v_tilde(int n, int q, double T) const;
	
};



#endif //_TEST6_
