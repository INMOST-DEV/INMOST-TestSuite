#ifndef _TEST13_
#define _TEST13_

#include "inmost.h"
#include "abstract_test.h"


class Test13 : public AbstractTest
{
	double WI[2], pbhp[2];
	double fWI[1], fpbhp[1];
	INMOST::Cell cc[2];
	INMOST::Face fc[1];
	INMOST::rMatrix E0,E1,K0,K1,B0,B1,Kf;
	double phi0,phi1, M0, M1, Mf;
	INMOST::MarkerType darcy;
public:
	Test13();
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
	INMOST::MarkerType DarcyDomain() const {return darcy;}
	bool Analytic() const {return false;}
	int Number() const {return 13;}
	
	void Init(INMOST::Mesh & m);
	void SetBC(INMOST::Mesh & m, double T, INMOST::MarkerType boundary) const;
	void SetForce(INMOST::Mesh & m,const INMOST::dynamic_variable & p,  double T) const;
	void SetForceDarcy(INMOST::Mesh & m,const INMOST::dynamic_variable & p,  double T) const;
	void SetForceDarcyFracture(INMOST::Mesh & m,const INMOST::dynamic_variable & p,const INMOST::dynamic_variable & af,  double T) const;
	void SetProperty(INMOST::Mesh & m) const;
	
	std::string ProblemName() const {return "Poromechanics-darcy";}
};



#endif //_TEST13_
