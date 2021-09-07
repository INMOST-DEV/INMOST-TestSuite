#ifndef _TEST12_
#define _TEST12_

#include "inmost.h"
#include "abstract_test.h"

// [1] Y. ABOUSLEIMAN, A. H.-D. CHENG, L. CUI, E. DETOURNAY and J.-C. ROEGIERS
// Mandel's problem revisited

class Test12 : public AbstractTest
{
	const double a; //Size of the box in x-y direction
	const double b; //Size of the box in z direction
	const double Ex, Ey, Ez; //Young modulus
	const double nuxy, nuxz, nuyz; // Poisson ratio
	const double Gxy, Gxz, Gyz; // Shear moduli
	const double k1, k2, k3; //fluid permeability
	const double mu; //fluid viscosity
	const double Ks, Kf; // solid and fluid bulk modulus
	const double phi; // porosity
	const double nuzx, nuyx;
	const double M11, M12, M13, M33, M55; // formulas 3(a-e) in [1]
	const double alpha1, alpha2, alpha3; //Biot modulus
	const double M; //specific storage coefficient
	const double c1; //consolidation coefficient, formula (21) in [1]
	const double A1, A2; //formulas (25) and (30) in [1]
	const double F; // force
	const bool forcebc; // set top force or displacement
	INMOST::rMatrix C, K, B;
	std::vector<double> wi; //trigonometric equation roots
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
