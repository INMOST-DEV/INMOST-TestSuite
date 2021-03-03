#ifndef _ATEST_
#define _ATEST_

#include "inmost.h"
//Tests
// test1 : analytic linear solution
// test2 : analytic nonlinear solution
// test3 : analytic steady-state with discontinuous tensors
// test4 : analytic solution with heterogeneous tensors defined by coordinates
// test5 : Terzaghi column (not implemented)
// test6 : Barry and Mercer (not implemented) 
// test7 : analytic anisotropic locking issue with neumann and dirichlet bc
// test8 : sphere with infinite tensor in the middle (not implemented)
// test9 : sector - real experement from paper of Nicola (not implemented)
// test10 : norne field
// test11 : two wells stability
// test12 : mandel problem (not implemented)
// test13 : two wells, left part - poromech, right part - darcy


class AbstractTest
{
public:
	bool bdf2;
	bool staggered;
	bool gravity_in_flux;
	bool first_order_flux;
	bool cranknicholson;
	bool func_output;
	AbstractTest();
	//media characteristics, tensors
	virtual INMOST::vMatrix ElasticTensor(double x, double y, double z) const = 0; //matrix 6 x 6
	virtual INMOST::vMatrix BiotTensor(double x, double y, double z) const = 0; //matrix 3 x 3
	virtual INMOST::vMatrix PermTensor(double x, double y, double z) const = 0; //matrix 3 x 3
	//media characteristics, scalars
	virtual double Porosity(double x, double y, double z) const = 0;
	virtual double InverseBiotModulus(double x, double y, double z) const = 0;
	virtual double FluidDensity() const = 0;
	virtual double FluidViscosity() const = 0;
	virtual double SolidDensity() const = 0;
	virtual double Gravity() const = 0;
	//solution
	virtual INMOST::hMatrix Displacement(double x, double y, double z, double t) const = 0; //vector 3 x 1: u,v,w
	virtual INMOST::hessian_variable Pressure(double x, double y, double z, double t) const = 0; // p
	//boundary condition coefficients
	virtual INMOST::rMatrix BCMech_coef(double x, double y, double z, double t) const = 0;
	virtual INMOST::rMatrix BCFlow_coef(double x, double y, double z, double t) const = 0;
	//boundary conditions
	INMOST::rMatrix BCMech(double x, double y, double z, double t, double nrm[3]) const;
	INMOST::rMatrix BCFlow(double x, double y, double z, double t, double nrm[3]) const;
	//force and flux
	INMOST::rMatrix Force(double x, double y, double z, double t) const;
	INMOST::rMatrix Flux(double x, double y, double z, double t, double nrm[3], bool bndcond=false) const;
	//mechanical stress, 6x1 vector in voigt notation
	// full - includes biot term
	INMOST::rMatrix Stress(double x, double y, double z, double t, bool full = false) const;
	virtual bool Analytic() const = 0;
	virtual int Number() const = 0;
	virtual ~AbstractTest() {}
	virtual void Init(INMOST::Mesh & m) {}
	virtual void SetProperty(INMOST::Mesh & m) const;
	virtual void SetInitial(INMOST::Mesh & m, double T, double Told, INMOST::MarkerType orient) const;
	virtual void SetForce(INMOST::Mesh & m, const INMOST::dynamic_variable & p, double T) const;
	virtual void SetForceDarcy(INMOST::Mesh & m, const INMOST::dynamic_variable & p, double T) const {}
	virtual void SetForceDarcyFracture(INMOST::Mesh & m, const INMOST::dynamic_variable & p, const INMOST::dynamic_variable & af, double T) const {}
	virtual void SetBC(INMOST::Mesh & m, double T, INMOST::MarkerType boundary) const;
	virtual INMOST::MarkerType DarcyDomain() const = 0;
	virtual std::string ProblemName() const = 0;
};

void Report(const AbstractTest * test, INMOST::Mesh & m, double T, double dT, double time_shift,INMOST::MarkerType orient);
AbstractTest * MakeTest(int number);

int GenTensor(INMOST::rMatrix & C, int size, int seed, bool print = false);
INMOST::rMatrix GenMatrix1(int size);
bool rCheckEigen(const INMOST::rMatrix & A, bool print=true);
bool vCheckEigen(const INMOST::vMatrix & A, bool print=true);

INMOST::rMatrix GenIsotropicTensor(double E, double nu);
INMOST::rMatrix GenAnisotropicTensor(double E1, double E2, double E3, double nu12, double nu13, double nu23, double G12, double G23, double G13);

/// Decomposes elastic tensor into 9 co-normal 3x3 tensors and computes
/// a matrix of projections of co-normal tensors onto normal component
/// and respective tangential parts.
void SplitConormals(const INMOST::Storage::real_array & E,
					const INMOST::rMatrix & n,
					INMOST::rMatrix & L, // 3x3
					INMOST::rMatrix & T, // 3x12
					bool transversal);

void varPrint(const INMOST::variable & v, double tol = 0);
void vMatrixPrint(const INMOST::vMatrix & m, double tol = 0);

void SVD2Eigen(const INMOST::rMatrix & U, INMOST::rMatrix & S, INMOST::rMatrix & V);

#endif //_ATEST_
