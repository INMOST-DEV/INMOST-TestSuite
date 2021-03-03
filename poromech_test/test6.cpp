#include "test6.h"

using namespace INMOST;

//shortcut for types
typedef Storage::real real;
typedef Storage::enumerator enumerator;
typedef Storage::var_array var_array;
typedef Storage::real_array real_array;
typedef Storage::reference_array ref_array;

const int nqmax = 250;

double Test6::Gravity() const {return 0.0;}
double Test6::FluidDensity() const {return 1;}
double Test6::FluidViscosity() const {return 1;}
double Test6::SolidDensity() const {return 1;}
double Test6::Porosity(double _x, double _y, double _z) const {return 1;}
double Test6::InverseBiotModulus(double _x, double _y, double _z) const {return c0;}

Test6::Test6()
:
E(1.0e+3),
nu(0.1),
a(1),b(1),
K(1.0e-2),
alpha(1),
x0(0.25),y0(0.25),
c0(0.0)
{
	lambda = E*nu/(1+nu)/(1-2*nu);
	mu = E/(1+nu)/2.0;
	beta = (lambda+2*mu)*K/a/b;
	std::cout << "lambda " << lambda << " mu " << mu << " beta " << beta << std::endl;
	B0 = rMatrix::Unit(3)*alpha;
	K0 = rMatrix::Unit(3)*K;
	C0 = GenIsotropicTensor(E,nu);
}
const double pi = 3.1415926535897932384626433832795;
double Test6::p_tilde(int n, int q, double T) const
{
	double gamman = n*pi;
	double gammaq = q*pi;
	double gammanq = gamman*gamman+gammaq*gammaq;
	double div = gammanq*gammanq+1;
	double ss = sin(gamman*x0)*sin(gammaq*y0);
	double s1 = (1.0/(gammanq + 1.0/gammanq))*sin(beta*T);
	double c2 = cos(beta*T)/div;
	double e3 = exp(-gammanq*beta*T)/div;
	std::cout << " n " << n << " q " << q << " ss " << ss << " s1 " << s1 << " c2 " << c2 << " e3 " << e3 << " div " << div << " cs1 " << (1.0/(gammanq + 1.0/gammanq)) << std::endl;
	return -2*ss*(s1-c2+e3);
}

double Test6::u_tilde(int n, int q, double T) const
{
	double gamman = n*pi;
	double gammanq = (n*n+q*q)*pi*pi;
	return gamman/gammanq*p_tilde(n,q,T);
}

double Test6::v_tilde(int n, int q, double T) const
{
	double gammaq = q*pi;
	double gammanq = (n*n+q*q)*pi*pi;
	return gammaq/gammanq*p_tilde(n,q,T);
}

vMatrix Test6::ElasticTensor(double _x, double _y, double _z) const {return C0;}
vMatrix Test6::BiotTensor(double _x, double _y, double _z) const { return B0; }
vMatrix Test6::PermTensor(double _x, double _y, double _z) const {return K0;}

hMatrix Test6::Displacement(double _x, double _y, double _z, double _t) const
{
	unknown x(_x,0), y(_y,1), z(_z,2), t(_t,3);
	hMatrix sol(3,1);
	sol.Zero();
	bool stop = false;
	
	for(int n = 1; n < nqmax && !stop; ++n)
	{
		stop = true;
		for(int q = 1; q < nqmax; q++)
		{
			double gamman = n*pi;
			double gammaq = q*pi;
			double u = u_tilde(n,q,_t);
			double v = v_tilde(n,q,_t);
			if( fabs(u) > 1.0e-9 && fabs(v) > 1.0e-9 ) stop = false;
			sol(0,0) += 4*u*cos(gamman*x)*sin(gammaq*y);
			sol(1,0) += 4*v*sin(gamman*x)*cos(gammaq*y);
			//std::cout << " n " << n << " q " << q << " u " << u << " v " << v << std::endl;
		}
		//if( stop ) std::cout << "at i " << i << std::endl;
	}
	return sol;
}

hessian_variable Test6::Pressure(double _x, double _y, double _z, double _t) const
{
	unknown x(_x,0), y(_y,1), z(_z,2), t(_t,3);
	hessian_variable sol = 0;
	bool stop = false;
	_t = pi/2.0/beta;
	
	for(int n = 1; n < nqmax && !stop; ++n)
	{
		stop = true;
		for(int q = 1; q < nqmax*4; q++)
		{
			double gamman = n*pi;
			double gammaq = q*pi;
			double p = p_tilde(n,q,_t);
			//if( fabs(p) > 1.0e-9 ) stop = false;
			sol -= 4*(lambda+2*mu)*p*sin(gamman*x)*sin(gammaq*y);
			std::cout << " n " << n << " q " << q << " p " << p << " sol " << get_value(sol) << " c " << 4*(lambda+2*mu) << std::endl;
		}
		//if( stop ) std::cout << "at i " << i << std::endl;
	}
	std::cout << "x " << _x << " y " << _y << " t " << _t << std::endl;
	return sol;
}

rMatrix Test6::BCFlow_coef(double _x, double _y, double _z, double _t) const
{
	rMatrix ret(2,1);
	ret(0,0) = 1;
	ret(1,0) = 0;
	return ret;
}

rMatrix Test6::BCMech_coef(double _x, double _y, double _z, double _t) const
{
	rMatrix ret(4,1);
	ret(0,0) = 0; //dirichlet perp
	ret(1,0) = 1; //neumann perp
	ret(2,0) = 1; //dirichlet parallel
	ret(3,0) = 0; //neumann parallel
	return ret;
}


void Test6::Init(Mesh & m)
{
	double ccnt[3] = {x0,y0,0.5};
	for(Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it) if( it->GetStatus() != Element::Ghost )
	{
		if( it->Inside(ccnt) )
		{
			source = it->self();
			double cnt[3];
			it->Centroid(cnt);
			std::cout << "proc " << m.GetProcessorRank();
			std::cout << " found source cell " << source.LocalID();
			std::cout << " at " << cnt[0] << " " << cnt[1] << " " << cnt[2];
			std::cout << std::endl;
		}
	}
}


void Test6::SetForce(Mesh & m, const INMOST::dynamic_variable & p, double T) const
{
	TagVariableArray tag_F = m.GetTag("FORCE");
	if(source.isValid())
		tag_F[source][3] = 2*beta*sin(beta*T)/source.Volume();
}





void Test6::SetBC(Mesh & m, double T, MarkerType boundary) const
{
	Automatizator * aut = Automatizator::GetCurrent();
	Automatizator::RemoveCurrent();
	TagRealArray tag_BC_flow = m.GetTag("BOUNDARY_CONDITION_FLOW");
	TagRealArray tag_BC_mech = m.GetTag("BOUNDARY_CONDITION_ELASTIC");
#if defined(USE_OMP)
#pragma omp parallel for
#endif
	for(int i = 0; i < m.FaceLastLocalID(); ++i) if( m.isValidFace(i) )
	{
		real cnt[3], nrm[3];
		Face f = m.FaceByLocalID(i);
		if( f.GetMarker(boundary) )
		{
			f.Barycenter(cnt);
			f.UnitNormal(nrm);
			if( fabs(fabs(nrm[2])-1) < 1.0e-6 ) // z-normal
			{
				tag_BC_flow[f][0] = 0;
				tag_BC_flow[f][1] = 1;
				tag_BC_flow[f][2] = 0;
				tag_BC_mech[f][0] = 1;
				tag_BC_mech[f][1] = 0;
				tag_BC_mech[f][2] = 0;
				tag_BC_mech[f][3] = 1;
				tag_BC_mech[f][4] = 0;
				tag_BC_mech[f][5] = 0;
				tag_BC_mech[f][6] = 0;
			}
			else
			{
				tag_BC_flow(f,3,1)(0,2,0,1) = BCFlow_coef(cnt[0],cnt[1],cnt[2],T);
				tag_BC_flow(f,3,1)(2,3,0,1).Zero();
				tag_BC_mech(f,7,1)(0,4,0,1) = BCMech_coef(cnt[0],cnt[1],cnt[2],T);
				tag_BC_mech(f,7,1)(4,7,0,1).Zero();
			}
		}
	}
	Automatizator::MakeCurrent(aut);
}




