#include "test6.h"

// [1] Barry, Mercer 
// Exact Solutions for Two Dimensional Time-Dependent 
// Flow and Deformation Within a Poroelastic Medium
//
// [2] Peichao Li et al
// Analytical solutions of a finite two-dimensional fluid-saturated
// poroelastic medium with compressible constituents

using namespace INMOST;


//shortcut for types
typedef Storage::real real;
typedef Storage::enumerator enumerator;
typedef Storage::var_array var_array;
typedef Storage::real_array real_array;
typedef Storage::reference_array ref_array;

const int nqmax = 5000;
const double pi = 3.1415926535897932384626433832795;
const double errtol = 1e-7;

double Test6::Gravity() const {return 0.0;}
double Test6::FluidDensity() const {return 1;}
double Test6::FluidViscosity() const {return 1;}
double Test6::SolidDensity() const {return 1;}
double Test6::Porosity(double _x, double _y, double _z) const {return 1;}
double Test6::InverseBiotModulus(double _x, double _y, double _z) const { return M ? 1.0/M : 0.0; }

Test6::Test6()
: E(1.0e+4), nu(0.25), a(1), b(1), K(0.001), alpha(1.0), aQ(2), 
  x0(0.25), y0(0.25), M(0), //q0(0.03), 
  //G(E/(2.0*(1+nu))),m(1.0/(1-2*nu)),lambda_f(K/mu),chi(lambda_f*M)
	lambda(E* nu / (1 + nu) / (1 - 2 * nu)), mu(E / (1 + nu) / 2.0),
	omega(1)
{
	//lambda = E*nu/(1+nu)/(1-2*nu);
	//mu = E/(1+nu)/2.0;
	//beta = (lambda+2*mu)*K/a/b;
	std::cout << "E: " << E << " nu: " << nu << std::endl;
	//std::cout << "G: " << G << " m " << m << std::endl;
	std::cout << "K: " << K << " M " << M << std::endl;
	//std::cout << "lambda_f: " << lambda_f << " chi " << chi << std::endl;
	std::cout << "alpha: " << alpha << std::endl;
	std::cout << "lambda " << lambda << " mu " << mu  << std::endl;
	std::cout << "(lambda+2*mu)*K:" << (lambda + 2 * mu) * K << std::endl;
	B0 = rMatrix::Unit(3)*alpha;
	K0 = rMatrix::Unit(3)*K;
	C0 = GenIsotropicTensor(E,nu);
}



vMatrix Test6::ElasticTensor(double _x, double _y, double _z) const {return C0;}
vMatrix Test6::BiotTensor(double _x, double _y, double _z) const { return B0; }
vMatrix Test6::PermTensor(double _x, double _y, double _z) const {return K0;}

hMatrix Test6::Displacement(double _x, double _y, double _z, double _t) const
{
	unknown x(_x,0), y(_y,1), z(_z,2), t(_t,3);
	rMatrix sol(3,1);
	double u_tilde, v_tilde, p_tilde;
	sol.Zero();
	//double t_coef0 = -1.0 / ( alpha * alpha / (G * lambda_f * (m + 1)) + (chi ? 1.0/chi : 0.0));
	double t_coef0 = (lambda + 2 * mu) * K;
	int stop = 3;
	int n = 1, q = 1;
	while(n+q < std::min(stop,nqmax) )
	{
		double lambda_n = n*pi;
		double lambda_q = q*pi;
		double lambda_nq = (n * n / a + q * q / b) * pi * pi;
		//double uv_coef = alpha * q0 * sin(lambda_n * x0) * sin(lambda_q * y0)
		//	/ (G * lambda_f * (m + 1) * lambda_nq * lambda_nq);
		//double u_coef = uv_coef * lambda_n;
		//double v_coef = uv_coef * lambda_q;
		//double t_coef = t_coef0 * lambda_nq;
		//u_tilde = u_coef * (1.0 - exp(t_coef * _t));
		//v_tilde = v_coef * (1.0 - exp(t_coef * _t));
		p_tilde = aQ*sin(lambda_n * x0) * sin(lambda_q * y0) 
			/ (t_coef0 * (lambda_nq * lambda_nq + omega * omega))
			* (lambda_nq * sin(omega * t_coef0 * _t) + omega * (exp(-lambda_nq * t_coef0 * _t) - cos(omega * t_coef0 * _t)));
		u_tilde = -lambda_n * p_tilde / lambda_nq;
		v_tilde = -lambda_q * p_tilde / lambda_nq;
		//std::cout << "n=" << n << " q=" << q << " ut=" << get_value(u_tilde) << " vt=" << get_value(v_tilde) << std::endl;
		if (std::fabs(get_value(u_tilde)) > errtol * std::fabs(get_value(sol(0, 0))) &&
			std::fabs(get_value(v_tilde)) > errtol * std::fabs(get_value(sol(1, 0))))
			stop = n+q+2;
		sol(0,0) += u_tilde*cos(lambda_n*_x)*sin(lambda_q*_y);
		sol(1,0) += v_tilde*sin(lambda_n*_x)*cos(lambda_q*_y);
		n++; q--;
		if (q == 0) { q = n; n = 1; }
	}
//#pragma omp critical
//	std::cout << " uv stop at n=" << n << " q=" << q << " stop=" << stop << std::endl;
	if (sol.CheckNansInfs())
	{
		std::cout << "Bad sol" << std::endl;
		sol.Print();
	}
	return 4.0 / (a * b) * sol;
}

hessian_variable Test6::Pressure(double _x, double _y, double _z, double _t) const
{
	unknown x(_x,0), y(_y,1), z(_z,2), t(_t,3);
	double sol = 0;
	double p_tilde;
	//double t_coef0 = -1.0 / (alpha * alpha / (G * lambda_f * (m + 1)) + (chi ? 1.0 / chi : 0.0));
	double t_coef0 = (lambda + 2 * mu) * K;
	int stop = 3;
	int n = 1, q = 1;
	while (n+q < std::min(stop,nqmax) )
	{
		double lambda_n = n*pi/a;
		double lambda_q = q*pi/b;
		double lambda_nq = (n * n / a + q * q / b) * pi * pi;
		//double p_coef = - q0 * sin(lambda_n * x0) * sin(lambda_q * y0) / (lambda_f * lambda_nq);
		//double t_coef = t_coef0*lambda_nq;
		//p_tilde = p_coef*(1.0 - exp(t_coef * _t));
		p_tilde = aQ*sin(lambda_n * x0) * sin(lambda_q * y0) 
			/ (t_coef0 * (lambda_nq * lambda_nq + omega * omega))
			* (lambda_nq * sin(omega * t_coef0 * _t) - omega * cos(omega * t_coef0 * _t) + omega * exp(-lambda_nq * t_coef0 * _t));

		//std::cout << "n=" << std::setw(4) << n << " q=" << std::setw(4) << q << " pt=" << std::setw(14) << get_value(p_tilde) << " sol=" << std::setw(14) << get_value(sol) << " x=" << _x << " y=" << _y << std::endl;
		if ( std::fabs(get_value(p_tilde)) > errtol * std::fabs(get_value(sol)))
			stop = n+q+2;
		sol += p_tilde * sin(lambda_n * _x) * sin(lambda_q * _y);
		n++; q--;
		if (q == 0) { q = n; n = 1; }
	}
//#pragma omp critical
//	std::cout << " p stop at n=" << n << " q=" << q << " stop=" << stop << std::endl;
	if (check_nans_infs(sol))
	{
		std::cout << "Bad sol " << get_value(sol) << std::endl;
	}
	return 4.0 / (a * b) * sol * (lambda + 2 * mu);
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
	double t_coef0 = (lambda + 2 * mu) * K;
	if (source.isValid())
		//tag_F[source][3] = -q0 / source.Volume();
		tag_F[source][3] = aQ*sin(omega*t_coef0*T) / source.Volume();
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




