#include "test5.h"

// [1] Herminio T.Honorio, Clovis R.Maliska, Massimiliano Ferronato, Carlo Janna
// A stabilized element-based finite volume method for poroelastic problems

using namespace INMOST;

const int istop = 1000000;
const double epsuv = 1.0e-11;
const double epsp = 1.0e-11;
const double epsnln = 1.0e-13;
//shortcut for types
typedef Storage::real real;
typedef Storage::enumerator enumerator;
typedef Storage::var_array var_array;
typedef Storage::real_array real_array;
typedef Storage::reference_array ref_array;

const double pi = 3.1415926535897932384626433832795;
const double errtol = 1e-7;

double Test5::Gravity() const {return 0.0;}
double Test5::FluidDensity() const {return 998.2;}
double Test5::FluidViscosity() const {return mu;}
double Test5::SolidDensity() const {return 2000;}
double Test5::Porosity(double _x, double _y, double _z) const {return 0.3;}
double Test5::InverseBiotModulus(double _x, double _y, double _z) const { return M ? 1.0/M : 0.0; }

Test5::Test5()
: 
	// M = (phi*cf+(alpha-phi)*cs)^{-1}
	// cs = 0
	// cf = 4.59 x 10^{-4}
	// phi = 0.3
	// M = (4.59 * 0.4 x 10^{-4})^{-1} = 10^4 / (4.59*0.4) = 5446,623
	// G = 0.819
	// lambda = 1.227
	// E = G*(3*lambda+2*G)/(lambda+G) 
	// nu = lambda/(2*(lambda+G))
	//a(1), b(1), E(1.0e+3), nu(0.25),
	//K(0.01), alpha(1.0), M(0.1), mu(1), F(10),
	a(1), b(2), 
	//E(2.129159824046921), nu(0.2998533724340176),
	E(2.5), nu(0.45),
	//K(9.98e-10), 
	K(1.0e-9),
	alpha(1.0), 
	//M(5446.623), 
	M(10000),
	//mu(1.002e-3), 
	mu(0.01),
	F(10),
	forcebc(true),
	//derived parameters
	lambda(E* nu / (1 + nu) / (1 - 2 * nu)), 
	G(E / (1 + nu) / 2.0),
	Ku(lambda + 2 * G / 3.0 + alpha * alpha * M),
	B(alpha*M/Ku),
	nuu( (3*nu + alpha*B*(1-2*nu))/(3 - alpha*B*(1-2*nu)) ),
	c( 2*K*G*(1-nu)*(nuu-nu)/(mu*alpha*alpha*(1-nuu)*(1-2*nu)*(1-2*nu)))
{
	std::cout << "E: " << E << " nu: " << nu << std::endl;
	std::cout << "K: " << K << " M " << M << std::endl;
	std::cout << "alpha: " << alpha << std::endl;
	std::cout << "lambda " << lambda << " G " << G  << std::endl;
	std::cout << "consolidation coefficient: " << c << std::endl;
	std::cout << "undrained Poisson: " << nuu << std::endl;
	std::cout << "Skempton's coefficient: " << B << std::endl;
	//double B2 = 3.0 * (nuu - nu) / (alpha * (1 - 2 * nu) * (1 + nuu));
	//double c2 = 2.0 * K * B2 * B2 * G * (1 - nu) * (1 + nuu) * (1 + nuu) / (9.0 * (1 - nuu) * (nuu - nu));
	//std::cout << "Skempton's 2: " << B2 << " consolidation 2: " << c2 << std::endl;
	//std::cout << "(lambda+2*mu)*K:" << (lambda + 2 * mu) * K << std::endl;
	B0 = rMatrix::Unit(3)*alpha;
	K0 = rMatrix::Unit(3)*K;
	C0 = GenIsotropicTensor(E,nu);
}



vMatrix Test5::ElasticTensor(double _x, double _y, double _z) const {return C0;}
vMatrix Test5::BiotTensor(double _x, double _y, double _z) const { return B0; }
vMatrix Test5::PermTensor(double _x, double _y, double _z) const {return K0;}

static double solvetri(double coef, double x)
{
	// j dx = -r, dx = -r/j, x += dx, x -= r/j
//#pragma omp critical
	{
		double r, j;
		int nit = 0;
		do
		{
			r = cos(x) - coef * sin(x) / x;
			j = coef * sin(x) / x / x - coef * cos(x) / x - sin(x);
			x -= r / j;
			//std::cout << nit << " x " << x << " r " << r << " j " << j << " coef " << coef << std::endl;
			nit++;
		} while (std::fabs(r) > epsnln*std::fabs(x));
	}
	return x;
}

hMatrix Test5::Displacement(double _x, double _y, double _z, double _t) const
{
	unknown x(_x,0), y(_y,1), z(_z,2), t(_t,3);
	hMatrix sol(3,1);
	double du, dv;
	const double cnuu = (nuu - nu) / (1 - nu);
	double wi;
	sol.Zero();
	int stop = -1;
	if (_t)
	{
		sol(0, 0) = F / (G * a) * nu * _x / 2.0;
		sol(1, 0) = F / (G * a) * (nu - 1) * _y / 2.0;
		for (int i = 0; i < istop; ++i)
		{
			//cos(wi) = cnuu * sin(wi) / wi
			wi = solvetri(cnuu, 0.5 * pi + pi * i);
			du = F / (G * a) * (a * sin(wi * _x / a) - nuu * _x * sin(wi)) * cos(wi) / (wi - sin(wi) * cos(wi)) * exp(-wi * wi * c * _t / a / a);
			dv = F / (G * a) * (1 - nuu) * _y * sin(wi) * cos(wi) / (wi - sin(wi) * cos(wi)) * exp(-wi * wi * c * _t / a / a);
			sol(0, 0) += du;
			sol(1, 0) += dv;
			if (std::fabs(get_value(du)) < epsuv * std::fabs(get_value(sol(0, 0))) &&
				std::fabs(get_value(dv)) < epsuv * std::fabs(get_value(sol(1, 0))))
				break;
		}
	}
	else
	{
		sol(0, 0) = F / (G * a) * nuu * _x / 2.0;
		sol(1, 0) = F / (G * a) * (nuu - 1) * _y / 2.0;
	}
	return sol;
}

hessian_variable Test5::Pressure(double _x, double _y, double _z, double _t) const
{
	unknown x(_x,0), y(_y,1), z(_z,2), t(_t,3);
	double sol = 0, dp;
	const double cnuu = (nuu - nu) / (1 - nu);
	const double p0 = F * B * (1 + nuu) / (3.0 * a);
	if (_t)
	{
		for (int i = 0; i < istop; ++i)
		{
			//cos(wi) = cnuu * sin(wi) / wi
			//sin(wi) = wi*cos(wi)/cnuu
			double wi = solvetri(cnuu, 0.5 * pi + pi * i);
			dp = 2.0 * p0 * sin(wi) / (wi - sin(wi) * cos(wi)) * (cos(wi * _x / a) - cos(wi)) * exp(-wi * wi * c * _t / a / a);
			sol += dp;
			if (std::fabs(get_value(dp)) < epsp * std::fabs(get_value(sol)))
				break;
		}
	}
	else sol = p0;
	//else sol = F * B * (1 + nuu) / (6.0 * a);
	return sol;
}

rMatrix Test5::BCFlow_coef(double _x, double _y, double _z, double _t) const
{
	rMatrix ret(2,1);
	//impermiable
	ret(0,0) = 0;
	ret(1,0) = 1;
	return ret;
}

rMatrix Test5::BCMech_coef(double _x, double _y, double _z, double _t) const
{
	rMatrix ret(4,1);
	//other sides - roller
	{
		ret(0, 0) = 1; //dirichlet perp
		ret(1, 0) = 0; //neumann perp
		ret(2, 0) = 0; //dirichlet parallel
		ret(3, 0) = 1; //neumann parallel
	}
	return ret;
}




void Test5::Init(Mesh& m)
{
	double cmin[3] = { +1.0e+20,+1.0e+20,+1.0e+20 }, cmax[3] = { -1.0e+20,-1.0e+20,-1.0e+20 };
	for (Mesh::iteratorNode it = m.BeginNode(); it != m.EndNode(); ++it) 
	{
		real_array& c = it->Coords();
		for (int k = 0; k < c.size(); ++k)
		{
			cmin[k] = std::min(cmin[k], c[k]);
			cmax[k] = std::max(cmax[k], c[k]);
		}
	}
	m.AggregateMax(cmax, 3);
	m.AggregateMin(cmin, 3);
	if (!m.GetProcessorRank())
		std::cout << "Mesh in " << cmin[0] << ":" << cmax[0] << " " << cmin[1] << ":" << cmax[1] << " " << cmin[2] << ":" << cmax[2] << std::endl;
	double s[3] = { a,b,a };
	for (Mesh::iteratorNode it = m.BeginNode(); it != m.EndNode(); ++it)
	{
		real_array& c = it->Coords();
		for (int k = 0; k < c.size(); ++k)
			c[k] = (c[k] - cmin[k]) / (cmax[k] - cmin[k]) * s[k];
	}
	//m.RecomputeGeometricData();
}

void Test5::SetBC(Mesh & m, double T, MarkerType boundary) const
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
			if( fabs(nrm[1]-1) < 1.0e-6 ) // y-normal up
			{
				//fluid
				//double P = Pressure(cnt[0], cnt[1], cnt[2], T).GetValue();
				tag_BC_flow[f][0] = 0;
				tag_BC_flow[f][1] = 1;
				tag_BC_flow[f][2] = 0;
				//mech
				if (forcebc) //uniform force 
				{
					tag_BC_mech[f][0] = 0; //D-perp
					tag_BC_mech[f][1] = 1; //N-perp
					tag_BC_mech[f][2] = 0; //D-parallel
					tag_BC_mech[f][3] = 1; //N-parallel
					tag_BC_mech[f][4] = -nrm[0] * F;
					tag_BC_mech[f][5] = -nrm[1] * F;
					tag_BC_mech[f][6] = -nrm[2] * F;
				}
				else //analytical displacement
				{
					rMatrix UVW = Displacement(cnt[0], cnt[1], cnt[2], T);
					double nUVW = nrm[0] * UVW(0, 0) + nrm[1] * UVW(1, 0) + nrm[2] * UVW(2, 0);
					tag_BC_mech[f][0] = 1; //D-perp
					tag_BC_mech[f][1] = 0; //N-perp
					tag_BC_mech[f][2] = 0; //D-parallel
					tag_BC_mech[f][3] = 1; //N-parallel
					tag_BC_mech[f][4] = nrm[0] * nUVW;// -alpha * P);
					tag_BC_mech[f][5] = nrm[1] * nUVW;// -alpha * P);
					tag_BC_mech[f][6] = nrm[2] * nUVW;// -alpha * P);
				}
			}
			else if (fabs(nrm[0] - 1) < 1.0e-6) // x-normal up
			{
				//fluid
				tag_BC_flow[f][0] = 1;
				tag_BC_flow[f][1] = 0;
				tag_BC_flow[f][2] = 0;
				//mech
				tag_BC_mech[f][0] = 0; //D-perp
				tag_BC_mech[f][1] = 1; //N-perp
				tag_BC_mech[f][2] = 0; //D-parallel
				tag_BC_mech[f][3] = 1; //N-parallel
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



void Test5::SetForce(Mesh& m, const INMOST::dynamic_variable& p, double T) const
{
	//TagVariableArray tag_F = m.GetTag("FORCE");
}


