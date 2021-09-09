#include "test8.h"

// [1] В.Е. Борисов, А.В. Иванов, Б.В. Критский, И.С.Меньшов, Е.Б.Савенков
// Численное моделирование задач пороупругости

using namespace INMOST;

const int istop = 1000000;

//shortcut for types
typedef Storage::real real;
typedef Storage::enumerator enumerator;
typedef Storage::var_array var_array;
typedef Storage::real_array real_array;
typedef Storage::reference_array ref_array;

const double pi = 3.1415926535897932384626433832795;
const double errtol = 1e-12;

double Test8::Gravity() const {return 0.0;}
double Test8::FluidDensity() const {return 1;}
double Test8::FluidViscosity() const {return mu;}
double Test8::SolidDensity() const {return 1;}
double Test8::Porosity(double _x, double _y, double _z) const {return 0.5;}
double Test8::InverseBiotModulus(double _x, double _y, double _z) const { return M ? 1.0/M : 0.0; }

Test8::Test8()
: 
	//a(1), b(1), E(1.0e+3), nu(0.25),
	//K(0.01), alpha(1.0), M(0.1), mu(1), F(10),
	a(10), b(1),  //a - column height, b - depth and width
	//E(2.129159824046921), nu(0.2998533724340176),
	E(100), nu(0.45),
	K(1.0e-11),
	//K(9.98e-10), 
	alpha(1.0), 
	//M(5446.623), 
	M(100000),
	mu(1.0e-3), 
	F(10),
	//derived parameters
	lambda(E* nu / (1 + nu) / (1 - 2 * nu)), 
	G(E / (1 + nu) / 2.0),
	Ku((lambda + 2 * G / 3.0 + alpha * alpha * M)),
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
	//std::cout << "(lambda+2*mu)*K:" << (lambda + 2 * mu) * K << std::endl;
	B0 = rMatrix::Unit(3)*alpha;
	K0 = rMatrix::Unit(3)*K;
	C0 = GenIsotropicTensor(E,nu);
}



vMatrix Test8::ElasticTensor(double _x, double _y, double _z) const {return C0;}
vMatrix Test8::BiotTensor(double _x, double _y, double _z) const { return B0; }
vMatrix Test8::PermTensor(double _x, double _y, double _z) const {return K0;}

hMatrix Test8::Displacement(double _x, double _y, double _z, double _t) const
{
	unknown x(_x,0), y(_y,1), z(_z,2), t(_t,3);
	hMatrix sol(3,1);
	//const double u0 = F * (a - _x) * (1 - 2 * nuu) / (2 * G *(1 - nuu));
	//const double cmg = (nuu - nu) / (2 * G * (1 - nu) * (1 - nuu));
	double du;
	double cu = -4.0 * F * a * (nuu - nu) / (G * pi * pi * (1 - nu) * (1 - nuu));
	sol.Zero();
	sol(0, 0) = F * (1 - 2*nu) * (a-_x) / (2 * G * (1-nu));
	//sol(0, 0) = u0 + cmg*F*(a-_x);
	for (int i = 0; i < istop; ++i)
	{
		du = cu * exp(-std::pow(1 + 2 * i, 2) * pi * pi * c * _t / 4.0 / a / a) * cos((1+2*i)*pi*_x/2.0/a)/ std::pow(1 + 2 * i, 2);
		//du = -8 * a / pi / pi * cmg * F * exp(-std::pow(1 + 2 * i, 2) * pi * pi * c * _t / 4.0 / a / a) * cos((1 + 2 * i) * pi * _x / 2.0 / a) / std::pow(1 + 2 * i, 2);
		sol(0, 0) += du;
		if (std::fabs(get_value(du)) < errtol * std::fabs(get_value(sol(0, 0))))
			break;
	}
//#pragma omp critical
//	std::cout << "z " << (a - _x) << " w " << get_value(sol(0, 0)) << " w/z " << get_value(sol(0, 0))/(a-_x) <<  " t " << _t << std::endl;
	if (sol.CheckNansInfs())
	{
#pragma omp critical
		{
			std::cout << "Bad sol" << std::endl;
			sol.Print();
			std::cout << "du: " << std::endl;
			std::cout << get_value(du) << std::endl;
		}
	}
	return sol;
}

static double _erfc(double x)
{
	const double a = 8 * (3 - pi) / (pi - 4.0) / (3.0 * pi);
	return 1.0 - x / sqrt(x * x + 1.0e-12) * sqrt(1.0 - exp(-x * x * (4.0 / pi + a * x * x) / (1.0 + a * x * x)));
}

hessian_variable Test8::Pressure(double _x, double _y, double _z, double _t) const
{
	unknown x(_x,0), y(_y,1), z(_z,2), t(_t,3);
	const double p0 = F * B * (1 + nuu) / 3.0 / (1 - nuu);
	double sol = p0, dp, ep, em;
	for (int i = 0; i < istop; ++i)
	{
		ep = std::erfc(((1.0 + 2.0 * i) * a + (_x - a)) / std::sqrt(4.0 * c * _t));
		em = std::erfc(((1.0 + 2.0 * i) * a - (_x - a)) / std::sqrt(4.0 * c * _t));
		dp = p0*std::pow(-1, i + 1) * (em + ep);
		sol += dp;
		if (std::fabs(get_value(dp)) < errtol * std::fabs(get_value(sol)))
			break;
	}
	if (check_nans_infs(sol))
	{
		std::cout << "Bad pressure " << sol << std::endl;
	}
	return sol;
}

rMatrix Test8::BCFlow_coef(double _x, double _y, double _z, double _t) const
{
	rMatrix ret(2,1);
	//impermiable
	ret(0,0) = 0;
	ret(1,0) = 1;
	return ret;
}

rMatrix Test8::BCMech_coef(double _x, double _y, double _z, double _t) const
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






void Test8::SetBC(Mesh & m, double T, MarkerType boundary) const
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
			if( fabs(nrm[0]+1) < 1.0e-6 ) // x-normal down
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
				tag_BC_mech[f][4] = -nrm[0] * F;
				tag_BC_mech[f][5] = -nrm[1] * F;
				tag_BC_mech[f][6] = -nrm[2] * F;
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




void Test8::Init(Mesh& m)
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
	double s[3] = { a,b,b };
	for (Mesh::iteratorNode it = m.BeginNode(); it != m.EndNode(); ++it)
	{
		real_array& c = it->Coords();
		for (int k = 0; k < c.size(); ++k)
			c[k] = (c[k] - cmin[k]) / (cmax[k] - cmin[k]) * s[k];
	}
	//m.RecomputeGeometricData();
}


void Test8::SetForce(Mesh& m, const INMOST::dynamic_variable& p, double T) const
{
	/*
	TagVariableArray tag_F = m.GetTag("FORCE");
#if defined(USE_OMP)
#pragma omp parallel for
#endif
	for (int i = 0; i < m.CellLastLocalID(); ++i) if (m.isValidCell(i))
		tag_F(m.CellByLocalID(i), 4, 1).Zero();
	*/
}


