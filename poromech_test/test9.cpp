#include "test9.h"

// [1] Herminio T.Honorio, Clovis R.Maliska, Massimiliano Ferronato, Carlo Janna
// A stabilized element-based finite volume method for poroelastic problems

using namespace INMOST;

const int istop = 100000;

//shortcut for types
typedef Storage::real real;
typedef Storage::enumerator enumerator;
typedef Storage::var_array var_array;
typedef Storage::real_array real_array;
typedef Storage::reference_array ref_array;

const double pi = 3.1415926535897932384626433832795;
const double errtol = 1e-20;
const double epsnln = 1e-15;

double Test9::Gravity() const {return 0.0;}
double Test9::FluidDensity() const {return 1;}
double Test9::FluidViscosity() const {return mu;}
double Test9::SolidDensity() const {return 1;}
double Test9::Porosity(double _x, double _y, double _z) const {return 0.5;}
double Test9::InverseBiotModulus(double _x, double _y, double _z) const 
{ 
	if (_x < a1) //top
		return S1;
	else //top
		return S2;
}


static double trig_r(double c1, double c2, double x)
{
	return cos(x) + c1 * cos(c2 * x);
}

static double trig_j(double c1, double c2, double x)
{
	return -sin(x) - c1 * c2 * sin(c2 * x);
}


static double solvetri(double beta, double theta, double x, bool print)
{
	const double ratio = (std::sqrt(5.0) - 1.0) / 2.0;
	const double c1 = (1 - beta) / (1 + beta);
	const double c2 = (1 - theta) / (1 + theta);
	// j dx = -r, dx = -r/j, x += dx, x -= r/j
	{
		double r, j, nx, lsx;
		int nit = 0, lsit;
		do
		{
			r = trig_r(c1,c2,x);
			j = trig_j(c1,c2,x);
			nx = x - r / j;
			//line search
			{
				double pos[4] = { 0.0, 1.0 - ratio, ratio, 1.0 };
				double val[4] = { r, 0.0, 0.0, 0.0 };
				val[1] = trig_r(c1,c2,nx * pos[1] + x * (1.0 - pos[1]));
				val[2] = trig_r(c1,c2,nx * pos[2] + x * (1.0 - pos[2]));
				val[3] = trig_r(c1,c2,nx);
				lsit = 0;
				while (std::fabs(val[1] - val[2]) > epsnln * std::fabs(nx))
				{
					if (print) 
						std::cout << "pos " << pos[0] << " " << pos[1] << " " << pos[2] << " " << pos[3] << " val " << val[0] << " " << val[1] << " " << val[2] << " " << val[3] << std::endl;
					if (std::fabs(val[1]) < std::fabs(val[2]))
					{
						pos[3] = pos[2]; val[3] = val[2];
						pos[2] = pos[1]; val[2] = val[1];
						pos[1] = pos[3] - ratio * (pos[3] - pos[0]);
						val[1] = trig_r(c1,c2,nx * pos[1] + x * (1.0 - pos[1]));
					}
					else
					{
						pos[0] = pos[1]; val[0] = val[1];
						pos[1] = pos[2]; val[1] = val[2];
						pos[2] = pos[0] + ratio * (pos[3] - pos[0]);
						val[2] = trig_r(c1,c2,nx * pos[2] + x * (1.0 - pos[2]));
					}
					lsit++;
				}
				int kmin = 0;
				for (int k = 1; k < 4; ++k)
					if (std::fabs(val[k]) < std::fabs(val[kmin])) kmin = k;
				lsx = nx * pos[kmin] + x * (1.0 - pos[kmin]);
			}
			if( print )
				std::cout << nit << " x " << x << " r " << r << " j " << j << " new x " << nx << " new lsx " << lsx << " lsit " << lsit << " c1 " << c1 << " c2 " << c2 << std::endl;
			x = lsx;
			nit++;
		} while (std::fabs(r) > epsnln * std::fabs(x));
	}
	return x;
}

Test9::Test9()
	:
	//a(1), b(1), E(1.0e+3), nu(0.25),
	//K(0.01), alpha(1.0), M(0.1), mu(1), F(10),
	a1(7.5), a2(2.5), b(1),
	//E(2.129159824046921), nu(0.2998533724340176),
	E1(10), E2(2),
	nu1(0.45), nu2(0.15),
	//k1(1.0e-8), k2(1.0e-5),
	k1(1.0e-5), k2(1.0e-8),
	alpha1(1.0), alpha2(0.8),
	M1(1000),
	//M2(5000),
	mu(1.0e-2),
	F(10),
	//derived parameters
	S1(M1 ? 1.0 / M1 : 0.0),
	lambda1(E1* nu1 / (1 + nu1) / (1 - 2 * nu1)),
	G1(E1 / (1 + nu1) / 2.0),
	lambda2(E2* nu2 / (1 + nu2) / (1 - 2 * nu2)),
	G2(E2 / (1 + nu2) / 2.0),
	m1(1.0/(lambda1 + 2*G1)),
	m2(1.0/(lambda2 + 2*G2)),
	//M2(alpha1/alpha2*m1/m2*M1),//alpha2 m2 M2 = alpha1 m1 M1
	M2(alpha1*m1*M1 / (alpha2*m2*(1.0 + alpha1*m1*(alpha1 - alpha2)*M1))),
	S2(M2 ? 1.0 / M2 : 0.0),
	//alpha2(M1 / M2 * m1 / m2), //a2 m2 S1 = a1 m1 S2
	//Ku1((lambda1 + 2 * G1 / 3.0 + alpha1 * alpha1 * M1)),
	//Ku2((lambda2 + 2 * G2 / 3.0 + alpha2 * alpha2 * M2)),
	//B1(alpha1* M1 / Ku1),
	//B2(alpha2* M2 / Ku2),
	B1(alpha1 * m1 / (S1 + alpha1 * alpha1 * m1)),
	B2(alpha2 * m2 / (S2 + alpha2 * alpha2 * m2)),
	//nuu1((3 * nu1 + alpha1 * B1 * (1 - 2 * nu1)) / (3 - alpha1 * B1 * (1 - 2 * nu1))),
	//nuu2((3 * nu2 + alpha2 * B2 * (1 - 2 * nu2)) / (3 - alpha2 * B2 * (1 - 2 * nu2))),
	//c1(2 * k1 * G1 * (1 - nu1) * (nuu1 - nu1) / (mu * alpha1 * alpha1 * (1 - nuu1) * (1 - 2 * nu1) * (1 - 2 * nu1))),
	//c2(2 * k2 * G2 * (1 - nu2) * (nuu2 - nu2) / (mu * alpha2 * alpha2 * (1 - nuu2) * (1 - 2 * nu2) * (1 - 2 * nu2))),
	c1(k1/(mu* (S1 + alpha1 * alpha1 * m1))),
	c2(k2/(mu* (S2 + alpha2 * alpha2 * m2))),
	theta(a1 / a2 * std::sqrt(c2 / c1)),
	beta(k2 / k1 * std::sqrt(c1 / c2))
{
	std::cout << "Top layer: " << std::endl;
	std::cout << "E: " << E1 << " nu: " << nu1 << std::endl;
	std::cout << "K: " << k1 << " M " << M1 << std::endl;
	std::cout << "alpha: " << alpha1 << std::endl;
	std::cout << "lambda " << lambda1 << " G " << G1 << std::endl;
	std::cout << "consolidation coefficient: " << c1 << std::endl;
	//std::cout << "undrained Poisson: " << nuu1 << std::endl;
	std::cout << "Skempton's coefficient: " << B1 << std::endl;
	std::cout << "Specific storativity coefficient: " << S1 << std::endl;
	std::cout << "Confined compressibility: " << m1 << std::endl;

	std::cout << "Bottom layer: " << std::endl;
	std::cout << "E: " << E2 << " nu: " << nu2 << std::endl;
	std::cout << "K: " << k2 << " M " << M2 << std::endl;
	std::cout << "alpha: " << alpha2 << std::endl;
	std::cout << "lambda " << lambda2 << " G " << G2 << std::endl;
	std::cout << "consolidation coefficient: " << c2 << std::endl;
	//std::cout << "undrained Poisson: " << nuu2 << std::endl;
	std::cout << "Skempton's coefficient: " << B2 << std::endl;
	std::cout << "Specific storativity coefficient: " << S2 << std::endl;
	std::cout << "Confined compressibility: " << m2 << std::endl;

	std::cout << "theta " << theta << " beta " << beta << std::endl;
	//std::cout << "(lambda+2*mu)*K:" << (lambda + 2 * mu) * K << std::endl;
	bB1 = rMatrix::Unit(3) * alpha1;
	K1 = rMatrix::Unit(3) * k1;
	C1 = GenIsotropicTensor(E1, nu1);
	bB2 = rMatrix::Unit(3) * alpha2;
	K2 = rMatrix::Unit(3) * k2;
	C2 = GenIsotropicTensor(E2, nu2);

	std::cout << "beta: " << beta << std::endl;
	std::cout << "theta: " << theta << std::endl;
	std::cout << "Searching roots of trigonometric equation:" << std::endl;
	std::cout << (1-beta) / (1+beta) << "*cos(" << (1-theta)/(1+theta) << "*x)+cos(x) = 0, wi = x/" << 1+theta << std::endl;
	wi.resize(istop);
#pragma omp parallel for
	for (int i = 0; i < istop; ++i)
	{
		wi[i] = solvetri(beta, theta, pi * 0.5 + pi * i, false) / (1 + theta);
		if (std::fabs(wi[i] * (1 + theta) - (pi * 0.5 + pi * i)) > pi)
		{
#pragma omp critical
			{
				std::cout << "i " << i << " x " << wi[i] * (1 + theta) << " guess " << pi * 0.5 + pi * i << " farther than period!" << std::endl;
				solvetri(beta, theta, pi * 0.5 + pi * i, true);
			}
		}
	}
}



vMatrix Test9::ElasticTensor(double _x, double _y, double _z) const 
{
	if (_x < a1) //top layer
		return C1;
	else //bottom layer
		return C2;
}
vMatrix Test9::BiotTensor(double _x, double _y, double _z) const 
{ 
	if (_x < a1) //top layer
		return bB1;
	else //bottom layer
		return bB2;
}
vMatrix Test9::PermTensor(double _x, double _y, double _z) const 
{
	if (_x < a1) //top layer
		return K1;
	else //bottom layer
		return K2;
}

double Test9::Qfunc(double wi, double t) const
{
	const double p0 = B1 * F;
	return 2.0 * p0 * exp(-wi * wi * c2 * t / a2 / a2) / wi / ((1 + beta * theta) * cos(theta * wi) * sin(wi) + (beta + theta) * sin(theta * wi) * cos(wi));
}

double Test9::P1func(double wi, double x, double t) const
{
	return Qfunc(wi,t) * (cos(wi) * cos(theta * wi * (a1 - x) / a1) - beta * sin(wi) * sin(theta * wi * (a1 - x) / a1));
}

double Test9::intP1func(double wi, double x, double t) const
{
	return Qfunc(wi,t) / wi / theta * a1 * (cos(wi) * sin(theta * wi * (a1 - x) / a1) + beta * sin(wi) * cos(theta * wi * (a1 - x) / a1));
}

double Test9::P2func(double wi, double x, double t) const
{
	return Qfunc(wi,t) * (cos(wi) * cos(wi * (a1 - x) / a2) - sin(wi) * sin(wi * (a1 - x) / a2));
}

double Test9::intP2func(double wi, double x, double t) const
{
	return Qfunc(wi,t) / wi * a2 * (cos(wi) * sin(wi * (a1 - x) / a2) + sin(wi) * cos(wi * (a1 - x) / a2));
}



hMatrix Test9::Displacement(double _x, double _y, double _z, double _t) const
{
	unknown x(_x,0), y(_y,1), z(_z,2), t(_t,3);
	hMatrix sol(3,1);
	double du;
	sol.Zero();
	
	if (_t)
	{
		
		if (_x < a1)
		{
			sol(0, 0) = (a1 - _x) * m1 * F + a2 * m2 * F;
			for (int i = 0; i < istop; ++i)
			{
				du = 0.0;
				//du += a2 * (alpha2 - alpha1) * m2 * P2func(wi[i], a1, _t);
				du -= alpha2 * m2 * (intP2func(wi[i], a1, _t) - intP2func(wi[i], a1 + a2, _t));
				du -= alpha1 * m1 * (intP1func(wi[i], _x, _t) - intP1func(wi[i], a1, _t));
				sol(0, 0) += du;
				if (std::fabs(get_value(du)) < errtol * std::fabs(get_value(sol(0, 0))))
					break;
			}
		}
		else
		{
			sol(0, 0) = (a1 + a2 - _x) * m2 * F;
			for (int i = 0; i < istop; ++i)
			{
				du = 0.0;
				//du += (a1 + a2 - _x) * (alpha2 - alpha1) * m2 * P2func(wi[i], a1, _t);
				du -= alpha2 * m2 * (intP2func(wi[i], _x, _t) - intP2func(wi[i], a1 + a2, _t));
				sol(0, 0) += du;
				if (std::fabs(get_value(du)) < errtol * std::fabs(get_value(sol(0, 0))))
					break;
			}
		}
	}
	else
	{
		const double p0 = B1 * F;
		if (_x < a1)
			sol(0, 0) = (a1 - _x) * m1 * (F - alpha1*p0) + a2 * m2 * (F - alpha2*p0);
		else
			sol(0, 0) = (a1 + a2 - _x) * m2 * (F - alpha2*p0);
	}
	
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


hessian_variable Test9::Pressure(double _x, double _y, double _z, double _t) const
{
	unknown x(_x,0), y(_y,1), z(_z,2), t(_t,3);
	const double p0 = B1*F;
	double sol = p0, dp;
	//int its = 0;
	if (_t)
	{
		sol = 0.0;
		const double zz = a1 - _x;
		for (int i = 0; i < istop; ++i)
		{
			if (_x < a1) //top layer
				dp = P1func(wi[i],_x,_t);
			else //bottom layer
				dp = P2func(wi[i],_x,_t);
			sol += dp;
			//its = i;
			if (std::fabs(get_value(dp)) < errtol * std::fabs(get_value(sol)))
				break;
		}
//#pragma omp critical
//		std::cout << "z " << zz << " p " << sol << " p0 " << p0 << " its " << its << " dp " << dp << " Q " << Q << " t " << _t << std::endl;
	}

	if (check_nans_infs(sol))
	{
		std::cout << "Bad pressure " << sol << std::endl;
	}
	return sol;
}

rMatrix Test9::BCFlow_coef(double _x, double _y, double _z, double _t) const
{
	rMatrix ret(2,1);
	//impermiable
	ret(0,0) = 0;
	ret(1,0) = 1;
	return ret;
}

rMatrix Test9::BCMech_coef(double _x, double _y, double _z, double _t) const
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






void Test9::SetBC(Mesh & m, double T, MarkerType boundary) const
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

void Test9::Init(Mesh& m)
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
	double s[3] = { a1+a2,b,b };
	for (Mesh::iteratorNode it = m.BeginNode(); it != m.EndNode(); ++it)
	{
		real_array& c = it->Coords();
		for (int k = 0; k < c.size(); ++k)
			c[k] = (c[k] - cmin[k]) / (cmax[k] - cmin[k]) * s[k];
	}
	//m.RecomputeGeometricData();
}




void Test9::SetForce(Mesh& m, const INMOST::dynamic_variable& p, double T) const
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


