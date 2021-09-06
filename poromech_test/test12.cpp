#include "test12.h"

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

double Test12::Gravity() const {return 0.0;}
double Test12::FluidDensity() const {return 998.2;}
double Test12::FluidViscosity() const {return mu;}
double Test12::SolidDensity() const {return 2000;}
double Test12::Porosity(double _x, double _y, double _z) const {return phi;}
double Test12::InverseBiotModulus(double _x, double _y, double _z) const { return M ? 1.0/M : 0.0; }

Test12::Test12()
: 
	a(1), b(1), c(1),
	k(1.0e-8),
	mu(0.01),
	Ks(1.0), 
	Kf(1.0),
	phi(0.5),
	F(10),
	forcebc(true)
{
	const double vC0[36] =
	{
		1,  2,  3,  4,  5,  6,
		2,  7,  8,  9, 10, 11,
		3,  8, 12, 13, 14, 15,
		4,  9, 13, 16, 17, 18,
		5, 10, 14, 17, 19, 20,
		6, 11, 15, 18, 20, 21
	};
	C = rMatrix(vC0, 6, 6);
	B(0, 0) = 1.0 - (C(0, 0) + C(0, 1) + C(0, 2)) / (3.0 * Ks);
	B(1, 1) = 1.0 - (C(0, 1) + C(1, 1) + C(1, 2)) / (3.0 * Ks);
	B(2, 2) = 1.0 - (C(0, 2) + C(1, 2) + C(2, 2)) / (3.0 * Ks);
	B(1, 2) = B(2, 1) = -(C(0, 3) + C(1, 3) + C(2, 3)) / (3.0 * Ks);
	B(0, 2) = B(2, 0) = -(C(0, 4) + C(1, 4) + C(2, 4)) / (3.0 * Ks);
	B(0, 1) = B(1, 0) = -(C(0, 5) + C(1, 5) + C(2, 5)) / (3.0 * Ks);
	M = Ks / (1.0 - (C(0, 0) + C(1, 1) + C(2, 2) + 2 * (C(0, 1) + C(0, 2) + C(1, 2))) / (9.0 * Ks) - phi * (1 - Ks / Kf));
	K = raMatrix::Unit(3) * k;

	rMatrix S = C.Invert();
	rMatrix alpha(6, 1);
	alpha(0, 0) = B(0, 0); // a11
	alpha(1, 0) = B(1, 1); // a22
	alpha(2, 0) = B(2, 2); // a33
	alpha(3, 0) = B(1, 2); // a23
	alpha(4, 0) = B(0, 2); // a13
	alpha(5, 0) = B(0, 1); // a12
	rMatrix beta = S * alpha;
	double C24 = S(1, 3) ? S(1, 3) : S(1, 1);
	double C34 = S(2, 3) ? S(2, 3) : S(1, 2);
	A11 = C24 * S(2, 2) - S(1, 2) * C34;
	A12 = C24 * beta(2, 2) - S(1, 2) * beta(1, 3);
	A13 = S(2, 2) * beta(1, 2) - C34 * beta(2, 2);
	A13 = C24 / A11;
	A15 = -S(2, 3) / A11;
	B11 = 1.0 / M + (alpha.DotProduct(beta) - alpha(1, 0) * beta(1, 0));
	B12 = alpha(0, 0) * S(0, 1) + alpha(2, 0) * S(1, 2) + alpha(3, 0) * S(1, 3) + alpha(4, 0) * S(1, 4) + alpha(5, 0) * S(1, 5);
	B13 = alpha(0, 0) * S(0, 2) + alpha(2, 0) * S(2, 2) + alpha(3, 0) * S(2, 3) + alpha(4, 0) * S(2, 4) + alpha(5, 0) * S(2, 5);
	p0 = (B12 * S(1, 2) - B13 * S(1, 1)) / (B12 * beta(1, 1) - B11 * S(1, 1)) * F / (2.0 * a);
	c1 = A11 / (A11 * B11 - A13 * B12 - A12 * B13) * k / mu;
	c2 = A11 * (A15 * B12 + A14 * B13) / (A11 * B11 - A13 * B12 - A12 * B13);
	c3 = (B13 * S(1, 2) - B12 * S(2, 2)) / (A11 * B11 - A13 * B12 - A12 * B13);
	Dmcoef = (B12 * S(1, 2) - B13 * S(1, 1)) / (B12 * beta(1, 0) - B11 * S(1, 1)) * F / a;
	lambda0 = 1 + A11 * A14 / (c2 * A12);

}



vMatrix Test12::ElasticTensor(double _x, double _y, double _z) const {return C0;}
vMatrix Test12::BiotTensor(double _x, double _y, double _z) const { return B0; }
vMatrix Test12::PermTensor(double _x, double _y, double _z) const {return K0;}

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

hMatrix Test12::Displacement(double _x, double _y, double _z, double _t) const
{
	unknown x(_x,0), y(_y,1), z(_z,2), t(_t,3);
	hMatrix sol(3,1);
	double du, dv;
	const double cnuu = (nuu - nu) / (1 - nu);
	double wi;
	sol.Zero();
	int stop = -1;
	sol(0, 0) = F / (G * a) * nu * _x / 2.0;
	sol(1, 0) = F / (G * a) * (nu - 1) * _y / 2.0;
	for (int i = 0; i < istop; ++i)
	{
		//cos(wi) = cnuu * sin(wi) / wi
		wi = solvetri(cnuu,0.5 * pi + pi * i);
		du = F / (G * a) * (a * sin(wi * _x / a) - nuu * _x * sin(wi)) * cos(wi) / (wi - sin(wi) * cos(wi)) * exp(-wi * wi * c * _t / a / a);
		dv = F / (G * a) * (1 - nuu) * _y * sin(wi) * cos(wi) / (wi - sin(wi) * cos(wi)) * exp(-wi * wi * c * _t / a / a);
		sol(0, 0) += du;
		sol(1, 0) += dv;
		if (std::fabs(get_value(du)) < epsuv * std::fabs(get_value(sol(0, 0))) && 
			std::fabs(get_value(dv)) < epsuv * std::fabs(get_value(sol(1, 0))))
			break;
	}
	return sol;
}

hessian_variable Test12::Pressure(double _x, double _y, double _z, double _t) const
{
	unknown x(_x,0), y(_y,1), z(_z,2), t(_t,3);
	double sol = 0, dp;
	const double cnuu = (nuu - nu) / (1 - nu);
	const double p0 = F * B * (1 + nuu) / (3.0 * a);
//#pragma omp critical
	if (_t)
	{
		for (int i = 0; i < istop; ++i)
		{
			//cos(wi) = cnuu * sin(wi) / wi
			//sin(wi) = wi*cos(wi)/cnuu
			double wi = solvetri(cnuu, 0.5 * pi + pi * i);
			dp = 2.0 * p0 * sin(wi) / (wi - sin(wi) * cos(wi)) * (cos(wi * _x / a) - cos(wi)) * exp(-wi * wi * c * _t / a / a);
			//dp = 2 * F * B * (1 + nuu) / (3.0 * a) * sin(wi) / (wi - sin(wi) * cos(wi)) * (cos(wi * _x / a) - cos(wi)) * exp(-wi * wi * c * _t / a / a);
			//dp = 2 * F * B * (1 + nuu) / (3.0 * a) * sin(wi) / (wi - sin(wi) * cos(wi)) * (cos(wi * _x / a) ) * exp(-wi * wi * c * _t / a / a);
			sol += dp;
			if (std::fabs(get_value(dp)) < epsp * std::fabs(get_value(sol)))
			{
//#pragma omp critical
//				std::cout << i << " p " << sol << " p0 " << F * B * (1 + nuu) / (3.0 * a) << std::endl;
				break;
			}
		}
	}
	else sol = p0;
	//else sol = F * B * (1 + nuu) / (6.0 * a);
	return sol;
}

rMatrix Test12::BCFlow_coef(double _x, double _y, double _z, double _t) const
{
	rMatrix ret(2,1);
	//impermiable
	ret(0,0) = 0;
	ret(1,0) = 1;
	return ret;
}

rMatrix Test12::BCMech_coef(double _x, double _y, double _z, double _t) const
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




void Test12::Init(Mesh& m)
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

void Test12::SetBC(Mesh & m, double T, MarkerType boundary) const
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



void Test12::SetForce(Mesh& m, const INMOST::dynamic_variable& p, double T) const
{
	//TagVariableArray tag_F = m.GetTag("FORCE");
}


