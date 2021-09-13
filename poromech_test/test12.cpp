#include "test12.h"

// [1] Y. ABOUSLEIMAN, A. H.-D. CHENG, L. CUI, E. DETOURNAY and J.-C. ROEGIERS
// Mandel's problem revisited

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

static double trig_r(double coef, double x)
{
	return cos(x) - coef * sin(x) / x;
}

static double trig_j(double coef, double x)
{
	return -((x*x - coef) * sin(x) + coef * x * cos(x)) / (x*x);
}


static double solvetri(double coef, double x, bool print)
{
	const double ratio = (std::sqrt(5.0) - 1.0) / 2.0;
	// j dx = -r, dx = -r/j, x += dx, x -= r/j
	{
		double r, j, nx, lsx;
		int nit = 0, lsit;
		do
		{
			r = trig_r(coef,x);
			j = trig_j(coef,x);
			nx = x - r / j;
			//line search
			{
				double pos[4] = { 0.0, 1.0 - ratio, ratio, 1.0 };
				double val[4] = { r, 0.0, 0.0, 0.0 };
				val[1] = trig_r(coef, nx * pos[1] + x * (1.0 - pos[1]));
				val[2] = trig_r(coef, nx * pos[2] + x * (1.0 - pos[2]));
				val[3] = trig_r(coef, nx);
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
						val[1] = trig_r(coef, nx * pos[1] + x * (1.0 - pos[1]));
					}
					else
					{
						pos[0] = pos[1]; val[0] = val[1];
						pos[1] = pos[2]; val[1] = val[2];
						pos[2] = pos[0] + ratio * (pos[3] - pos[0]);
						val[2] = trig_r(coef, nx * pos[2] + x * (1.0 - pos[2]));
					}
					lsit++;
				}
				int kmin = 0;
				for (int k = 1; k < 4; ++k)
					if (std::fabs(val[k]) < std::fabs(val[kmin])) kmin = k;
				lsx = nx * pos[kmin] + x * (1.0 - pos[kmin]);
			}
			if (print)
				std::cout << nit << " x " << x << " r " << r << " j " << j << " new x " << nx << " new lsx " << lsx << " lsit " << lsit << " coef " << coef << std::endl;
			x = lsx;
			nit++;
		} while (std::fabs(r) > epsnln * std::fabs(x));
	}
	return x;
}

Test12::Test12()
	:
	a(1), b(2),
	//transverse isotropic material
	//Ex(10), Ey(2), Ez(Ex),
	//nuxy(0.25), nuxz(0.3), nuyz(Ex/Ey*nuxy),
	//Gxy(2), Gxz(0.5*Ex/(1+nuxz)), Gyz(Gxy),
	Ex(2.5), Ey(Ex), Ez(1.5),
	nuxy(0.15), nuxz(0.35), nuyz(nuxz),
	Gxy(0.5 * Ex / (1 + nuxy)), Gxz(2), Gyz(Gxz),
	k1(1.0e-6), k2(k1), k3(1.0e-8),
	mu(0.01),
	Ks(3.0/((1.0-2*nuxy-2*nuxz)/Ex + (1.0-2*nuyz)/Ey + 1.0/Ez)), 
	//Ks(1.0/(1.0/Ex+1.0/Ey+1.0/Ez)),
	Kf(1.5),
	phi(0.5),
	nuzx(Ez/Ex*nuxz), nuyx(Ey/Ex*nuxy),
	M11(Ex* (Ez - Ex * nuzx * nuzx) / (1 + nuyx) / (Ez - Ez * nuyx - 2 * Ex * nuzx * nuzx)),
	M12(Ex* (nuyx* Ez + Ex * nuzx * nuzx) / (1 + nuyx) / (Ez - Ez * nuyx - 2 * Ex * nuzx * nuzx)),
	M13(Ex*Ez*nuzx/ (Ez - Ez * nuyx - 2 * Ex * nuzx * nuzx)),
	M33(Ez*Ez*(1-nuyx)/ (Ez - Ez * nuyx - 2 * Ex * nuzx * nuzx)),
	M55(Gxz),
	alpha1(1.0-(M11+M12+M13)/(3.0*Ks)), alpha2(alpha1), alpha3(1-(2*M13+M33)/(3.0*Ks)),
	M(Ks/(1 + phi*(Ks/Kf-1.0)-(2*M11+M33+2*M12+4*M13)/(9.0*Ks))),
	c1(k1*M*M11/mu/(M11+alpha1*alpha1*M)),
	//c3(k3*M*M33/mu/(M33+alpha3*alpha3*M)),
	A1((alpha1*alpha1*M33-2*alpha1*alpha3*M13+alpha3*alpha3*M11)/(alpha3*M11-alpha1*M13)+(M11*M33-M13*M13)/(M*(alpha3*M11-alpha1*M13))),
	A2((alpha3*M11-alpha1*M13)/M11),
	F(10),
	forcebc(false)
{
	const double vK[3] = { k1,k2,k3 }, vB[3] = { alpha1,alpha2,alpha3 };
	C = GenAnisotropicTensor(Ex, Ey, Ez, nuxy, nuxz, nuyz, Gxy, Gxz, Gyz);
	K = rMatrix::FromDiagonal(vK, 3);
	B = rMatrix::FromDiagonal(vB, 3);
	wi.resize(istop);
#pragma omp parallel for
	for (int i = 0; i < istop; ++i)
	{
		wi[i] = solvetri(A2/A1, pi * 0.5 + pi * i, false);
		if (std::fabs(wi[i] - (pi * 0.5 + pi * i)) > pi)
		{
#pragma omp critical
			{
				std::cout << "i " << i << " x " << wi[i] << " guess " << pi * 0.5 + pi * i << " farther than period!" << std::endl;
				solvetri(A2/A1, pi * 0.5 + pi * i, true);
			}
		}
	}
	std::cout << "C:" << std::endl;
	C.Print();
	std::cout << "K:" << std::endl;
	K.Print();
	std::cout << "B:" << std::endl;
	B.Print();
	std::cout << "M: " << M << std::endl;
	std::cout << "nuzx: " << nuzx << " nuyx: " << nuyx << std::endl;
	std::cout << "Ks: " << Ks << " Kf " << Kf << std::endl;
	std::cout << "M11: " << M11 << " M12 " << M12 << " M13 " << M13 << " M33 " << M33 << " M55 " << M55 << " c1 " << c1 << " A1 " << A1 << " A2 " << A2 << std::endl;
	std::cout << "C11: " << C(0,0) << " C12 " << C(0,1) << " C13 " << C(0,2) << " C33 " << C(2,2) << " C55 " << C(4,4) << std::endl;
}



vMatrix Test12::ElasticTensor(double _x, double _y, double _z) const {return C;}
vMatrix Test12::BiotTensor(double _x, double _y, double _z) const { return B; }
vMatrix Test12::PermTensor(double _x, double _y, double _z) const {return K;}

hMatrix Test12::Displacement(double _x, double _y, double _z, double _t) const
{
	unknown x(_x,0), y(_y,1), z(_z,2), t(_t,3);
	hMatrix sol(3,1);
	double du1, du2, dw;
	const double uc1 = F / a * M13 / (M11 * M33 - M13 * M13);
	const double uc2 = -2 * F / a * (alpha1 * alpha3 * M + M13) / (A1*M*(alpha3*M11-alpha1*M13));
	const double uc3 = 2 * F * alpha1 / A2 / M11;
	const double wc1 = -F / a * M11 / (M11 * M33 - M13 * M13);
	const double u0 = (uc1 + uc2*A2/(2.0*(A1-A2)) + uc3*A2/(2*(A1-A2))) * _x;
	const double w0 = wc1 * (1.0 - A2 / A1) * _z;
	double cu0 = uc1 * _x, cu1 = 0, cu2 = 0, cw0 = wc1 * _z, cw1 = 0, S1, S2;
	int stop = 0;
	sol.Zero();
	if (_t)
	{
		sol(0, 0) = uc1 * _x;
		sol(2, 0) = wc1 * _z;
		
		for (int i = 0; i < istop; ++i)
		{
			S1 = cos(wi[i]) * sin(wi[i]) / (wi[i] - sin(wi[i]) * cos(wi[i])) * exp(-wi[i] * wi[i] * c1 * _t / a / a);
			S2 = cos(wi[i]) * sin(wi[i] * _x / a) / (wi[i] - sin(wi[i]) * cos(wi[i])) * exp(-wi[i] * wi[i] * c1 * _t / a / a);
			du1 = uc2 * S1 * _x;
			du2 = uc3 * S2;
			dw = wc1 * 2.0 * (A2 / A1 - 1.0) * S1 * _z;
			sol(0, 0) += du1 + du2;
			sol(2, 0) += dw;
			cu1 += du1;
			cu2 += du2;
			cw1 += dw;
			if (std::fabs(get_value(du1 + du2)) < epsuv * std::fabs(get_value(sol(0, 0))) &&
				std::fabs(get_value(dw)) < epsuv * std::fabs(get_value(sol(2, 0))))
			{
				stop = i;
				break;
			}
		}
	}
	else
	{
		sol(0, 0) = u0;
		sol(2, 0) = w0;
	}
//#pragma omp critical
//	std::cout << "u: " << get_value(sol(0, 0)) << "; u0: " << u0 << "; w: " << get_value(sol(2, 0)) << "; w0: " << w0 << "; iters: " << stop << "; x: " << _x << "; z: " << _z << "; t: " << _t << std::endl;
	return sol;
}

hessian_variable Test12::Pressure(double _x, double _y, double _z, double _t) const
{
	unknown x(_x,0), y(_y,1), z(_z,2), t(_t,3);
	double sol = 0, dp, Di;
	const double p0 = F / (a * A1);
	if (_t)
	{
		for (int i = 0; i < istop; ++i)
		{
			Di = 2 * F * (A2 - A1) / a / A1 * sin(wi[i]) / (wi[i] - sin(wi[i]) * cos(wi[i])); //formula (39) in [1]
			dp = Di/(A2-A1)*(cos(wi[i]*_x/a)-cos(wi[i]))*exp(-wi[i]*wi[i]*c1*_t/a/a); //formula (34) in [1]
			sol += dp;
			if (std::fabs(get_value(dp)) < epsp * std::fabs(get_value(sol)))
				break;
		}
	}
	else sol = p0;
//#pragma omp critical
//	std::cout << "p " << sol << " p0 " << p0 << " x " << _x << " t " << _t << std::endl;
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
	double s[3] = { a, a, b };
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
			if( fabs(nrm[2]-1) < 1.0e-6 ) // z-normal up
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


