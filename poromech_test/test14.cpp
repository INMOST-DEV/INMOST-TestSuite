#include "test14.h"

using namespace INMOST;

//shortcut for types
typedef Storage::real real;
typedef Storage::enumerator enumerator;
typedef Storage::var_array var_array;
typedef Storage::real_array real_array;
typedef Storage::reference_array ref_array;

double Test14::Gravity() const { return 0.0; }
double Test14::FluidDensity() const {return 20;}
double Test14::FluidViscosity() const {return 0.1;}
double Test14::SolidDensity() const {return 120;}
double Test14::Porosity(double _x, double _y, double _z) const
{
	throw -1;
	return 0.0;
}
double Test14::InverseBiotModulus(double _x, double _y, double _z) const
{
	throw -1;
	return 0.0;
}



Test14::Test14() : pbhp(400), WI(500000), trans_scale(1.0) {}

vMatrix Test14::ElasticTensor(double _x, double _y, double _z) const
{
	throw -1;
	return rMatrix::Unit(6);
}

vMatrix Test14::BiotTensor(double _x, double _y, double _z) const
{
	throw -1;
	return rMatrix::Unit(3);
}

vMatrix Test14::PermTensor(double _x, double _y, double _z) const
{
	throw -1;
	return rMatrix::Unit(3);
}

hMatrix Test14::Displacement(double _x, double _y, double _z, double _t) const
{
	hMatrix sol(3,1);
	sol(0,0) = 0; //u
	sol(1,0) = 0; //v
	sol(2,0) = 0; //w
	return sol;
}

hessian_variable Test14::Pressure(double _x, double _y, double _z, double _t) const
{
	return 250;
}

rMatrix Test14::BCFlow_coef(double _x, double _y, double _z, double _t) const
{
	throw -1;
	rMatrix ret(2,1);
	ret(0,0) = 0; //dirichlet parallel
	ret(1,0) = 1; //neumann parallel
	return ret;
}

rMatrix Test14::BCMech_coef(double _x, double _y, double _z, double _t) const
{
	throw -1;
	rMatrix ret(4,1);
	ret(0,0) = 1; //dirichlet perp
	ret(1,0) = 0; //neumann perp
	ret(2,0) = 0; //dirichlet parallel
	ret(3,0) = 1; //neumann parallel
	return ret;
}

void Test14::Init(Mesh & m)
{
	//prepare wells
	int nx = 6, ny = 22, tot = nx*ny;
	double wx = 240, wy = 440, wz = 340;
	double hx = 240 / 60, hy = 440 / 220;
	double cntc[3], cntw[3];
	wellcells.SetMeshLink(&m);
	SearchKDTree t(&m);
	TagReal tag_WI = m.CreateTag("WI", DATA_REAL, CELL, CELL, 1);
	for (int j = 0; j < ny; ++j)
		for (int i = 0; i < nx; ++i)
		{
			cntw[0] = 0.5*hx + i * (wx - hx) / (nx - 1);
			cntw[1] = 0.5*hy + j * (wy - hy) / (ny - 1);
			cntw[2] = 0.5 * wz;
			//std::cout << "i " << i << " j " << j << " ind " << (i + j * nx) << " center " << cntw[0] << " " << cntw[1] << " " << cntw[2] << std::endl;
			Cell c = t.SearchCell(cntw,false);
			if (c.isValid() && c.GetStatus() != Element::Ghost)
			{
				c->Centroid(cntc);
				std::cout << "proc " << m.GetProcessorRank() << " found w" << (i + j * nx) << " cell " << c.LocalID() << " center " << cntc[0] << " " << cntc[1] << " " << cntc[2] << std::endl;
				tag_WI[c] = WI;
				wellcells.push_back(c);
			}
		}
	
	if(m.HaveTag("PERM"))
	{
		TagRealArray perm = m.GetTag("PERM");
		int perm_size = perm.GetSize();
		for (Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it)
			for (int k = 0; k < perm_size; ++k) perm[*it][k] *= trans_scale;
	}
	
}

void Test14::SetBC(Mesh & m, double T, MarkerType boundary) const
{
	TagRealArray tag_BC_flow = m.GetTag("BOUNDARY_CONDITION_FLOW");
	TagRealArray tag_BC_mech = m.GetTag("BOUNDARY_CONDITION_ELASTIC");
	real BCf[3] = {0,1,0};
	real BCm[7] = { 1,0,0,1,0,0,0 };
	rMatrix mBCf(BCf, 3,1);
	rMatrix mBCm(BCm, 7,1);
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
			tag_BC_flow(f,3,1) = mBCf;
			tag_BC_mech(f,7,1) = mBCm;
		}
	}
}

void Test14::SetForce(Mesh & m, const INMOST::dynamic_variable & p, double T) const
{
	TagVariableArray tag_F = m.GetTag("FORCE");
#if defined(USE_OMP)
#pragma omp parallel for
#endif
	for (int i = 0; i < m.CellLastLocalID(); ++i) if (m.isValidCell(i))
	{
		Cell c = m.CellByLocalID(i);
		if (c.GetStatus() == Element::Ghost) continue;
		if (c.HaveData(tag_F)) c.DelData(tag_F);
	}
	double mu = FluidViscosity();
	for(ElementArray<Cell>::const_iterator it = wellcells.begin(); it != wellcells.end(); ++it)
		tag_F[*it][3] -= WI/mu*(p[it->self()] - pbhp)/it->Volume();
}

void Test14::SetProperty(Mesh & m) const
{
	const bool avgprop = false;
	//supposed to be defined on the mesh
	TagReal poro = m.GetTag("PORO");
	TagRealArray perm = m.GetTag("PERM");
	if( !poro.isValid() )
	{
		std::cout << "Poro was not defined" << std::endl;
		exit(-1);
	}
	if( !perm.isValid() )
	{
		std::cout << "Perm was not defined" << std::endl;
		exit(-1);
	}
	//this should be calculated
	
	TagRealArray tag_E = m.GetTag("ELASTIC_TENSOR");
	TagRealArray tag_B = m.GetTag("BIOT_COEFFICIENT");
	TagReal tag_M    = m.GetTag("INVERSE_BIOT_MODULUS");
	double S0min = 1.0e+20, S0max = -1.0e+20;
	double S1min = 1.0e+20, S1max = -1.0e+20;
	double S2min = 1.0e+20, S2max = -1.0e+20;
	rMatrix K(3, 3);
	for (int it = 0; it < m.CellLastLocalID(); ++it) if (m.isValidCell(it))
	{
		Cell c = m.CellByLocalID(it);
		K = rMatrix::FromTensor(perm[c].data(), perm[c].size(), 3);
		double S0 = K(0,0);
		double S1 = K(1,1);
		double S2 = K(2,2);
		S0min = std::min(S0, S0min);
		S0max = std::max(S0, S0max);
		S1min = std::min(S1, S1min);
		S1max = std::max(S1, S1max);
		S2min = std::min(S2, S2min);
		S2max = std::max(S2, S2max);
	}
#if defined(USE_OMP)
#pragma omp parallel
#endif
	{
		//rMatrix C(6, 6, 0.0);
		rMatrix K(3, 3);
#if defined(USE_OMP)
#pragma omp for
#endif
		for(int it = 0; it < m.CellLastLocalID(); ++it) if( m.isValidCell(it) )
		{
			Cell c = m.CellByLocalID(it);
			K = rMatrix::FromTensor(perm[c].data(), perm[c].size(), 3);
			double S0 = K(0,0);
			double S1 = K(1,1);
			double S2 = K(2,2);
			double avg_poro = poro[c];


			if (avgprop)
			{
				ElementArray<Cell> around;// = c->BridgeAdjacencies2Cell(NODE);
				double tot_vol = c.Volume();
				S0 *= tot_vol;
				S1 *= tot_vol;
				S2 *= tot_vol;
				avg_poro *= tot_vol;
				//tag_B(c,3,3) = perm(c,3,3)*c.Volume();
				if (!around.empty())
					for (ElementArray<Cell>::iterator kt = around.begin(); kt != around.end(); ++kt)
					{
						K = rMatrix::FromTensor(perm[*kt].data(), perm[*kt].size(), 3);
						//tag_B(c,3,3) += perm(kt->self(),3,3)*kt->Volume();
						avg_poro += poro[kt->self()] * kt->Volume();
						S0 += K(0,0) * kt->Volume();
						S1 += K(1,1) * kt->Volume();
						S2 += K(2,2) * kt->Volume();
						tot_vol += kt->Volume();
					}
				//tag_B(c,3,3) /= tot_vol;
				avg_poro /= tot_vol;
				S0 /= tot_vol;
				S1 /= tot_vol;
				S2 /= tot_vol;
			}

			
			
			tag_B(c,3,3) = rMatrix::Unit(3)*((avg_poro+1)*0.5);
			tag_B(c,3,3).CheckNans();
			
			real E1,E2,E3,nu12,nu13,nu23,G23,G13,G12;

			
			/*
			E1 = S0 * 1000;
			E2 = S1 * 1000;
			E3 = S2 * 1000;
			G23 = 500 * S0;
			G13 = 500 * S1;
			G12 = 500 * S2;
			nu12 = 0.1;// / E2 * E1;
			nu13 = 0.1;// / E3 * E1;
			nu23 = 0.2;// / E3 * E2;
			*/
			/*
			E1 = E2 = E3 = 1000;
			nu12 = nu13 = nu23 = 0.35;
			G12 = G13 = G23 = 0.5 * E1 / (1 + nu12);
			*/
			/*
			E1 = (0.75 * (S0 - S0min) / (S0max - S0min) + 0.25) * 1000;
			E2 = (0.75 * (S1 - S1min) / (S1max - S1min) + 0.25) * 1000;
			E3 = (0.75 * (S2 - S2min) / (S2max - S2min) + 0.25) * 1000;
			G23 = 500 * (0.75 * (S0 - S0min) / (S0max - S0min) + 0.25);
			G13 = 500 * (0.75 * (S1 - S1min) / (S1max - S1min) + 0.25);
			G12 = 500 * (0.75 * (S2 - S2min) / (S2max - S2min) + 0.25);
			nu12 = 0.1;// / E2 * E1;
			nu13 = 0.1;// / E3 * E1;
			nu23 = 0.2;// / E3 * E2;
			*/

			E1 = (0.5 * (S0 - S0min) / (S0max - S0min) + 0.5) * 10000;
			E2 = (0.5 * (S1 - S1min) / (S1max - S1min) + 0.5) * 10000;
			E3 = (0.5 * (S2 - S2min) / (S2max - S2min) + 0.5) * 10000;
			G23 = 5000 * (0.5 * (S0 - S0min) / (S0max - S0min) + 0.5);
			G13 = 5000 * (0.5 * (S1 - S1min) / (S1max - S1min) + 0.5);
			G12 = 5000 * (0.5 * (S2 - S2min) / (S2max - S2min) + 0.5);
			nu12 = 0.1;// / E2 * E1;
			nu13 = 0.1;// / E3 * E1;
			nu23 = 0.2;// / E3 * E2;

			//C.Print();
			if (tag_E[c].size() == 21)
			{
				raSymmetricMatrix E = raSymmetricMatrixMake(tag_E[c].data(), 6);
				E = GenAnisotropicTensor(E1, E2, E3, nu12, nu13, nu23, G12, G13, G23);
				E.CheckNans();

				if (!rCheckEigen(E, false))
				{
#pragma omp critical
					{
						real nu21, nu31, nu32;
						nu21 = nu12 / E1 * E2;
						nu31 = nu13 / E1 * E3;
						nu32 = nu23 / E2 * E3;
						std::cout << "K:" << std::endl;
						perm(c, 3, 3).Print();
						std::cout << "S: " << S0 << " " << S1 << " " << S2 << std::endl;
						
						std::cout << "E1: " << E1 << "; E2: " << E2 << "; E3: " << E3 << ";" << std::endl;
						std::cout << "nu12: " << nu12 << "; nu13: " << nu13 << "; nu23: " << nu23 << std::endl;// << nu23 << std::endl;
						std::cout << "G12: " << G12 << "; G13: " << G13 << "; G23: " << G23 << "; " << std::endl;
						std::cout << "nu21: " << nu21 << " nu31: " << nu31 << " nu32: " << nu32 << std::endl;

						std::cout << "E:" << std::endl;
						E.Print();

						std::cout << "negative!" << std::endl;
						rCheckEigen(tag_E(c, 6, 6));
					}

				}
			}
			else if (tag_E[c].size() == 36)
			{
				tag_E(c, 6, 6) = GenAnisotropicTensor(E1, E2, E3, nu12, nu13, nu23, G12, G13, G23);
				tag_E(c, 6, 6).CheckNans();

				if (!rCheckEigen(tag_E(c, 6, 6), false))
				{
					std::cout << "K:" << std::endl;
					perm(c, 3, 3).Print();
					std::cout << "S: " << S0 << " " << S1 << " " << S2 << std::endl;
					


					std::cout << "E:" << std::endl;
					tag_E(c, 6, 6).Print();

					std::cout << "negative!" << std::endl;
					rCheckEigen(tag_E(c, 6, 6));

				}
			}
			else std::cout << __FILE__ << ":" << __LINE__ << " do not know what to do with elasticity tensor of size " << tag_E.GetSize() << std::endl;
			
			tag_M[c] = 0.5/(1.0 + avg_poro);
		}
	}
}











