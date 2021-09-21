#include "test14.h"

using namespace INMOST;

//shortcut for types
typedef Storage::real real;
typedef Storage::enumerator enumerator;
typedef Storage::var_array var_array;
typedef Storage::real_array real_array;
typedef Storage::reference_array ref_array;

double Test14::Gravity() const {return -9.8;}
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



Test14::Test14() : pbhp(400), WI(50000) {}

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
	double wx = 60, wy = 220, wz = 85;
	std::vector< double > cnt(tot*3);
	double cntc[3];
	for(int i = 0; i < nx; ++i)
		for (int j = 0; j < ny; ++j)
		{
			std::cout << "i " << i << " j " << j << " ind " << (i + j * nx) <<  " center " << 0.5 + i * (wx - 1.0) / (nx - 1) << " " << 0.5 + j * (wy - 1.0) / (ny - 1) << " " << 0.5 + 0.5 * wz << std::endl;
			cnt[(i + j * nx) * 3 + 0] = 0.5 + i * (wx - 1.0) / (nx - 1);
			cnt[(i + j * nx) * 3 + 1] = 0.5 + j * (wy - 1.0) / (ny - 1);
			cnt[(i + j * nx) * 3 + 2] = 0.5 + 0.5 * wz;
		}
	wellcells.SetMeshLink(&m);
	wellcells.resize(tot);
	TagReal tag_WI = m.CreateTag("WI",DATA_REAL,CELL,CELL,1);
	for(Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it) if( it->GetStatus() != Element::Ghost )
	{
		for(int q = 0; q < tot; ++q) if( !wellcells[q].isValid() && it->Inside(&cnt[q*3]) )
		{
			wellcells[q] = it->self();
			it->Centroid(cntc);
			std::cout << "proc " << m.GetProcessorRank() << " found w" << q << " center " << cnt[q * 3 + 0] << " " << cnt[q * 3 + 1] << " " << cnt[q * 3 + 2] << " cell " << wellcells[q].LocalID() << " center " << cntc[0] << " " << cntc[1] << " " << cntc[2] << std::endl;
			tag_WI[wellcells[q]] = WI;
		}
	}
	TagRealArray perm;
	if( m.HaveTag("PERM") )
		perm = m.GetTag("PERM");
	if( perm.isValid() && perm.GetSize() != 9 )
	{
		int perm_size = perm.GetSize();
		std::vector<real> perm_data;
		for(Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it)
			perm_data.insert(perm_data.end(),perm[it->self()].begin(),perm[it->self()].end());
		int k = 0;
		perm = m.DeleteTag(perm);
		perm = m.CreateTag("PERM",DATA_REAL,CELL,NONE,9);
		for(Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it)
		{
			perm(it->self(),3,3) = rMatrix::FromTensor(&perm_data[k*perm_size],perm_size,3);
			k++;
		}
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
#if defined(USE_OMP)
#pragma omp parallel
#endif
	{
		//rMatrix C(6, 6, 0.0);
#if defined(USE_OMP)
#pragma omp for
#endif
		for(int it = 0; it < m.CellLastLocalID(); ++it) if( m.isValidCell(it) )
		{
			real cnt[3];
			Cell c = m.CellByLocalID(it);
			c.Barycenter(cnt);
			
			double S0 = perm(c,3,3)(0,0);
			double S1 = perm(c,3,3)(1,1);
			double S2 = perm(c,3,3)(2,2);
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
						//tag_B(c,3,3) += perm(kt->self(),3,3)*kt->Volume();
						avg_poro += poro[kt->self()] * kt->Volume();
						S0 += perm(kt->self(), 3, 3)(0, 0) * kt->Volume();
						S1 += perm(kt->self(), 3, 3)(1, 1) * kt->Volume();
						S2 += perm(kt->self(), 3, 3)(2, 2) * kt->Volume();
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

			
			E1 = S0 * 1000;
			E2 = S1 * 1000;
			E3 = S2 * 1000;
			G23 = 500 * S0;
			G13 = 500 * S1;
			G12 = 500 * S2;
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











