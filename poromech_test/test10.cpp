#include "test10.h"

using namespace INMOST;

//shortcut for types
typedef Storage::real real;
typedef Storage::enumerator enumerator;
typedef Storage::var_array var_array;
typedef Storage::real_array real_array;
typedef Storage::reference_array ref_array;

double Test10::Gravity() const {return -9.8;}
double Test10::FluidDensity() const {return 20;}
double Test10::FluidViscosity() const {return 0.1;}
double Test10::SolidDensity() const {return 120;}
double Test10::Porosity(double _x, double _y, double _z) const
{
	throw -1;
	return 0.0;
}
double Test10::InverseBiotModulus(double _x, double _y, double _z) const
{
	throw -1;
	return 0.0;
}



Test10::Test10()
//: trans_scale(0.0085)
	: trans_scale(1.0)
{
	WI[0] = 500000;
	WI[1] = 500000;
	WI[2] = 500000;
	pbhp[0] = 300000;
	pbhp[1] = 1000;
	pbhp[2] = 0;
}

vMatrix Test10::ElasticTensor(double _x, double _y, double _z) const
{
	throw -1;
	return rMatrix::Unit(6);
}

vMatrix Test10::BiotTensor(double _x, double _y, double _z) const
{
	throw -1;
	return rMatrix::Unit(3);
}

vMatrix Test10::PermTensor(double _x, double _y, double _z) const
{
	throw -1;
	return rMatrix::Unit(3);
}

hMatrix Test10::Displacement(double _x, double _y, double _z, double _t) const
{
	hMatrix sol(3,1);
	sol(0,0) = 0; //u
	sol(1,0) = 0; //v
	sol(2,0) = 0; //w
	return sol;
}

hessian_variable Test10::Pressure(double _x, double _y, double _z, double _t) const
{
	return 250;
}

rMatrix Test10::BCFlow_coef(double _x, double _y, double _z, double _t) const
{
	throw -1;
	rMatrix ret(2,1);
	ret(0,0) = 0; //dirichlet parallel
	ret(1,0) = 1; //neumann parallel
	return ret;
}

rMatrix Test10::BCMech_coef(double _x, double _y, double _z, double _t) const
{
	throw -1;
	rMatrix ret(4,1);
	ret(0,0) = 1; //dirichlet perp
	ret(1,0) = 0; //neumann perp
	ret(2,0) = 0; //dirichlet parallel
	ret(3,0) = 1; //neumann parallel
	return ret;
}

void Test10::Init(Mesh & m)
{
	//prepare wells
	real ccnt[3][3] =
	{
		{4.567151e+05, 7.321079e+06, 2.767665e+03},
		{4.609346e+05, 7.323503e+06, 2.597767e+03},
		{4.595400e+05, 7.326078e+06, 2.803586e+03}
	};
	TagReal tag_WI = m.CreateTag("WI",DATA_REAL,CELL,CELL,1);
	for(Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it) if( it->GetStatus() != Element::Ghost )
	{
		for(int q = 0; q < 3; ++q) if( it->Inside(ccnt[q]) )
		{
			cc[q] = it->self();
			std::cout << "proc " << m.GetProcessorRank() << " found c" << q << " " << cc[q].LocalID() << std::endl;
			tag_WI[cc[q]] = WI[q];
		}
	}
	
	if (m.HaveTag("PERM"))
	{
		TagRealArray perm = m.GetTag("PERM");
		int perm_size = perm.GetSize();
		for (Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it)
			for (int k = 0; k < perm_size; ++k) perm[*it][k] *= trans_scale;
	}
	rMatrix I(3,3), E(3,3),n(3,1),x(3,1);
	I = rMatrix::Unit(3);
	int bad = 0;
	for(Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it) if( it->GetStatus() != Element::Ghost )
	{
		ElementArray<Face> faces = it->getFaces();
		E.Zero();
		for(ElementArray<Face>::iterator jt = faces.begin(); jt != faces.end(); ++jt)
		{
			jt->OrientedNormal(it->self(),n.data());
			jt->Barycenter(x.data());
			E += x*n.Transpose();
		}
		E = (E/it->Volume()-I);
		if( E.FrobeniusNorm() > 1.0e-4 || it->Volume() < 0)
			bad++;
	}
	bad = m.Integrate(bad);
	if( m.GetProcessorRank() == 0 )
		std::cout << "Cells with bad divergence: " << bad << std::endl;
	
	int num_outside = 0, num_total = 0;
#if defined(USE_OMP)
#pragma omp parallel for reduction(+:num_outside) reduction(+:num_total)
#endif
	for (int it = 0; it < m.CellLastLocalID(); ++it) if (m.isValidCell(it))
	{
		Cell c = m.CellByLocalID(it);
		if( c.GetStatus() == Element::Ghost ) continue;
		Storage::real cnt[3];
		c->Barycenter(cnt);
		if (!c->Inside(cnt)) num_outside++;
		num_total++;
	}
	num_outside = m.Integrate(num_outside);
	num_total = m.Integrate(num_total);
	if( m.GetProcessorRank() == 0 )
		std::cout << "Barycenter is outside for: " << num_outside << " of " << num_total << " cells " << std::endl;
}

void Test10::SetBC(Mesh & m, double T, MarkerType boundary) const
{
	TagRealArray tag_BC_flow = m.GetTag("BOUNDARY_CONDITION_FLOW");
	TagRealArray tag_BC_mech = m.GetTag("BOUNDARY_CONDITION_ELASTIC");
	TagBulk      tag_BC_type = m.GetTag("BCFACEDIR");
	real BCf[3] = {0,1,0};
	//real BCm_bottom[7] = {1,0,1,0,0,0,0};
	real BCm_bottom[7] = { 1,0,0,1,0,0,0 };
	real BCm_side[7] = {1,0,0,1,0,0,0};
	//real BCm_top[7] = {0,1,0,1,0,0,0};
	real BCm_top[7] = { 1,0,0,1,0,0,0 };
	//real BCm_top[7] = { 1,0,1,0,0,0,0 };
	rMatrix mBCf(BCf, 3,1);
	rMatrix mBCm_bottom(BCm_bottom, 7,1);
	rMatrix mBCm_side(BCm_side, 7,1);
	rMatrix mBCm_top(BCm_top, 7,1);
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
			//std::cout << "Old bc:" << std::endl;
			//std::cout << "p "; tag_BC_flow(f,3,1).Transpose().Print();
			//std::cout << "u "; tag_BC_mech(f,7,1).Transpose().Print();
			tag_BC_flow(f,3,1) = mBCf;
			tag_BC_mech(f,7,1) = mBCm_bottom;
			if( tag_BC_type[f] == 0 ) //vertical face
				tag_BC_mech(f,7,1) = mBCm_side;
			else if( tag_BC_type[f] == 1 ) //top face
				tag_BC_mech(f,7,1) = mBCm_top;
			else if( tag_BC_type[f] == 2 ) //bottom face
				tag_BC_mech(f,7,1) = mBCm_bottom;
			//std::cout << "New bc:" << std::endl;
			//std::cout << "p "; tag_BC_flow(f,3,1).Transpose().Print();
			//std::cout << "u "; tag_BC_mech(f,7,1).Transpose().Print();
		}
	}
}

void Test10::SetForce(Mesh & m, const INMOST::dynamic_variable & p, double T) const
{
	TagVariableArray tag_F = m.GetTag("FORCE");
#if defined(USE_OMP)
#pragma omp parallel for
#endif
	for(int i = 0; i < m.CellLastLocalID(); ++i) if( m.isValidCell(i) )
	{
		Cell c = m.CellByLocalID(i);
		if( c.GetStatus() == Element::Ghost ) continue;
		if( c.HaveData(tag_F) ) c.DelData(tag_F);
		//tag_F(c,4,1).Zero();// = Force(cnt[0],cnt[1],cnt[2],T);
	}
	//c1 4.567151e+05 7.321079e+06 2.767665e+03 9018
	//c2 4.609346e+05 7.323503e+06 2.597767e+03 20756
	//c3 4.595400e+05 7.326078e+06 2.803586e+03 6698
	
	double mu = FluidViscosity();
	
	for(int q = 0; q < 3; ++q) if( cc[q].isValid() )
	{
		tag_F[cc[q]][3] -= WI[q]/mu*(p[cc[q]] - pbhp[q])/cc[q].Volume();
		//std::cout << "Pressure on " << cc[q].LocalID() << " is " << get_value(p[cc[q]]) << " bhp is " << pbhp[q] << " WI " << WI[q] << std::endl;
		//varPrint(tag_F[cc[q]][3]);
		//std::cout << std::endl;
	}
}

void Test10::SetProperty(Mesh & m) const
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
	for (int it = 0; it < m.CellLastLocalID(); ++it) if (m.isValidCell(it))
	{
		Cell c = m.CellByLocalID(it);
		double S0 = perm[c][0];
		double S1 = perm[c][1];
		double S2 = perm[c][2];
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
#if defined(USE_OMP)
#pragma omp for
#endif
		for(int it = 0; it < m.CellLastLocalID(); ++it) if( m.isValidCell(it) )
		{
			real cnt[3];
			Cell c = m.CellByLocalID(it);
			c.Barycenter(cnt);
			
			double S0 = perm[c][0];
			double S1 = perm[c][1];
			double S2 = perm[c][2];
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
						S0 += perm[kt->self()][0] * kt->Volume();
						S1 += perm[kt->self()][1] * kt->Volume();
						S2 += perm[kt->self()][2] * kt->Volume();
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
			
			
			// transversal isotropy:
			// free: E1, E3, nu12, nu13, G13
			// E2 = E1
			// nu23 = nu13
			// G13 = G23
			// G12 = 0.5*E1/(1+nu12)

			/*
			E1 = 0.1 + S0*100.0;
			E2 = 0.1 + S1*100.0;
			E3 = 0.1 + S2*100.0;
			nu12 = 0.15;
			nu13 = nu23 = 0.05;
			G23 = 1 + S0*150.0;
			G13 = 1 + S1*150.0;
			//if( E1 == E2 && nu13 == nu23 )
			//	G12 = 0.5 * E1 / (1 + nu12);
			//else
			G12 = 1 + S2*150.0;
			*/
			/*
			//works
			E1 = S0*400;
			E2 = S1*400;
			E3 = S2*400;
			G23 = 300 * S0;
			G13 = 300 * S1;
			G12 = 300 * S2;
			nu12 = 0.1 / E2 * E1;
			nu13 = 0.1 / E3 * E1;
			nu23 = 0.2 / E3 * E2;
			*/

			/*
			E1 = S0 * 1000;
			E2 = S1 * 1000;
			E3 = S2 * 1000;
			G23 = 500 * S0;
			G13 = 500 * S1;
			G12 = 500 * S2;
			nu12 = 0.1 / E2 * E1;
			nu13 = 0.1 / E3 * E1;
			nu23 = 0.2 / E3 * E2;
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
						std::cout << "nu12: " << nu12 << "; nu13: " << nu13 << "; nu23: nu13; " << std::endl;// << nu23 << std::endl;
						//std::cout << "nu21: " << nu21 << " nu31: " << nu31 << " nu32: " << nu32 << std::endl;
						std::cout << "G12: " << G12 << "; G13: " << G13 << "; G23: " << G23 << "; " << std::endl;

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
			//tag_M[c] = 0.5;
			//tag_M[c] = 2;
		}
	}
	std::cout << "Saving props.vtk" << std::endl;
	m.Save("props.vtk");

	{
		rMatrix I(3, 3), E(3, 3), n(3, 1), x(3, 1);
		I = rMatrix::Unit(3);
		int bad = 0;
		for (Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it) if (it->GetStatus() != Element::Ghost)
		{
			ElementArray<Face> faces = it->getFaces();
			E.Zero();
			for (ElementArray<Face>::iterator jt = faces.begin(); jt != faces.end(); ++jt)
			{
				jt->OrientedNormal(it->self(), n.data());
				jt->Barycenter(x.data());
				E += x * n.Transpose();
			}
			E = (E / it->Volume() - I);
			if (E.FrobeniusNorm() > 1.0e-4 || it->Volume() < 0)
				bad++;
		}
		bad = m.Integrate(bad);
		if (m.GetProcessorRank() == 0)
			std::cout << "Cells with bad divergence: " << bad << std::endl;

		int num_outside = 0, num_total = 0, negvol = 0;
#if defined(USE_OMP)
#pragma omp parallel for reduction(+:num_outside) reduction(+:num_total) reduction(+:negvol)
#endif
		for (int it = 0; it < m.CellLastLocalID(); ++it) if (m.isValidCell(it))
		{
			Cell c = m.CellByLocalID(it);
			if (c.GetStatus() == Element::Ghost) continue;
			Storage::real cnt[3];
			c->Barycenter(cnt);
			if (!c->Inside(cnt)) num_outside++;
			num_total++;
			if (c->Volume() < 0) negvol++;
		}
		num_outside = m.Integrate(num_outside);
		num_total = m.Integrate(num_total);
		negvol = m.Integrate(negvol);
		if (m.GetProcessorRank() == 0)
			std::cout << "Barycenter is outside for: " << num_outside << " of " << num_total << " cells, negative volume in " << negvol << " cells " << std::endl;
	}
}











