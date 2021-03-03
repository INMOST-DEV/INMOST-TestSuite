#include "test13.h"

using namespace INMOST;

//shortcut for types
typedef Storage::real real;
typedef Storage::enumerator enumerator;
typedef Storage::var_array var_array;
typedef Storage::real_array real_array;
typedef Storage::reference_array ref_array;

double Test13::Gravity() const {return -9.8;}
//double Test13::Gravity() const {return 0;}
double Test13::FluidDensity() const {return 20;}
double Test13::FluidViscosity() const {return 0.1;}
double Test13::SolidDensity() const {return 120;}
double Test13::Porosity(double _x, double _y, double _z) const
{
	throw -1;
	return 0.0;
}
double Test13::InverseBiotModulus(double _x, double _y, double _z) const
{
	throw -1;
	return 0.0;
}



Test13::Test13()
{
	WI[0] = 100;
	WI[1] = 100;
	//pbhp[0] = 265;
	pbhp[0] = 110;
	pbhp[1] = 110;
	//pbhp[1] = 265;
	
	fWI[0] = 100;
	fpbhp[0] = 475;
	
	real vE0[36] =
	{
		93,           46,           22,           13,           72,           35,
		46,           95,           41,           62,           56,           24,
		22,           41,           89,           25,           33,           21,
		13,           62,           25,           87,           13,           25,
		72,           56,           33,           13,           99,           57,
		35,           24,           21,           25,           57,           78
	};
	real vE1[36] =
	{
		81,           33,            5,           14,           37,           32,
		33,          100,           26,           58,           58,           73,
		5,           26,           90,           23,           60,           28,
		14,           58,           23,           77,           29,            3,
		37,           58,           60,           29,           84,           60,
		32,           73,           28,            3,           60,           99
		
	};
	real vK0[9] =
	{
		25,            2,           39,
		2,           42,            7,
		39,            7,          100
	};
	real vK1[9] =
	{
		1,            2,            2,
		2,           68,           42,
		2,           42,           49
	};
	real vB0[9] =
	{
		1,            6,            5,
		6,           67,           27,
		5,           27,           76
	};
	real vB1[9] =
	{
		1,            1,            1,
		1,           68,           13,
		1,           13,           73
	};
	//E0 = rMatrix(vE0,6,6)*0.01+rMatrix::Unit(6)*10;
	//E0 = rMatrix(vE0,6,6)/10.0;
	E0 = GenIsotropicTensor(30,0.15); //close to clay
	//K0 = rMatrix(vK0,3,3);
	K0 = rMatrix::Unit(3)*5;
	//B0 = rMatrix(vB0,3,3)/1000.0;
	B0 = rMatrix::Unit(3)*0.2;
	//E1 = rMatrix(vE1,6,6)*0.01+rMatrix::Unit(6)*10;
	//E1 = rMatrix(vE1,6,6)/10.0;
	E1 = GenIsotropicTensor(400,0.17); //close to
	//K1 = rMatrix(vK1,3,3);
	K1 = rMatrix::Unit(3)*5.5;
	//B1 = rMatrix(vB1,3,3)/1000.0;
	B1 = rMatrix::Unit(3)*1.35;
	
	Kf = (K0+K1)*0.5*100;
	
	//rCheckEigen(K1,true);
	
	std::cout << "B0: " << std::endl;
	B0.Print();
	rCheckEigen(B0,true);
	
	std::cout << "B1: " << std::endl;
	B1.Print();
	rCheckEigen(B1,true);
	
	phi0 = 0.25;
	phi1 = phi0;//0.35;
	M0 = 1.0;
	M1 = M0;//1.5;
	Mf = (M0+M1)*0.5*10;
}

vMatrix Test13::ElasticTensor(double _x, double _y, double _z) const
{
	throw -1;
	return rMatrix::Unit(6);
}

vMatrix Test13::BiotTensor(double _x, double _y, double _z) const
{
	throw -1;
	return rMatrix::Unit(3);
}

vMatrix Test13::PermTensor(double _x, double _y, double _z) const
{
	throw -1;
	return rMatrix::Unit(3);
}

hMatrix Test13::Displacement(double _x, double _y, double _z, double _t) const
{
	hMatrix sol(3,1);
	sol(0,0) = 0; //u
	sol(1,0) = 0; //v
	sol(2,0) = 0; //w
	return sol;
}

hessian_variable Test13::Pressure(double _x, double _y, double _z, double _t) const
{
	return 250;
}

rMatrix Test13::BCFlow_coef(double _x, double _y, double _z, double _t) const
{
	throw -1;
	rMatrix ret(2,1);
	ret(0,0) = 0; //dirichlet parallel
	ret(1,0) = 1; //neumann parallel
	return ret;
}

rMatrix Test13::BCMech_coef(double _x, double _y, double _z, double _t) const
{
	throw -1;
	rMatrix ret(4,1);
	ret(0,0) = 1; //dirichlet perp
	ret(1,0) = 0; //neumann perp
	ret(2,0) = 0; //dirichlet parallel
	ret(3,0) = 1; //neumann parallel
	return ret;
}

void Test13::Init(Mesh & m)
{
	//prepare wells
	darcy = m.CreateMarker();
	real ccnt[2][3] =
	{
		{0.125, 0.5, 0.5},
		{0.875, 0.5, 0.5}
	};
	real fcnt[1][3] =
	{
		{0.5,0.5,0.5}
		//{0.5,0.5,100}
	};
	real cnt[3] = {0,0,0};
	TagInteger tag_well_mark = m.CreateTag("WELL_MARK",DATA_INTEGER,CELL|FACE,NONE,1);
	TagReal tag_WI = m.CreateTag("WI",DATA_REAL,CELL|FACE,CELL|FACE,1);
	for(Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it) 
	{
		it->Centroid(cnt);
		//if( cnt[0] > 0.65 )
		//if( cnt[0] > 0.5 )
		//	it->SetMarker(darcy);
		if( it->GetStatus() != Element::Ghost )
		{
			for(int q = 0; q < 2; ++q) if( it->Inside(ccnt[q]) )
			{
				tag_well_mark[*it] = 1;
				cc[q] = it->self();
				std::cout << "proc " << m.GetProcessorRank() << " found c" << q << " " << cc[q].LocalID() << std::endl;
				tag_WI[cc[q]] = WI[q];
			}
		}
	}
	int nfracs = 0;
	for(Mesh::iteratorFace it = m.BeginFace(); it != m.EndFace(); ++it)
	{
		it->Centroid(cnt);
		if( 
			//false
			(fabs(cnt[0]-0.5) < 1.0e-5 && fabs(cnt[1] - 0.5) < 0.3 && fabs(cnt[2] - 0.5) < 0.2)
			//||
			//(fabs(cnt[0]-0.5) < 0.2 && fabs(cnt[1] - 0.4) < 1.0e-5 && fabs(cnt[2] - 0.5) < 0.3)
			)
		{
			it->SetMarker(darcy);
			std::cout << "face " << it->LocalID() << " cnt " << cnt[0] << " " << cnt[1] << " " << cnt[2] << std::endl;
			nfracs++;
		}
		if( it->GetStatus() != Element::Ghost )
		{
			for(int q = 0; q < 1; ++q) //if( it->Inside(fcnt[q]) )
			{
				double l = 0;
				for(int g = 0; g < 3; ++g) l += (cnt[g]-fcnt[q][g])*(cnt[g]-fcnt[q][g]);
				l = sqrt(l);
				if( l < 1.0e-7 )
				{
					tag_well_mark[*it] = 1;
					fc[q] = it->self();
					std::cout << "proc " << m.GetProcessorRank() << " found f" << q << " " << fc[q].LocalID() << std::endl;
					tag_WI[fc[q]] = fWI[q];
				}
			}
		}
	}
	std::cout << "fracture faces: " << nfracs << std::endl;
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
}

void Test13::SetBC(Mesh & m, double T, MarkerType boundary) const
{
	TagRealArray tag_BC_flow = m.GetTag("BOUNDARY_CONDITION_FLOW");
	TagRealArray tag_BC_mech = m.GetTag("BOUNDARY_CONDITION_ELASTIC");
	real BCf[3] = {0,1,0};
	real BCm_bottom[7] = {1,0,1,0,0,0,0};
	real BCm_side[7] = {1,0,0,1,0,0,0};
	real BCm_top[7] = {0,1,0,1,0,0,0};
	real BCm_fix[7] = {1,0,1,0,0,0,0};
	rMatrix mBCf(BCf, 3,1);
	rMatrix mBCm_bottom(BCm_bottom, 7,1);
	rMatrix mBCm_side(BCm_side, 7,1);
	rMatrix mBCm_top(BCm_top, 7,1);
	rMatrix mBCm_fix(BCm_fix, 7,1);
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
			Cell c = f.BackCell();
			//std::cout << "Old bc:" << std::endl;
			//std::cout << "p "; tag_BC_flow(f,3,1).Transpose().Print();
			//std::cout << "u "; tag_BC_mech(f,7,1).Transpose().Print();
			tag_BC_flow(f,3,1) = mBCf;
			if( !c.GetMarker(darcy) )
			{
				tag_BC_mech(f,7,1) = mBCm_fix;
				/*
				tag_BC_mech(f,7,1) = mBCm_bottom;
				if( cnt[2] > 1.0-1.0e-9 ) //top face
					tag_BC_mech(f,7,1) = mBCm_top;
					//tag_BC_mech(f,7,1) = mBCm_fix;
				else if( cnt[2] < 0.0+1.0e-9 ) //bottom face
					tag_BC_mech(f,7,1) = mBCm_bottom;
				else
					tag_BC_mech(f,7,1) = mBCm_side;
				*/
			}
			//std::cout << "New bc:" << std::endl;
			//std::cout << "p "; tag_BC_flow(f,3,1).Transpose().Print();
			//std::cout << "u "; tag_BC_mech(f,7,1).Transpose().Print();
		}
	}
}

void Test13::SetForce(Mesh & m, const INMOST::dynamic_variable & p, double T) const
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
	}
	
	double mu = FluidViscosity();
	
	for(int q = 0; q < 2; ++q) if( cc[q].isValid() )
	{
		if( cc[q].GetMarker(darcy) ) continue;
		tag_F[cc[q]][3] -= WI[q]/mu*(p[cc[q]] - pbhp[q])/cc[q].Volume();
		std::cout << "Pressure on " << cc[q].LocalID() << " is " << get_value(p[cc[q]]) << " bhp is " << pbhp[q] << " WI " << WI[q] << " mu " << mu << " volume " << cc[q].Volume() << " ";
		varPrint(tag_F[cc[q]][3]);
		std::cout << " (poromech) " << std::endl;
	}
}

void Test13::SetForceDarcy(Mesh & m, const INMOST::dynamic_variable & p,  double T) const
{
	TagVariableArray tag_FD = m.GetTag("FORCE_DARCY");
#if defined(USE_OMP)
#pragma omp parallel for
#endif
	for(int i = 0; i < m.CellLastLocalID(); ++i) if( m.isValidCell(i) )
	{
		Cell c = m.CellByLocalID(i);
		if( c.GetStatus() == Element::Ghost ) continue;
		if( c.HaveData(tag_FD) ) c.DelData(tag_FD);
	}
	
	double mu = FluidViscosity();
	
	for(int q = 0; q < 2; ++q) if( cc[q].isValid() )
	{
		if( !cc[q].GetMarker(darcy) ) continue;
		tag_FD[cc[q]][0] -= WI[q]/mu*(p[cc[q]] - pbhp[q])/cc[q].Volume();
		std::cout << "Pressure on " << cc[q].LocalID() << " is " << get_value(p[cc[q]]) << " bhp is " << pbhp[q] << " WI " << WI[q] << " mu " << mu << " volume " << cc[q].Volume() << " ";
		varPrint(tag_FD[cc[q]][0]);
		std::cout << " (darcy) " << std::endl;
	}
	
}

void Test13::SetForceDarcyFracture(Mesh & m, const INMOST::dynamic_variable & p, const INMOST::dynamic_variable & af, double T) const
{
	TagVariableArray tag_FD = m.GetTag("FORCE_DARCY");
#if defined(USE_OMP)
#pragma omp parallel for
#endif
	for(int i = 0; i < m.FaceLastLocalID(); ++i) if( m.isValidFace(i) )
	{
		Face f = m.FaceByLocalID(i);
		if( f.GetStatus() == Element::Ghost ) continue;
		if( f.HaveData(tag_FD) ) f.DelData(tag_FD);
	}
	
	double mu = FluidViscosity();
	
	for(int q = 0; q < 1; ++q) if( fc[q].isValid() )
	{
		if( !fc[q].GetMarker(darcy) ) continue;
		tag_FD[fc[q]][0] -= fWI[q]/mu*(p[fc[q]] - fpbhp[q])/(fc[q].Area()*af[fc[q]]);
		std::cout << "Pressure on face " << fc[q].LocalID() << " is " << get_value(p[fc[q]]) << " bhp is " << fpbhp[q] << " WI " << fWI[q] << " mu " << mu << " volume " << get_value(fc[q].Area()*af[fc[q]]) << " ";
		varPrint(tag_FD[fc[q]][0]);
		std::cout << " (fracture darcy) " << std::endl;
	}
}

void Test13::SetProperty(Mesh & m) const
{
	//this should be calculated
	TagRealArray tag_K = m.GetTag("PERM");
	TagRealArray tag_E = m.GetTag("ELASTIC_TENSOR");
	TagRealArray tag_B = m.GetTag("BIOT_COEFFICIENT");
	TagReal tag_M    = m.GetTag("INVERSE_BIOT_MODULUS");
	TagReal tag_poro = m.GetTag("PORO");
#if defined(USE_OMP)
#pragma omp parallel
#endif
	{
#if defined(USE_OMP)
#pragma omp for
#endif
		for(int it = 0; it < m.CellLastLocalID(); ++it) if( m.isValidCell(it) )
		{
			real cnt[3];
			Cell c = m.CellByLocalID(it);
			c.Barycenter(cnt);
			
			if( cnt[0] < 0.5 )
			{
				if( tag_K.GetSize() == 6 )
					raSymmetricMatrixMake(tag_K[c].data(),3) = K0;
				else if( tag_K.GetSize() == 9 )
					tag_K(c,3,3) = K0;
				else std::cout << __FILE__ << ":" << __LINE__ << " PERM has size " << tag_K.GetSize() << std::endl;
				if( tag_E.GetSize() == 21 )
					raSymmetricMatrixMake(tag_E[c].data(),6) = E0;
				else if( tag_E.GetSize() == 36 )
					tag_E(c,6,6) = E0;
				else std::cout << __FILE__ << ":" << __LINE__ << " ELASTIC_TENSOR has size " << tag_E.GetSize() << std::endl;
				if( tag_B.GetSize() == 6 )
					raSymmetricMatrixMake(tag_B[c].data(),3) = B0;
				else if( tag_B.GetSize() == 9 )
					tag_B(c,3,3) = B0;
				else std::cout << __FILE__ << ":" << __LINE__ << " BIOT_COEFFICIENT has size " << tag_B.GetSize() << std::endl;
				tag_M[c] = M0;
				tag_poro[c] = phi0;
			}
			else
			{
				if( tag_K.GetSize() == 6 )
					raSymmetricMatrixMake(tag_K[c].data(),3) = K1;
				else if( tag_K.GetSize() == 9 )
					tag_K(c,3,3) = K1;
				else std::cout << __FILE__ << ":" << __LINE__ << " PERM has size " << tag_K.GetSize() << std::endl;
				if( tag_E.GetSize() == 21 )
					raSymmetricMatrixMake(tag_E[c].data(),6) = E1;
				else if( tag_E.GetSize() == 36 )
					tag_E(c,6,6) = E1;
				else std::cout << __FILE__ << ":" << __LINE__ << " ELASTIC_TENSOR has size " << tag_E.GetSize() << std::endl;
				if( tag_B.GetSize() == 6 )
					raSymmetricMatrixMake(tag_B[c].data(),3) = B1;
				else if( tag_B.GetSize() == 9 )
					tag_B(c,3,3) = B1;
				else std::cout << __FILE__ << ":" << __LINE__ << " BIOT_COEFFICIENT has size " << tag_B.GetSize() << std::endl;
				tag_M[c] = M1;
				tag_poro[c] = phi1;
			}
	
		}
#if defined(USE_OMP)
#pragma omp for
#endif
		for(int it = 0; it < m.FaceLastLocalID(); ++it) if( m.isValidFace(it) )
		{
			real cnt[3];
			Face f = m.FaceByLocalID(it);
			if( f.GetMarker(darcy) )
			{
				if( tag_K.GetSize() == 6 )
					raSymmetricMatrixMake(tag_K[f].data(),3) = Kf;
				else if( tag_K.GetSize() == 9 )
					tag_K(f,3,3) = Kf;
				else std::cout << __FILE__ << ":" << __LINE__ << " PERM has size " << tag_K.GetSize() << std::endl;
				
				if( tag_M.isDefined(FACE) && Mf != 0 )
				{
					std::cout << "set " << Mf << std::endl;
					tag_M[f] = Mf;
				}
			}
		}
	}
	std::cout << "Saving props.vtk" << std::endl;
	m.Save("props.vtk");
}












