#include "abstract_test.h"
//#include "options.h"
using namespace INMOST;

//shortcut for types
typedef Storage::real real;
typedef Storage::enumerator enumerator;
typedef Storage::var_array var_array;
typedef Storage::real_array real_array;
typedef Storage::reference_array ref_array;



void SVD2Eigen(const rMatrix & U, rMatrix & S, rMatrix & V)
{
	for (int i = 0; i < V.Cols(); ++i)
	{
		double dot = 0.0;
		for (int j = 0; j < V.Rows(); ++j)
			dot += U(j, i)*V(j, i);
		if (dot < 0.0)
		{
			S(i, i) *= -1;
			for (int j = 0; j < V.Rows(); ++j)
				V(j, i) *= -1;
		}
	}
	//check
	//if ((U - V).FrobeniusNorm() > 1.0e-8)
	//	(U - V).Print();
}


bool CheckEigen(const rMatrix & U, rMatrix & S, rMatrix & V, bool print)
{
	//SVD2Eigen
	SVD2Eigen(U,S,V);
	bool negative = false;
	int iend = S.Rows();
	double sum = 0;
	for (int i = 0; i < iend; ++i) sum += fabs(S(i, i));
	for (int i = 0; i < iend; ++i) if (S(i, i) < -1.0e-9*sum) negative = true;

	//if (negative && print)
	if( print )
	{
		for (int i = 0; i < iend; ++i)
			std::cout << S(i, i) << " ";
		std::cout << std::endl;
	}

	return !negative;
}

template<typename T>
bool CheckEigen(const Matrix<T> & A, bool print)
{
	rMatrix Aval = A, U, S, V;
	Aval.SVD(U,S,V);
	return CheckEigen(U,S,V,print);
}

bool rCheckEigen(const rMatrix & A, bool print) {return CheckEigen(A,print);}
bool vCheckEigen(const vMatrix & A, bool print) {return CheckEigen(A,print);}

int GenTensor(rMatrix & C, int size, int seed, bool print)
{
	rMatrix U(size, size), S(size, size), V(size, size);
	C.Resize(size,size);
	for (int k = seed; k < INT_MAX; ++k)
	{
		std::cout << k << "\r";
		std::cout.flush();
		srand(k);
		//construct some symmetric C
		for (int i = 0; i < size; ++i)
		for (int j = i; j < size; ++j)
			C(i, j) = C(j, i) = ceil(rand() / (1.*RAND_MAX)*100.0);
		//perform SVD
		C.SVD(U, S, V);
		//check column of U is the same sign as row of V^T
		if (CheckEigen(U, S, V,print))
		{
			SVD2Eigen(U, S, V);
			std::cout << "eigenvalues: ";
			for (int q = 0; q < size; ++q)
				std::cout << S(q, q) << " ";
			std::cout << std::endl;
			return k;
		}
	}
	return -1;
}




rMatrix GenMatrix1(int size)
{
	bool print = false;
	rMatrix C(size,size), U(size,size), S(size,size), V(size,size);
	for(int i = 0; i < size; ++i)
		for(int j = i; j < size; ++j)
			C(i,j) = C(j,i) = rand()/(1.*RAND_MAX)*2.0-1.0;
	if( print )
	{
		std::cout << "Original matrix: " << std::endl;
		C.Print();
	}
	//perform SVD
	C.SVD(U,S,V);
	//check column of U is the same sign as row of V^T
	//SVD2Eigen(U,S,V);
	for(int i = 0; i < V.Cols(); ++i)
	{
		double dot = 0.0;
		for(int j = 0; j < V.Rows(); ++j)
			dot += U(j,i)*V(j,i);
		if( dot < 0.0 )
		{
			S(i,i) *= -1;
			for(int j = 0; j < V.Rows(); ++j)
				V(j,i) *= -1;
		}
	}
	if( print )
	{
		std::cout << "Original U: " << std::endl;
		U.Print();
		std::cout << "Original S: " << std::endl;
		S.Print();
		std::cout << "Original V: " << std::endl;
		V.Print();
		std::cout << "Original U-V: " << std::endl;
		(U-V).Print();
	}
	//make singular values positive
	for(int i = 0; i < size; ++i)
		if( S(i,i) < 0.0 )
			S(i,i) = std::abs(S(i,i)) + 1.0e-5;
	if( print )
	{
		std::cout << "Modified S: " << std::endl;
		S.Print();
	}
	//assemble matrix back
	C = U*S*V.Transpose();
	
	if( print )
	{
		std::cout << "Matrix: " << std::endl;
		C.Print();
	}
	
	if( print )
	{
		C.SVD(U,S,V);
		std::cout << "Modified U: " << std::endl;
		U.Print();
		std::cout << "Modified S: " << std::endl;
		S.Print();
		std::cout << "Modified V: " << std::endl;
		V.Print();
		std::cout << "Modified U-V: " << std::endl;
		(U-V).Print();
	}
	CheckEigen(C,true);
	return C;
}


rMatrix AbstractTest::Flux(double _x, double _y, double _z, double _t, double nrm[3], bool bndcond) const
{
	vMatrix sol = Displacement(_x,_y,_z,_t);
	rMatrix N(3,6), flux(4,1), n(nrm,3,1);
	rMatrix C = ElasticTensor(_x,_y,_z);
	rMatrix B = BiotTensor(_x,_y,_z);
	rMatrix K = PermTensor(_x,_y,_z);
	rMatrix epsilon(6,1);
	variable p = Pressure(_x,_y,_z,_t);
	rMatrix gradp(3,1);
	real mu = FluidViscosity();
	rMatrix soldt(3,1);
	rMatrix z(3,1);
	real g = Gravity();
	real rho_f = FluidDensity();
	real rho_s = SolidDensity();
	real phi = Porosity(_x,_y,_z);
	real rho = rho_f*phi+rho_s*(1-phi);
	z(0,0) = 0;
	z(1,0) = 0;
	z(2,0) = 1;
	
	
	epsilon(0,0) = sol(0,0).GetRow()[0]; //u_x
	epsilon(1,0) = sol(1,0).GetRow()[1]; //v_y
	epsilon(2,0) = sol(2,0).GetRow()[2]; //w_z
	epsilon(3,0) = sol(1,0).GetRow()[2] + sol(2,0).GetRow()[1]; //v_z + w_y
	epsilon(4,0) = sol(0,0).GetRow()[2] + sol(2,0).GetRow()[0]; //u_z + w_x
	epsilon(5,0) = sol(0,0).GetRow()[1] + sol(1,0).GetRow()[0]; //u_y + v_x
	

	//we need N*sigma
	//     | nx           nz ny |
	// N = |    ny     nz    nx |
	//     |       nz  ny nx    |
	N(0,0) = nrm[0];
	N(1,1) = nrm[1];
	N(2,2) = nrm[2];
	
	N(0,4) = nrm[2];
	N(0,5) = nrm[1];
	
	N(1,3) = nrm[2];
	N(1,5) = nrm[0];
	
	N(2,3) = nrm[1];
	N(2,4) = nrm[0];
	
	flux(0,3,0,1) = N*C*epsilon - B*p.GetValue()*n;
	
	if( gravity_in_flux )
		flux(0,3,0,1) += rho*g*_z*n;
	
	gradp(0,0) = p.GetRow()[0];
	gradp(1,0) = p.GetRow()[1];
	gradp(2,0) = p.GetRow()[2];
	
	
	
	soldt(0,0) = sol(0,0).GetRow()[3];
	soldt(1,0) = sol(1,0).GetRow()[3];
	soldt(2,0) = sol(2,0).GetRow()[3];
	
	
	flux(3,0) = (gradp-rho_f*g*z).DotProduct(K*n)/mu;
	if(!bndcond) flux(3,0)-=soldt.DotProduct(B*n);
	
	
	return flux;
}


rMatrix AbstractTest::Force(double _x, double _y, double _z, double _t) const
{
	hMatrix sol = Displacement(_x,_y,_z,_t);
	hessian_variable p = Pressure(_x,_y,_z,_t);
	vMatrix C = ElasticTensor(_x,_y,_z);
	vMatrix epsilon(6,1), sigma(6,1);
	vMatrix B = BiotTensor(_x,_y,_z);
	vMatrix K = PermTensor(_x,_y,_z);
	vMatrix w(6,1);
	rMatrix force(4,1);
	vMatrix gradp(3,1);
	vMatrix soldt(3,1);
	vMatrix u(3,1);
	rMatrix z(3,1);
	real g = Gravity();
	real mu = FluidViscosity();
	real im = InverseBiotModulus(_x,_y,_z);
	real rho_f = FluidDensity();
	real rho_s = SolidDensity();
	real phi = Porosity(_x,_y,_z);
	real rho = rho_f*phi+rho_s*(1-phi);
	real pdt = p.GetRow()[3];
	z(0,0) = 0;
	z(1,0) = 0;
	z(2,0) = 1;
	
	w(0,0) = B(0,0);
	w(1,0) = B(1,1);
	w(2,0) = B(2,2);
	w(3,0) = B(1,2);
	w(4,0) = B(0,2);
	w(5,0) = B(0,1);
	
	epsilon(0,0) = sol(0,0).GetVariable(0); //u_x
	epsilon(1,0) = sol(1,0).GetVariable(1); //v_y
	epsilon(2,0) = sol(2,0).GetVariable(2); //w_z
	epsilon(3,0) = sol(1,0).GetVariable(2) + sol(2,0).GetVariable(1); //v_z + w_y
	epsilon(4,0) = sol(0,0).GetVariable(2) + sol(2,0).GetVariable(0); //u_z + w_x
	epsilon(5,0) = sol(0,0).GetVariable(1) + sol(1,0).GetVariable(0); //u_y + v_x
	
	sigma = C*epsilon - w*variable(p);
	
	//sigma layot
	// s0 s5 s4
	// s5 s1 s3
	// s4 s3 s2
	// divergence is taken row-wise
	force(0,0) = -(sigma(0,0).GetRow()[0] + sigma(5,0).GetRow()[1] + sigma(4,0).GetRow()[2]);
	force(1,0) = -(sigma(5,0).GetRow()[0] + sigma(1,0).GetRow()[1] + sigma(3,0).GetRow()[2]);
	force(2,0) = -(sigma(4,0).GetRow()[0] + sigma(3,0).GetRow()[1] + sigma(2,0).GetRow()[2]);
	force(0,3,0,1) -= rho*g*z;
	
	gradp(0,0) = p.GetVariable(0);
	gradp(1,0) = p.GetVariable(1);
	gradp(2,0) = p.GetVariable(2);
	
	soldt(0,0) = sol(0,0).GetVariable(3);
	soldt(1,0) = sol(1,0).GetVariable(3);
	soldt(2,0) = sol(2,0).GetVariable(3);
	
	u = K*(gradp-rho_f*g*z)/mu - B*soldt;
	
	force(3,0) = im*pdt - (u(0,0).GetRow()[0] + u(1,0).GetRow()[1] + u(2,0).GetRow()[2]);
	
	//force.Transpose().Print();
	
	return force;
}

rMatrix AbstractTest::Stress(double x, double y, double z, double t, bool full) const
{
	vMatrix sol = Displacement(x,y,z,t);
	rMatrix G(3,3), e(3,3), ev(6,1);
	G(0,0) = sol(0,0).GetRow()[0];
	G(0,1) = sol(0,0).GetRow()[1];
	G(0,2) = sol(0,0).GetRow()[2];
	G(1,0) = sol(1,0).GetRow()[0];
	G(1,1) = sol(1,0).GetRow()[1];
	G(1,2) = sol(1,0).GetRow()[2];
	G(2,0) = sol(2,0).GetRow()[0];
	G(2,1) = sol(2,0).GetRow()[1];
	G(2,2) = sol(2,0).GetRow()[2];
	e = (G+G.Transpose())*0.5;
	ev(0,0) = e(0,0);
	ev(1,0) = e(1,1);
	ev(2,0) = e(2,2);
	ev(3,0) = e(1,2)*2;
	ev(4,0) = e(0,2)*2;
	ev(5,0) = e(0,1)*2;
	rMatrix C = ElasticTensor(x,y,z);
	if( full )
	{
		rMatrix B = BiotTensor(x,y,z), w(6,1);
		w(0,0) = B(0,0);
		w(1,0) = B(1,1);
		w(2,0) = B(2,2);
		w(3,0) = B(1,2);
		w(4,0) = B(0,2);
		w(5,0) = B(0,1);
		real p = Pressure(x,y,z,t).GetValue();
		return C*ev - w*p;
	}
	return C*ev;// - w*p;
}

rMatrix AbstractTest::BCFlow(double _x, double _y, double _z, double _t, double nrm[3]) const
{
	rMatrix ret(3,1);
	//ret(0,0) = 1; //dirichlet
	//ret(1,0) = 0;
	ret(0,2,0,1) = BCFlow_coef(_x,_y,_z,_t);
	ret(2,0) = ret(0,0)*Pressure(_x,_y,_z,_t).GetValue() + ret(1,0)*Flux(_x,_y,_z,_t,nrm,true)(3,0);
	return ret;
}

rMatrix AbstractTest::BCMech(double _x, double _y, double _z, double _t, double nrm[3]) const
{
	rMatrix ret(7,1), n(nrm,3,1), I = rMatrix::Unit(3), N = n*n.Transpose();
	//ret(0,0) = 1; //dirichlet perp
	//ret(1,0) = 0;
	//ret(2,0) = 1; //dirichlet parallel
	//ret(3,0) = 0;
	ret(0,4,0,1) = BCMech_coef(_x,_y,_z,_t);
	ret(4,7,0,1) = (N*ret(0,0)+(I-N)*ret(2,0))*Displacement(_x,_y,_z,_t) + (N*ret(1,0)+(I-N)*ret(3,0))*Flux(_x,_y,_z,_t,nrm,true)(0,3,0,1);
	return ret;
}


void AbstractTest::SetForce(Mesh & m, const INMOST::dynamic_variable & p, double T) const
{
	if( func_output && !m.GetProcessorRank() ) std::cout << __FUNCTION__ << " at time = " << T << std::endl;
	Automatizator * aut = Automatizator::GetCurrent();
	Automatizator::RemoveCurrent();
	TagVariableArray tag_F = m.GetTag("FORCE");
#if defined(USE_OMP)
#pragma omp parallel for
#endif
	for(int i = 0; i < m.CellLastLocalID(); ++i) if( m.isValidCell(i) )
	{
		real cnt[3];
		Cell c = m.CellByLocalID(i);
		if( c.GetStatus() == Element::Ghost ) continue;
		c.Barycenter(cnt);
		tag_F(c,4,1) = Force(cnt[0],cnt[1],cnt[2],T);
	}
	Automatizator::MakeCurrent(aut);
}


void AbstractTest::SetBC(Mesh & m, double T,MarkerType boundary) const
{
	if( func_output && !m.GetProcessorRank()) std::cout << __FUNCTION__ << " at time = " << T << std::endl;
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
			//std::cout << "Old bc:" << std::endl;
			//std::cout << "p "; tag_BC_flow(f,3,1).Transpose().Print();
			//std::cout << "u "; tag_BC_mech(f,7,1).Transpose().Print();
			tag_BC_flow(f,3,1) = BCFlow(cnt[0],cnt[1],cnt[2],T,nrm);
			tag_BC_mech(f,7,1) = BCMech(cnt[0],cnt[1],cnt[2],T,nrm);
			//std::cout << "New bc:" << std::endl;
			//std::cout << "p "; tag_BC_flow(f,3,1).Transpose().Print();
			//std::cout << "u "; tag_BC_mech(f,7,1).Transpose().Print();
		}
	}
	Automatizator::MakeCurrent(aut);
}

void AbstractTest::SetInitial(Mesh & m,double T, double Told, MarkerType orient) const
{
	if( func_output && !m.GetProcessorRank()) std::cout << __FUNCTION__ << " at time = " << T << " old time " << Told << std::endl;
	Automatizator * aut = Automatizator::GetCurrent();
	Automatizator::RemoveCurrent();
	if( m.HaveTag("UVWP") )
	{
		TagRealArray tag_UVWP = m.GetTag("UVWP");
		TagRealArray tag_UVWP0 = m.GetTag("UVWP0");
#if defined(USE_OMP)
#pragma omp parallel for
#endif
		for(int it = 0; it < m.CellLastLocalID(); ++it) if( m.isValidCell(it) )
		{
			real cnt[3];
			Cell c = m.CellByLocalID(it);
			c.Barycenter(cnt);
			tag_UVWP0[c][3] = tag_UVWP[c][3] = Pressure(cnt[0],cnt[1],cnt[2],T).GetValue();
			tag_UVWP0(c,4,1)(0,3,0,1) = tag_UVWP(c,4,1)(0,3,0,1) = Displacement(cnt[0],cnt[1],cnt[2],T);
		}
		if( m.HaveTag("PF") )
		{
			MarkerType darcy = DarcyDomain();
			TagReal tag_PF = m.GetTag("PF");
			for(Mesh::iteratorElement it = m.BeginElement(CELL|FACE); it != m.EndElement(); ++it)
			{
				real cnt[3];
				if( !it->GetMarker(darcy) ) continue;
				it->Barycenter(cnt);
				tag_PF[*it] = Pressure(cnt[0],cnt[1],cnt[2],T).GetValue();
			}
		}
		//std::cout << __FILE__ << ":" << __LINE__ << std::endl;
		//m.Save("init.vtk");
		return;
	}
	TagReal tag_P, tag_Pold, tag_nU;
	TagRealArray tag_UVW, tag_UVWold;
	tag_P    = m.GetTag("P");
	if( staggered )
		tag_nU   = m.GetTag("nU");
	else
		tag_UVW  = m.GetTag("UVW");
	if( bdf2 )
	{
		tag_Pold  = m.GetTag("Pold");
		tag_UVWold= m.GetTag("UVWold");
	}
#if defined(USE_OMP)
#pragma omp parallel for
#endif
	for(int it = 0; it < m.CellLastLocalID(); ++it) if( m.isValidCell(it) )
	{
		real cnt[3];
		Cell c = m.CellByLocalID(it);
		c.Barycenter(cnt);
		//initial state,pressure
		{
			tag_P[c] = Pressure(cnt[0],cnt[1],cnt[2],T).GetValue();
			if( !staggered )
				tag_UVW(c,3,1) = Displacement(cnt[0],cnt[1],cnt[2],T);
			if (bdf2)
			{
				tag_Pold[c] = Pressure(cnt[0], cnt[1], cnt[2], Told).GetValue();
				tag_UVWold(c, 3, 1) = Displacement(cnt[0], cnt[1], cnt[2], Told);
			}
		}
	}
	if( staggered )
	{
#if defined(USE_OMP)
#pragma omp parallel for
#endif
		for(int it = 0; it < m.FaceLastLocalID(); ++it) if( m.isValidFace(it) )
		{
			//initial state, displacement
			Face f = m.FaceByLocalID(it);
			rMatrix UVW(3,1);
			real nrm[3], cnt[3];
			real s = f.GetMarker(orient) ? -1 : 1;
			f.Barycenter(cnt);
			f.UnitNormal(nrm);
			UVW = Displacement(cnt[0],cnt[1],cnt[2],T);
			
			tag_nU[f] = s*(nrm[0]*UVW(0,0) + nrm[1]*UVW(1,0) + nrm[2]*UVW(2,0));
		}
	}
	Automatizator::MakeCurrent(aut);
	//std::cout << __FILE__ << ":" << __LINE__ << std::endl;
	//m.Save("init.vtk");
}


void AbstractTest::SetProperty(Mesh & m) const
{
	//media property, tensors
	TagRealArray tag_E = m.GetTag("ELASTIC_TENSOR");
	TagRealArray tag_B = m.GetTag("BIOT_COEFFICIENT");
	TagRealArray tag_K = m.GetTag("PERM");
	//media property, scalars
	TagReal tag_PORO = m.GetTag("PORO");
	TagReal tag_M    = m.GetTag("INVERSE_BIOT_MODULUS");
#if defined(USE_OMP)
#pragma omp parallel for
#endif
	for(int it = 0; it < m.CellLastLocalID(); ++it) if( m.isValidCell(it) )
	{
		real cnt[3];
		Cell c = m.CellByLocalID(it);
		c.Barycenter(cnt);
		if( tag_K.GetSize() == 9 )
			tag_K(c,3,3) = PermTensor(cnt[0],cnt[1],cnt[2]);
		else if( tag_K.GetSize() == 6 )
			raSymmetricMatrixMake(tag_K[c].data(),3) = PermTensor(cnt[0],cnt[1],cnt[2]);
		tag_B(c,3,3) = BiotTensor(cnt[0],cnt[1],cnt[2]);
		if( tag_E.GetSize() == 36 )
			tag_E(c,6,6) = ElasticTensor(cnt[0],cnt[1],cnt[2]);
		else if( tag_E.GetSize() == 21 )
			raSymmetricMatrixMake(tag_E[c].data(),6) = ElasticTensor(cnt[0],cnt[1],cnt[2]);
		tag_PORO[c] = Porosity(cnt[0],cnt[1],cnt[2]);
		tag_M[c] = InverseBiotModulus(cnt[0],cnt[1],cnt[2]);
	}
}


void Report(const AbstractTest * test,Mesh & m, double T, double dT, double time_shift, MarkerType orient)
{
	
	//reconstruct stress
	TagRealArray tag_E = m.GetTag("ELASTIC_TENSOR");
	TagRealArray tag_B = m.GetTag("BIOT_COEFFICIENT");
	TagRealArray tag_UVW;
	TagRealArray tag_FLUX  = m.GetTag("Fq");
	TagRealArray tag_UVWPF =  m.GetTag("UVWPF");
	TagRealArray tag_S = m.CreateTag("STRESS",DATA_REAL,CELL,NONE,6);
	TagVariableArray tag_rUVW;
	
	if( test->staggered )
		tag_rUVW = m.GetTag("UVW");
	else
		tag_UVW = m.GetTag("UVW");
	
	TagReal      tag_PORO = m.GetTag("PORO");
	TagReal      tag_P = m.GetTag("P");
	
	real rho_f = test->FluidDensity();
	real rho_s = test->SolidDensity();
	real gravity = test->Gravity();
	real energyU = 0; //strain energy
	
#if defined(USE_OMP)
#pragma omp parallel reduction(+:energyU)
#endif
	{
		rMatrix RR, NN, FF, x2(3,1), x1(3,1), n(3,1), r(3,1), sigma(6,1), F(3,1);
		rMatrix w(6,1), B(3,3), e(6,1), s(6,1), E(6,6), N(3,6), q(6,1);
		q(0,0) = 1;
		q(1,0) = 1;
		q(2,0) = 1;
		q(3,0) = 0;
		q(4,0) = 0;
		q(5,0) = 0;
		real pf, phi;
#if defined(USE_OMP)
#pragma omp for
#endif
		for(int i = 0; i < m.CellLastLocalID(); ++i) if( m.isValidCell(i) )
		{
			Cell c1 = m.CellByLocalID(i);
			if( c1.GetStatus() == Element::Ghost ) continue;
			ElementArray<Face> faces = c1.getFaces();
			int NF = (int)faces.size();
			RR.Resize(6,3*NF);
			NN.Resize(3*NF,6);
			FF.Resize(3*NF,1);
			c1.Barycenter(x1.data());
			E = tag_E(c1,6,6);
			RR.Zero();
			NN.Zero();
			FF.Zero();
			
			phi = tag_PORO[c1];
			
			B = tag_B(c1,3,3);
			
			w(0,0) = B(0,0);
			w(1,0) = B(1,1);
			w(2,0) = B(2,2);
			w(3,0) = B(1,2);
			w(4,0) = B(0,2);
			w(5,0) = B(0,1);
			
			for(int m = 0; m < faces.size(); ++m)
			{
				faces[m].Barycenter(x2.data());
				faces[m].OrientedNormal(c1,n.data());
				F = tag_FLUX(faces[m],4,1)*(faces[m].FaceOrientedOutside(c1)? 1.0 : -1.0)*faces[m].Area();
				if( faces[m].GetMarker(orient) ) F *= -1;
				
				r = x2-x1;
				
				RR(0,0+3*m) = r(0,0);
				RR(1,1+3*m) = r(1,0);
				RR(2,2+3*m) = r(2,0);
				
				RR(4,0+3*m) = r(2,0)*0.5;
				RR(5,0+3*m) = r(1,0)*0.5;
				
				RR(3,1+3*m) = r(2,0)*0.5;
				RR(5,1+3*m) = r(0,0)*0.5;
				
				RR(3,2+3*m) = r(1,0)*0.5;
				RR(4,2+3*m) = r(0,0)*0.5;
				
				N(0,0) = n(0,0);
				N(1,1) = n(1,0);
				N(2,2) = n(2,0);
				
				N(0,4) = n(2,0);
				N(0,5) = n(1,0);
				
				N(1,3) = n(2,0);
				N(1,5) = n(0,0);
				
				N(2,3) = n(1,0);
				N(2,4) = n(0,0);
				
				//NN(3*m,3*m+3,0,6).Print();
				
				//NQ = NN(3*m,3*m+3,0,6)*E;
				
				NN(3*m,3*m+3,0,6) = N;
				
				//NN(3*m,3*m+3,0,6).Print();
				if( test->first_order_flux )
					pf = tag_P[c1];
				else
					pf = tag_UVWPF[faces[m]][3];
				
				FF(3*m,3*m+3,0,1) = F(0,3,0,1) + N*w*pf;
				
				if( test->gravity_in_flux )
				{
					FF(3*m,3*m+3,0,1) -= N*q*(x2(2,0)*gravity*(rho_f*phi+rho_s*(1-phi)));
				}
			}
			
			//p = tag_P[c1];
			
			//FF += NN*w*p;
			
			//std::cout << "(R^T*N)^{-1}" <<std::endl;
			//(RR*NN).Invert(true).first.Print();
			//std::cout << "1.0/volume: " << 1.0/GetVolume(c1) << std::endl;
			
			//std::cout << "R:" << std::endl;
			//RR.Transpose().Print();
			//std::cout << "N:" << std::endl;
			//NN.Print();
			//std::cout << "E:" <<std::endl;
			//E.Print();
			
			
			//stress
			s = (RR*NN).Invert()*RR*FF;
			
			//strain
			e = E.Solve(s);
			
			energyU += 0.5*e.DotProduct(s)*c1.Volume();
			
			tag_S(c1,6,1) = s;
		}
	}
	energyU=m.Integrate(energyU);
	std::cout << "Strain energy: " << energyU << std::endl;
	m.ExchangeData(tag_S,CELL,0);
	
	
	std::fstream stat;
	
	if( m.GetProcessorRank() == 0 )
	{
		stat.open("stat.csv",std::ios::in);
		if( stat.fail() )
		{
			stat.clear();
			stat.open("stat.csv",std::ios::out);
			stat << "cells; faces; dT; T; ";
			stat << "u_C; v_C; w_C; p_C; uvw_C; uvwp_C; ";
			stat << "u_L2; v_L2; w_L2; p_L2; uvw_L2; uvwp_L2; ";
			stat << "uvw_mag_min; uvw_mag_max; uvw_ref_mag_min; uvw_ref_mag_max; ";
			stat << "p_min; p_max; p_ref_min; p_ref_max; ";
			stat << "Fx_C; Fy_C; Fz_C; q_C; F_C; Fq_C; ";
			stat << "Fx_L2; Fy_L2; Fz_L2; q_L2; F_L2; Fq_L2; ";
			stat << "Sxx_C; Syy_C; Szz_C; Syz_C; Sxz_C; Sxy_C; Smag_C; ";
			stat << "Sxx_L2; Syy_L2; Szz_L2; Syz_L2; Sxz_L2; Sxy_L2; Smag_L2; ";
			stat << "strain_energy; grid; problem; staggered; time-stepping;" << std::endl;
			stat.close();
		}
		else stat.close();
		
		stat.open("stat.csv",std::ios::app);
	}
	int ncells = m.TotalNumberOf(CELL);
	int nfaces = m.TotalNumberOf(FACE);
	if( m.GetProcessorRank() == 0 )
	{
		stat << ncells << ";" << nfaces << ";" << dT << ";" << T << ";";
	}
	if( test->Analytic() )
	{
		//TagRealArray tag_err  = m.CreateTag("ERR_UVWP",DATA_REAL,CELL,NONE,4);
		TagRealArray tag_erru = m.CreateTag("ERR_UVW", DATA_REAL, CELL, NONE, 4);
		TagRealArray tag_refu = m.CreateTag("REF_UVW", DATA_REAL, CELL, NONE, 4);
		//TagReal      tag_errm = m.CreateTag("ERR_UVWP_MAG",DATA_REAL,CELL,NONE,1);
		TagReal      tag_errp = m.CreateTag("ERR_P", DATA_REAL, CELL, NONE, 1);
		TagReal      tag_refp = m.CreateTag("REF_P", DATA_REAL, CELL, NONE, 1);

		real Cnorm[6] = {0,0,0,0,0,0}, L2norm[6] = {0,0,0,0,0,0}, Vnorm = 0;
		real uvw_mag_min = 1.0e20, uvw_mag_max = -1.0e20;
		real p_min = 1.0e20, p_max = -1.0e20;
		real uvw_ref_mag_min = 1.0e20, uvw_ref_mag_max = -1.0e20;
		real p_ref_min = 1.0e20, p_ref_max = -1.0e20;
		for(int i = 0; i < m.CellLastLocalID(); ++i) if( m.isValidCell(i) )
		{
			Cell c = m.CellByLocalID(i);
			if( c.GetStatus() == Element::Ghost ) continue;
			real UVWP[4], rUVWP[4], mag, cnt[3];
			c.Barycenter(cnt);
			if( test->staggered )
			{
				UVWP[0] = tag_rUVW[c][0].GetValue();
				UVWP[1] = tag_rUVW[c][1].GetValue();
				UVWP[2] = tag_rUVW[c][2].GetValue();
			}
			else
			{
				UVWP[0] = tag_UVW[c][0];
				UVWP[1] = tag_UVW[c][1];
				UVWP[2] = tag_UVW[c][2];
			}
			UVWP[3] = tag_P[c];
			
			Matrix<real,shell<real> >(shell<real>(rUVWP,4),4,1)(0,3,0,1) = test->Displacement(cnt[0],cnt[1],cnt[2],T);
			rUVWP[3] = test->Pressure(cnt[0],cnt[1],cnt[2],T).GetValue();
			tag_refu[c][0] = rUVWP[0];
			tag_refu[c][1] = rUVWP[1];
			tag_refu[c][2] = rUVWP[2];
			tag_refp[c] = rUVWP[3];
			real err[6], vol = c->Volume();
			
			mag = sqrt(UVWP[0]*UVWP[0]+UVWP[1]*UVWP[1]+UVWP[2]*UVWP[2]);
			if( uvw_mag_min > mag ) uvw_mag_min = mag;
			if( uvw_mag_max < mag ) uvw_mag_max = mag;
			mag = sqrt(rUVWP[0]*rUVWP[0]+rUVWP[1]*rUVWP[1]+rUVWP[2]*rUVWP[2]);
			if( uvw_ref_mag_min > mag ) uvw_ref_mag_min = mag;
			if( uvw_ref_mag_max < mag ) uvw_ref_mag_max = mag;
			mag = UVWP[3];
			if( p_min > mag ) p_min = mag;
			if( p_max < mag ) p_max = mag;
			mag = rUVWP[3];
			if( p_ref_min > mag ) p_ref_min = mag;
			if( p_ref_max < mag ) p_ref_max = mag;
			
			err[0] = UVWP[0] - rUVWP[0];
			err[1] = UVWP[1] - rUVWP[1];
			err[2] = UVWP[2] - rUVWP[2];
			err[3] = UVWP[3] - rUVWP[3];
			err[4] = sqrt(err[0]*err[0]+err[1]*err[1]+err[2]*err[2]);
			err[5] = sqrt(err[0]*err[0]+err[1]*err[1]+err[2]*err[2]+err[3]*err[3]);
			for(int k = 0; k < 6; ++k)
			{
				if( std::abs(err[k]) > Cnorm[k] )
					Cnorm[k] = std::abs(err[k]);
				L2norm[k] += err[k]*err[k]*vol;
			}
			Vnorm += vol;
			tag_erru[c][0] = /*tag_err[c][0] =*/ err[0];
			tag_erru[c][1] = /*tag_err[c][1] =*/ err[1];
			tag_erru[c][2] = /*tag_err[c][2] =*/ err[2];
			tag_errp[c] = /*tag_err[c][3] =*/ err[3];
			//tag_errm[c] = err[4];
		}
		Vnorm = m.Integrate(Vnorm);
		p_min = m.AggregateMin(p_min);
		p_max = m.AggregateMax(p_max);
		p_ref_min = m.AggregateMin(p_ref_min);
		p_ref_max = m.AggregateMax(p_ref_max);
		uvw_mag_min = m.AggregateMin(uvw_mag_min);
		uvw_mag_max = m.AggregateMax(uvw_mag_max);
		uvw_ref_mag_min = m.AggregateMin(uvw_ref_mag_min);
		uvw_ref_mag_max = m.AggregateMax(uvw_ref_mag_max);
		for(int k = 0; k < 6; ++k)
		{
			Cnorm[k] = m.AggregateMax(Cnorm[k]);
			L2norm[k] = m.Integrate(L2norm[k]);
			L2norm[k] = sqrt(L2norm[k]/Vnorm);
		}
		if( m.GetProcessorRank() == 0 )
		{
			std::string name[6] = {" U "," V "," W "," P ","UVWmag","mag"};
			for(int k = 0; k < 6; ++k)
			{
				std::cout << name[k] << " Cnorm " << Cnorm[k] << " L2norm " << L2norm[k] << std::endl;
			}
			std::cout << "Solution  magnetude bounds displacement: " << uvw_mag_min << ":" << uvw_mag_max << std::endl;
			std::cout << "Reference magnetude bounds displacement: " << uvw_ref_mag_min << ":" << uvw_ref_mag_max << std::endl;
			std::cout << "Solution  magnetude bounds pressure: " << p_min << ":" << p_max << std::endl;
			std::cout << "Reference magnetude bounds pressure: " << p_ref_min << ":" << p_ref_max << std::endl;
			for(int k = 0; k < 6; ++k) stat << Cnorm[k] << ";";
			for(int k = 0; k < 6; ++k) stat << L2norm[k] << "; ";
			stat << uvw_mag_min << "; " << uvw_mag_max << "; ";
			stat << uvw_ref_mag_min << "; " << uvw_ref_mag_max << "; ";
			stat << p_min << "; " << p_max << "; ";
			stat << p_ref_min << "; " << p_ref_max << "; ";
		}
	}
	else
	{
		real uvw_mag_min = 1.0e20, uvw_mag_max = -1.0e20;
		real p_min = 1.0e20, p_max = -1.0e20;
		real uvw_ref_mag_min = 1.0e20, uvw_ref_mag_max = -1.0e20;
		real p_ref_min = 1.0e20, p_ref_max = -1.0e20;
		for(int i = 0; i < m.CellLastLocalID(); ++i) if( m.isValidCell(i) )
		{
			Cell c = m.CellByLocalID(i);
			real UVWP[4], mag;
			if( test->staggered )
			{
				UVWP[0] = tag_rUVW[c][0].GetValue();
				UVWP[1] = tag_rUVW[c][1].GetValue();
				UVWP[2] = tag_rUVW[c][2].GetValue();
			}
			else
			{
				UVWP[0] = tag_UVW[c][0];
				UVWP[1] = tag_UVW[c][1];
				UVWP[2] = tag_UVW[c][2];
			}
			UVWP[3] = tag_P[c];
			
			
			mag = sqrt(UVWP[0]*UVWP[0]+UVWP[1]*UVWP[1]+UVWP[2]*UVWP[2]);
			if( uvw_mag_min > mag ) uvw_mag_min = mag;
			if( uvw_mag_max < mag ) uvw_mag_max = mag;
			
			mag = fabs(UVWP[3]);
			if( p_min > mag ) p_min = mag;
			if( p_max < mag ) p_max = mag;
		}
		p_min = m.AggregateMin(p_min);
		p_max = m.AggregateMax(p_max);
		p_ref_min = m.AggregateMin(p_ref_min);
		p_ref_max = m.AggregateMax(p_ref_max);
		uvw_mag_min = m.AggregateMin(uvw_mag_min);
		uvw_mag_max = m.AggregateMax(uvw_mag_max);
		uvw_ref_mag_min = m.AggregateMin(uvw_ref_mag_min);
		uvw_ref_mag_max = m.AggregateMax(uvw_ref_mag_max);
		if( m.GetProcessorRank() == 0 )
		{
			for(int k = 0; k < 6; ++k) stat  << "NA;";
			for(int k = 0; k < 6; ++k) stat  << "NA; ";
			stat << uvw_mag_min << "; " << uvw_mag_max << "; ";
			stat << "NA; " << "NA; ";
			stat << p_min << "; " << p_max << "; ";
			stat << "NA; " << "NA; ";
		}
	}
	
	if( test->Analytic() )
	{
		TagRealArray tag_err  = m.CreateTag("ERR_Fq",DATA_REAL,FACE,NONE,4);
		TagReal      tag_errm = m.CreateTag("ERR_Fq_MAG",DATA_REAL,FACE,NONE,1);
		real Cnorm[6] = {0,0,0,0,0,0}, L2norm[6] = {0,0,0,0,0,0}, Vnorm = 0;
		for(int i = 0; i < m.FaceLastLocalID(); ++i) if( m.isValidFace(i) )
		{
			Face c = m.FaceByLocalID(i);
			if( c.GetStatus() == Element::Ghost ) continue;
			real cnt[3],nrm[3];
			c.Barycenter(cnt);
			c.UnitNormal(nrm);
			Matrix<real,real_array> Fq(tag_FLUX[c],4,1);
			rMatrix rFq = test->Flux(cnt[0],cnt[1],cnt[2],T-dT*(1-time_shift),nrm);
			real err[6], vol = 0;
			ElementArray<Cell> fcells = c.getCells();
			for(int k = 0; k < fcells.size(); ++k)
				vol += fcells[k]->Volume()/(1.*fcells[k].nbAdjElements(FACE));
			err[0] = Fq(0,0) - rFq(0,0);
			err[1] = Fq(1,0) - rFq(1,0);
			err[2] = Fq(2,0) - rFq(2,0);
			err[3] = Fq(3,0) - rFq(3,0);
			err[4] = sqrt(err[0]*err[0]+err[1]*err[1]+err[2]*err[2]);
			err[5] = sqrt(err[0]*err[0]+err[1]*err[1]+err[2]*err[2]+err[3]*err[3]);
			for(int k = 0; k < 6; ++k)
			{
				if( std::abs(err[k]) > Cnorm[k] )
					Cnorm[k] = std::abs(err[k]);
				L2norm[k] += err[k]*err[k]*vol;
			}
			Vnorm += vol;
			tag_err[c][0] = err[0];
			tag_err[c][1] = err[1];
			tag_err[c][2] = err[2];
			tag_err[c][3] = err[3];
			tag_errm[c] = err[4];
		}
		Vnorm = m.Integrate(Vnorm);
		for(int k = 0; k < 6; ++k)
		{
			Cnorm[k] = m.AggregateMax(Cnorm[k]);
			L2norm[k] = m.Integrate(L2norm[k]);
			L2norm[k] = sqrt(L2norm[k]/Vnorm);
		}
		if( m.GetProcessorRank() == 0 )
		{
			std::string name[6] = {"Fx ","Fy ","Fz "," q ","Fmag","mag"};
			for(int k = 0; k < 6; ++k)
			{
				std::cout << name[k] << " Cnorm " << Cnorm[k] << " L2norm " << L2norm[k] << std::endl;
			}
			for(int k = 0; k < 6; ++k) stat << Cnorm[k] << ";";
			for(int k = 0; k < 6; ++k) stat << L2norm[k] << ";";
		}
	}
	else if( m.GetProcessorRank() == 0 )
	{
		for(int k = 0; k < 6; ++k) stat << "NA;";
		for(int k = 0; k < 6; ++k) stat << "NA;";
	}
	
	if( test->Analytic() )
	{
		//TagRealArray tag_ref  = m.CreateTag("rSTRESS",DATA_REAL,CELL,NONE,6);
		TagRealArray tag_err  = m.CreateTag("ERR_S",DATA_REAL,CELL,NONE,6);
		TagReal      tag_errm = m.CreateTag("ERR_S_MAG",DATA_REAL,CELL,NONE,1);
		real Cnorm[7] = {0,0,0,0,0,0,0}, L2norm[7] = {0,0,0,0,0,0,0}, Vnorm = 0;
		for(int i = 0; i < m.CellLastLocalID(); ++i) if( m.isValidCell(i) )
		{
			Cell c = m.CellByLocalID(i);
			if( c.GetStatus() == Element::Ghost ) continue;
			real cnt[3];
			c.Barycenter(cnt);
			Matrix<real,real_array> S(tag_S[c],6,1);
			rMatrix rS = test->Stress(cnt[0],cnt[1],cnt[2],T-dT*(1-time_shift));
			//tag_ref(c,6,1) = rS;
			real err[7], vol = c->Volume();
			err[0] = S(0,0) - rS(0,0);
			err[1] = S(1,0) - rS(1,0);
			err[2] = S(2,0) - rS(2,0);
			err[3] = S(3,0) - rS(3,0);
			err[4] = S(4,0) - rS(4,0);
			err[5] = S(5,0) - rS(5,0);
			err[6] = sqrt(err[0]*err[0]+err[1]*err[1]+err[2]*err[2]+err[3]*err[3]+err[4]*err[4]+err[5]*err[5]);
			for(int k = 0; k < 7; ++k)
			{
				if( std::abs(err[k]) > Cnorm[k] )
					Cnorm[k] = std::abs(err[k]);
				L2norm[k] += err[k]*err[k]*vol;
			}
			Vnorm += vol;
			tag_err[c][0] = err[0];
			tag_err[c][1] = err[1];
			tag_err[c][2] = err[2];
			tag_err[c][3] = err[3];
			tag_err[c][4] = err[4];
			tag_err[c][5] = err[5];
			tag_errm[c] = err[6];
		}
		Vnorm = m.Integrate(Vnorm);
		for(int k = 0; k < 7; ++k)
		{
			Cnorm[k] = m.AggregateMax(Cnorm[k]);
			L2norm[k] = m.Integrate(L2norm[k]);
			L2norm[k] = sqrt(L2norm[k]/Vnorm);
		}
		if( m.GetProcessorRank() == 0 )
		{
			std::string name[7] = {"Sxx ","Syy ","Szz ","Syz ","Sxz ","Sxy ","Smag"};
			for(int k = 0; k < 7; ++k)
			{
				std::cout << name[k] << " Cnorm " << Cnorm[k] << " L2norm " << L2norm[k] << std::endl;
			}
			for(int k = 0; k < 7; ++k) stat << Cnorm[k] << ";";
			for(int k = 0; k < 7; ++k) stat << L2norm[k] << ";";
		}
	}
	else if( m.GetProcessorRank() == 0 )
	{
		for(int k = 0; k < 7; ++k) stat << "NA;";
		for(int k = 0; k < 7; ++k) stat << "NA;";
	}
	
	if(m.GetProcessorRank() == 0 )
	{
		std::cout << "Strain energy: " << energyU << std::endl;
		
		stat << energyU << "; ";
		
		if( m.HaveTag("GRIDNAME") )
		{
			Tag tagname = m.GetTag("GRIDNAME");
			Storage::bulk_array name = m.self().BulkArray(tagname);
			stat << std::string(name.begin(),name.end()) << ";";
		}
		else stat << "NA;";
		stat << test->Number() << "; ";
		/*
		 if( m.HaveTag("PROBLEMNAME") )
		 {
		 Tag tagname = m.GetTag("PROBLEMNAME");
		 Storage::bulk_array name = m.self().BulkArray(tagname);
		 stat << std::string(name.begin(),name.end()) << ";";
		 }
		 else stat << "NA;";
		 */
		if( test->staggered )
			stat << "staggered;";
		else
			stat << "cell-centered;";
		if( test->bdf2 )
			stat << "BDF2;";
		else if( test->cranknicholson )
			stat << "CN;";
		else
			stat << "BE;";
		stat << std::endl;
	}
}

/// Decomposes elastic tensor into 9 co-normal 3x3 tensors and computes
/// a matrix of projections of co-normal tensors onto normal component
/// and respective tangential parts.
void SplitConormals(const real_array & E,
						   const rMatrix & n,
						   rMatrix & L, // 3x3
						   rMatrix & T, // 3x12
						   bool transversal)
{
	rMatrix C = rMatrix::FromTensor(E.data(),E.size(),6);
	rMatrix K(3,3); //co-normal matrix
	rMatrix G(3,1); //transversal direction
	// extract \mathbb{K}_1
	K(0,0) = C(0,0);
	K(1,1) = C(5,5);
	K(2,2) = C(4,4);
	K(0,1) = K(1,0) = C(0,5);
	K(0,2) = K(2,0) = C(0,4);
	K(1,2) = K(2,1) = C(4,5);
	
	L(0,0) = n.DotProduct(K*n);
	
	G = K*n;
	if( transversal ) G -= L(0,0)*n; // \mathbf{g}_1
	T(0,0+0) = G(0,0);
	T(0,1+0) = G(1,0);
	T(0,2+0) = G(2,0);
	
	// extract \mathbb{K}_2
	K(0,0) = C(5,5);
	K(1,1) = C(1,1);
	K(2,2) = C(3,3);
	K(0,1) = K(1,0) = C(1,5);
	K(0,2) = K(2,0) = C(3,5);
	K(1,2) = K(2,1) = C(1,3);
	
	L(1,1) = n.DotProduct(K*n);
	
	G = K*n;
	if( transversal ) G -= L(1,1)*n; // \mathbf{g}_2
	T(1,0+3) = G(0,0);
	T(1,1+3) = G(1,0);
	T(1,2+3) = G(2,0);
	
	// extract \mathbb{K}_3
	K(0,0) = C(4,4);
	K(1,1) = C(3,3);
	K(2,2) = C(2,2);
	K(0,1) = K(1,0) = C(3,4);
	K(0,2) = K(2,0) = C(2,4);
	K(1,2) = K(2,1) = C(2,3);
	
	L(2,2) = n.DotProduct(K*n);
	
	G = K*n;
	if( transversal ) G -= L(2,2)*n; // \mathbf{g}_3
	T(2,0+6) = G(0,0);
	T(2,1+6) = G(1,0);
	T(2,2+6) = G(2,0);
	
	// extract \mathbb{K}_4
	K(0,0) = C(0,5);
	K(0,1) = C(5,5);
	K(0,2) = C(4,5);
	K(1,0) = C(0,1);
	K(1,1) = C(1,5);
	K(1,2) = C(1,4);
	K(2,0) = C(0,3);
	K(2,1) = C(3,5);
	K(2,2) = C(3,4);
	
	L(0,1) = L(1,0) = n.DotProduct(K*n);
	
	G = K*n;
	if( transversal ) G -= L(0,1)*n; // \mathbf{g}_4
	T(0,0+3) = G(0,0);
	T(0,1+3) = G(1,0);
	T(0,2+3) = G(2,0);
	G = K.Transpose()*n;
	if( transversal ) G -= L(1,0)*n; // \mathbf{\widetilde{g}}_4
	T(1,0+0) = G(0,0);
	T(1,1+0) = G(1,0);
	T(1,2+0) = G(2,0);
	
	// extract \mathbb{K}_5
	K(0,0) = C(0,4);
	K(0,1) = C(4,5);
	K(0,2) = C(4,4);
	K(1,0) = C(0,3);
	K(1,1) = C(3,5);
	K(1,2) = C(3,4);
	K(2,0) = C(0,2);
	K(2,1) = C(2,5);
	K(2,2) = C(2,4);
	
	
	L(0,2) = L(2,0) = n.DotProduct(K*n);
	
	G = K*n;
	if( transversal ) G -= L(0,2)*n; // \mathbf{g}_5
	T(0,0+6) = G(0,0);
	T(0,1+6) = G(1,0);
	T(0,2+6) = G(2,0);
	G = K.Transpose()*n;
	if( transversal ) G -= L(2,0)*n; // \mathbf{\widetilde{g}}_5
	T(2,0+0) = G(0,0);
	T(2,1+0) = G(1,0);
	T(2,2+0) = G(2,0);
	
	// extract \mathbb{K}_6
	K(0,0) = C(4,5);
	K(0,1) = C(1,4);
	K(0,2) = C(3,4);
	K(1,0) = C(3,5);
	K(1,1) = C(1,3);
	K(1,2) = C(3,3);
	K(2,0) = C(2,5);
	K(2,1) = C(1,2);
	K(2,2) = C(2,3);
	
	L(1,2) = L(2,1) = n.DotProduct(K*n);
	
	G = K*n;
	if( transversal ) G -= L(1,2)*n; // \mathbf{g}_6
	T(1,0+6) = G(0,0);
	T(1,1+6) = G(1,0);
	T(1,2+6) = G(2,0);
	G = K.Transpose()*n;
	if( transversal ) G -= L(2,1)*n; // \mathbf{\widetilde{g}}_6
	T(2,0+3) = G(0,0);
	T(2,1+3) = G(1,0);
	T(2,2+3) = G(2,0);
}



#include "test1.h"
#include "test2.h"
#include "test3.h"
#include "test4.h"
#include "test5.h"
#include "test6.h"
#include "test7.h"
#include "test8.h"
#include "test9.h"
#include "test10.h"
#include "test11.h"
#include "test12.h"
#include "test13.h"
#include "test14.h"

AbstractTest * MakeTest(int testn)
{
	if (testn == 1)
		return new Test1();
	else if( testn == 2 )
		return new Test2();
	else if( testn == 3 )
		return new Test3();
	else if( testn == 4 )
		return new Test4();
	else if (testn == 5)
		return new Test5();
	else if( testn == 6 )
		return new Test6();
	else if( testn == 7 )
		return new Test7();
	else if (testn == 8)
		return new Test8();
	else if (testn == 9)
		return new Test9();
	else if( testn == 10 )
		return new Test10();
	else if( testn == 11 )
		return new Test11();
	else if (testn == 12)
		return new Test12();
	else if( testn == 13 )
		return new Test13();
	else if (testn == 14)
		return new Test14();
	return NULL;
}


AbstractTest::AbstractTest()
{
	bdf2 = 0;
	staggered = 0;
	gravity_in_flux = 0;
	cranknicholson = 0;
	first_order_flux = 0;
	func_output = 1;
}


void varPrint(const variable & v, double tol)
{
	if( fabs(v.GetValue()) > 0 )
		std::cout << v.GetValue() << " ";
	else
		std::cout << "0 ";
	if( !v.GetRow().Empty() )
	{
		std::stringstream s;
		bool first = true;
		for(int k = 0; k < v.GetRow().Size(); ++k) if( fabs(v.GetRow().GetValue(k)) > tol )
		{
			if( first ) first = false;
			else s << ",";
			s << "(" << v.GetRow().GetIndex(k) << "," << v.GetRow().GetValue(k) << ")";
		}
		if( !s.str().empty() )
			std::cout << "{" << s.str() << "}";
	}
}

void vMatrixPrint(const vMatrix & m, double tol)
{
	for(int i = 0; i < m.Rows(); ++i)
	{
		for(int j = 0; j < m.Cols(); ++j)
		{
			std::cout << "[" << i << "," << j << "]: ";
			varPrint(m(i,j),tol);
			std::cout << std::endl;
		}
	}
}


rMatrix GenIsotropicTensor(double E, double nu)
{
	rMatrix C(6,6,0.0);
	double lambda = E * nu / (1 + nu) / (1 - 2 * nu);
	double G = E/(2*(1+nu));
	C(0,0) = C(1,1) = C(2,2) = lambda + 2*G;
	C(0,1) = C(1,0) = lambda;
	C(0,2) = C(2,0) = lambda;
	C(1,2) = C(2,1) = lambda;
	C(3,3) = C(4,4) = C(5,5) = G;
	return C;
}

rMatrix GenAnisotropicTensor(double E1, double E2, double E3, double nu12, double nu13, double nu23, double G12, double G13, double G23)
{
	rMatrix C(6,6,0.0);
	double nu21 = E2/E1*nu12;
	double nu31 = E3/E1*nu13;
	double nu32 = E3/E2*nu23;
	
	C(0,0) =   1.0/E1;	C(0,1) = -nu21/E2;	C(0,2) = -nu31/E3;
	C(1,0) = -nu12/E1;	C(1,1) =   1.0/E2;	C(1,2) = -nu32/E3;
	C(2,0) = -nu13/E1;	C(2,1) = -nu23/E2;	C(2,2) =   1.0/E3;
	
	C(3,3) = 1.0/(2.0*G23);
	C(4,4) = 1.0/(2.0*G13);
	C(5,5) = 1.0/(2.0*G12);
	
	return C.Invert();
}
