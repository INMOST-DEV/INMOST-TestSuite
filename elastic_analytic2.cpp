#include "inmost.h"

using namespace INMOST;

bool print = false;
//const Storage::real pi = 3.1415926535897932384626433832795;
const Storage::real eps = 1.0e-6;
std::string problem_name = "elastic_analytic_test_2";

void SVD2Eigen(const rMatrix & U, rMatrix & S, rMatrix & V)
{
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
	//check
	if( (U-V).FrobeniusNorm() > 1.0e-8 )
		(U-V).Print();
}

int main(int argc, char ** argv)
{
	if( argc < 3 )
	{
		std::cout << "Usage: " << argv[0] << " mesh mesh_out" << std::endl;
		return 0;
	}
	
	Mesh * m = new Mesh;
	m->SetFileOption("VERBOSITY","2");
	try{m->Load(argv[1]);} catch(...) { std::cout << "Cannot load the mesh " << argv[1] << std::endl; return -1;}
	
	
	
	Mesh::GeomParam t;
	t[MEASURE] = FACE;
	t[ORIENTATION] = FACE;
	t[NORMAL] = FACE;
	t[BARYCENTER] = CELL | FACE | EDGE | NODE;
	m->PrepareGeometricData(t);
	
	
	double max[3] = {-1.0e20, -1.0e20, -1.0e20}, min[3] = {1.0e20,1.0e20,1.0e20};
	Storage::real c[3] = {0,0,0}, nrm[3] = {0,0,0};
	for(Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); ++it)
	{
		it->Barycenter(c);
		if( c[0] > max[0] ) max[0] = c[0];
		if( c[1] > max[1] ) max[1] = c[1];
		if( c[2] > max[2] ) max[2] = c[2];
		if( c[0] < min[0] ) min[0] = c[0];
		if( c[1] < min[1] ) min[1] = c[1];
		if( c[2] < min[2] ) min[2] = c[2];
	}
	
	if( max[0] <= min[0] ) {std::cout << "strange X " << min[0] << ":" << max[0] << std::endl; return -1;}
	if( max[1] <= min[1] ) {std::cout << "strange Y " << min[1] << ":" << max[1] << std::endl; return -1;}
	if( max[2] <= min[2] )
	{
		//2d mesh
		if( m->GetDimensions() == 3 )
		{
			//offset from z-plane
			min[2] -= 0.0001;
			max[2] += 0.0001;
		}
		else
		{
			min[2] = -0.0001;
			max[2] = +0.0001;
		}
	}
	
	std::cout << "Mesh bounds: " << min[0] << ":" << max[0] << " " << min[1] << ":" << max[1] << " " << min[2] << ":" << max[2] << std::endl;
	
	if( m->HaveTag("ELASTIC_TENSOR") ) m->DeleteTag(m->GetTag("ELASTIC_TENSOR"));
	if( m->HaveTag("MATERIAL") ) m->DeleteTag(m->GetTag("MATERIAL"));
	
	
	
	Tag tensor = m->CreateTag("ELASTIC_TENSOR",DATA_REAL,CELL,NONE,21);
	Tag solution_val = m->CreateTag("REFERENCE_SOLUTION",DATA_REAL,CELL|FACE,NONE,3);
	Tag stress_val = m->CreateTag("REFERENCE_STRESS",DATA_REAL,CELL,NONE,6);
	Tag solution_flux = m->CreateTag("REFERENCE_FLUX",DATA_REAL,FACE,NONE,3);
	Tag bndcond = m->CreateTag("BOUNDARY_CONDITION",DATA_REAL,FACE,FACE,6);
	Tag force = m->CreateTag("FORCE",DATA_REAL,CELL,NONE,3);
	Tag grad_val = m->CreateTag("REFERENCE_GRADIENT",DATA_REAL,CELL,NONE,9);

	{
		Storage::bulk_array name = m->self()->BulkArray(m->CreateTag("PROBLEMNAME",DATA_BULK,MESH,NONE));
		name.replace(name.begin(),name.end(),problem_name.begin(),problem_name.end());
	}
	
	
	for(Mesh::iteratorElement it = m->BeginElement(CELL|FACE); it != m->EndElement(); ++it)
	{
		it->Barycenter(c);
		unknown x(c[0],0);
		unknown y(c[1],1);
		unknown z(c[2],2);
		hMatrix sol(3,1);
		vMatrix epsilon(6,1), sigma(6,1);
		rMatrix N(3,6), sigma_value(6,1), flux(3,1);
		vMatrix C(6,6);
		rMatrix Cval(6,6), U(6,6), S(6,6), V(6,6);
		
		
		
		
		//sol(0,0) = (sin(((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) + (z-0.5)*(z-0.5))*16)+1)*(2*x-1);
		//sol(1,0) = (sin(((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) + (z-0.5)*(z-0.5))*16)+1)*(2*y-1);
		//sol(2,0) = (sin(((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) + (z-0.5)*(z-0.5))*16)+1)*(2*z-1);
		sol(0,0) = sin(((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) + (z-0.5)*(z-0.5))*8)+1;
		sol(1,0) = sin(((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) + (z-0.5)*(z-0.5))*8)+1;
		sol(2,0) = sin(((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) + (z-0.5)*(z-0.5))*8)+1;
		
		Storage::real a = 0.01, b = 0.1, c = 1, d = 0.002;
		
		C(0,0) = a*x*x + b*y*y + c*z*z+d;
		C(1,1) = c*x*x + a*y*y + b*z*z+d;
		C(2,2) = b*x*x + c*y*y + a*z*z+d;
		C(3,3) = c*x*x + b*y*y + a*z*z+d;
		C(4,4) = a*x*x + c*y*y + b*z*z+d;
		C(5,5) = b*x*x + a*y*y + c*z*z+d;
		C(0,1) = C(1,0) = -a*x*y;
		C(0,2) = C(2,0) = -b*x*z;
		C(0,3) = C(3,0) = -a*y*z;
		C(0,4) = C(4,0) = -a*(x*y+y*z);
		C(0,5) = C(5,0) = -a*(x*z+y*z);
		C(1,2) = C(2,1) = -b*(x*y+x*z);
		C(1,3) = C(3,1) = -a*(x*z);
		C(1,4) = C(4,1) = -a*(y*z);
		C(1,5) = C(5,1) = -b*(x*y);
		C(2,3) = C(3,2) = -a*(x*y+x*z);
		C(2,4) = C(4,2) = -a*(y*z+x*y);
		C(2,5) = C(5,2) = -b*(x*z+y*z);
		C(3,4) = C(4,3) = -a*(x*y*z);
		C(3,5) = C(5,3) = -a*(x*z-y*z);
		C(4,5) = C(5,4) = -a*(x*y*z);
		
		for(int i = 0; i < 6; ++i)
			for(int j = 0; j < 6; ++j)
				Cval(i,j) = get_value(C(i,j));
		
		Cval.SVD(U,S,V);
		SVD2Eigen(U,S,V);
		
		bool negative = false;
		
		for(int i = 0; i < 6; ++i) if( S(i,i) < 0.0 ) negative = true;
		
		if( negative )
		{
			for(int i = 0; i < 6; ++i)
				std::cout << S(i,i) << " ";
			std::cout << std::endl;
		}
		
		epsilon(0,0) = sol(0,0).GetVariable(0); //u_x
		epsilon(1,0) = sol(1,0).GetVariable(1); //v_y
		epsilon(2,0) = sol(2,0).GetVariable(2); //w_z
		epsilon(3,0) = sol(1,0).GetVariable(2) + sol(2,0).GetVariable(1); //v_z + w_y
		epsilon(4,0) = sol(0,0).GetVariable(2) + sol(2,0).GetVariable(0); //u_z + w_x
		epsilon(5,0) = sol(0,0).GetVariable(1) + sol(1,0).GetVariable(0); //u_y + v_x
		
		
		sigma = C*epsilon;
		
		
		
		
		if( it->GetElementType() == CELL )
		{
			//tensor layout
			//C_1  C_7  C_8  C_9  C_10 C_11
			//     C_2  C_12 C_13 C_14 C_15
			//          C_3  C_16 C_17 C_18
			//               C_4  C_19 C_20
			//                    C_5  C_21
			//                         C_6
			Storage::real_array perm = it->RealArray(tensor);
			perm[ 0] = Cval(0,0);
			perm[ 1] = Cval(0,1);
			perm[ 2] = Cval(0,2);
			perm[ 3] = Cval(0,3);
			perm[ 4] = Cval(0,4);
			perm[ 5] = Cval(0,5);
			perm[ 6] = Cval(1,1);
			perm[ 7] = Cval(1,2);
			perm[ 8] = Cval(1,3);
			perm[ 9] = Cval(1,4);
			perm[10] = Cval(1,5);
			perm[11] = Cval(2,2);
			perm[12] = Cval(2,3);
			perm[13] = Cval(2,4);
			perm[14] = Cval(2,5);
			perm[15] = Cval(3,3);
			perm[16] = Cval(3,4);
			perm[17] = Cval(3,5);
			perm[18] = Cval(4,4);
			perm[19] = Cval(4,5);
			perm[20] = Cval(5,5);
			
			rMatrix Cr = rMatrix::FromTensor(perm.data(),21,6);
			if( (C-Cr).FrobeniusNorm() > 1.0e-5 )
			{
				std::cout << "original matrix:" << std::endl;
				C.Print();
				std::cout << "stored matrix:" << std::endl;
				Cr.Print();
			}
			
			
			//sigma layot
			// s0 s5 s4
			// s5 s1 s3
			// s4 s3 s2
			// divergence is taken row-wise
			Storage::real_array f = it->RealArray(force);
			f[0] -= sigma(0,0).GetRow()[0];
			f[0] -= sigma(5,0).GetRow()[1];
			f[0] -= sigma(4,0).GetRow()[2];
			
			
			f[1] -= sigma(5,0).GetRow()[0];
			f[1] -= sigma(1,0).GetRow()[1];
			f[1] -= sigma(3,0).GetRow()[2];
			
			
			f[2] -= sigma(4,0).GetRow()[0];
			f[2] -= sigma(3,0).GetRow()[1];
			f[2] -= sigma(2,0).GetRow()[2];
			
			Storage::real_array rsol = it->RealArray(solution_val);
			rsol[0] = sol(0,0).GetValue();
			rsol[1] = sol(1,0).GetValue();
			rsol[2] = sol(2,0).GetValue();
			
			
			Storage::real_array s = it->RealArray(stress_val);
			s[0] = sigma(0,0).GetValue();
			s[1] = sigma(1,0).GetValue();
			s[2] = sigma(2,0).GetValue();
			s[3] = sigma(3,0).GetValue();
			s[4] = sigma(4,0).GetValue();
			s[5] = sigma(5,0).GetValue();
			
			Storage::real_array g = it->RealArray(grad_val);
			g[0] = sol(0,0).GetRow()[0];
			g[1] = sol(0,0).GetRow()[1];
			g[2] = sol(0,0).GetRow()[2];
			g[3] = sol(1,0).GetRow()[0];
			g[4] = sol(1,0).GetRow()[1];
			g[5] = sol(1,0).GetRow()[2];
			g[6] = sol(2,0).GetRow()[0];
			g[7] = sol(2,0).GetRow()[1];
			g[8] = sol(2,0).GetRow()[2];
		}
		
		
		
		else if( it->GetElementType() == FACE )
		{
			Storage::real_array rsol = it->RealArray(solution_val);
			rsol[0] = sol(0,0).GetValue();
			rsol[1] = sol(1,0).GetValue();
			rsol[2] = sol(2,0).GetValue();
			
			it->getAsFace()->UnitNormal(nrm);
			
			sigma_value(0,0) = sigma(0,0).GetValue();
			sigma_value(1,0) = sigma(1,0).GetValue();
			sigma_value(2,0) = sigma(2,0).GetValue();
			sigma_value(3,0) = sigma(3,0).GetValue();
			sigma_value(4,0) = sigma(4,0).GetValue();
			sigma_value(5,0) = sigma(5,0).GetValue();
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
			
			flux = N*sigma_value;
			
			
			Storage::real_array rflux = it->RealArray(solution_flux);
			
			rflux[0] = flux(0,0);
			rflux[1] = flux(1,0);
			rflux[2] = flux(2,0);
				
			if( it->getAsFace()->Boundary() )
			{
				Storage::real_array bc = it->RealArray(bndcond);
				// alpha*u + beta*(I-p*n*n^T)*N*sigma = gamma
				bc[0] = 1.0; //alpha
				bc[1] = 0.0; //beta
				bc[2] = 0.0; //projection
				//gamma
				bc[3] = sol(0,0).GetValue();
				bc[4] = sol(1,0).GetValue();
				bc[5] = sol(2,0).GetValue();
			}
			
		}
	}
	
	std::cout << "Saving output to " << argv[2] << std::endl;
	
	m->Save(argv[2]);
	
	delete m;
}


