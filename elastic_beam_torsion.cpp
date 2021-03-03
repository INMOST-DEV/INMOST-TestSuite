#include "inmost.h"

using namespace INMOST;

bool print = false;
const Storage::real pi = 3.1415926535897932384626433832795;
const Storage::real eps = 1.0e-6;
std::string problem_name = "elastic_beam_shear";


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
	
	double L = 10;
	double max[3] = {-1.0e20, -1.0e20, -1.0e20}, min[3] = {1.0e20,1.0e20,1.0e20};
	Storage::real c[3] = {0,0,0}, nrm[3] = {0,0,0};
	for(Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); ++it)
	{
		it->Centroid(c);
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
	
	//rescale mesh to the unit cube
	for(Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); ++it)
	{
		Storage::real_array c = it->Coords();
		c[0] = (c[0]-min[0])/(max[0]-min[0]);
		c[1] = (c[1]-min[1])/(max[1]-min[1]);
		c[2] = (c[2]-min[2])/(max[2]-min[2]);
	}
	
	
	//shift mesh be translation from (0.5,0.5,0) to (0,0,0)
	for(Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); ++it)
	{
		Storage::real_array c = it->Coords();
		c[0] -= 0.5;
		c[1] -= 0.5;
	}
	//rescale mesh into [-1,1]x[-1,1]x[0,10]
	for(Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); ++it)
	{
		Storage::real_array c = it->Coords();
		c[0] *= 2.0;
		c[1] *= 2.0;
		c[2] *= L;
	}
	
	
	Mesh::GeomParam t;
	t[MEASURE] = FACE |CELL;
	t[ORIENTATION] = FACE;
	t[NORMAL] = FACE;
	t[CENTROID] = CELL | FACE | EDGE | NODE;
	m->RemoveGeometricData(t);
	m->PrepareGeometricData(t);
	
	
	
	
	
	
	
	if( m->HaveTag("ELASTIC_TENSOR") ) m->DeleteTag(m->GetTag("ELASTIC_TENSOR"));
	//if( m->HaveTag("MATERIAL") ) m->DeleteTag(m->GetTag("MATERIAL"));
	
	
	
	Tag tensor = m->CreateTag("ELASTIC_TENSOR",DATA_REAL,CELL,NONE,21);
	Tag solution_val = m->CreateTag("REFERENCE_SOLUTION",DATA_REAL,CELL,NONE,3);
	Tag stress_val = m->CreateTag("REFERENCE_STRESS",DATA_REAL,CELL,NONE,6);
	Tag solution_flux = m->CreateTag("REFERENCE_FLUX",DATA_REAL,FACE,NONE,3);
	Tag bndcond = m->CreateTag("BOUNDARY_CONDITION",DATA_REAL,FACE,FACE,6);
	//Tag force = m->CreateTag("FORCE",DATA_REAL,CELL,NONE,3);
	//Tag grad_val = m->CreateTag("REFERENCE_GRADIENT",DATA_REAL,CELL,NONE,9);

	{
		Storage::bulk_array name = m->self()->BulkArray(m->CreateTag("PROBLEMNAME",DATA_BULK,MESH,NONE));
		name.replace(name.begin(),name.end(),problem_name.begin(),problem_name.end());
	}
	
	rMatrix C(6,6); //define elastic tensor with lame coefficients
	
	//Lame constants, young modulus E = 100000
	double E = 25;
	double nu = 0.3;
	double lambda = E*nu/((1+nu)*(1-2*nu));
	double mu = E/(2*(1+nu));
	double beta = 0.1;
	C.Zero();
	C(0,0) = C(1,1) = C(2,2) = lambda + 2*mu;
	C(0,1) = C(1,0) = lambda;
	C(0,2) = C(2,0) = lambda;
	C(1,2) = C(2,1) = lambda;
	C(3,3) = C(4,4) = C(5,5) = mu;
	
	//make it pseudo-2d
	//C(0,2) = C(2,0) = C(1,2) = C(2,1) = 0;
	//C(2,2) = C(3,3) = C(4,4) = 1;

	//0: u_x
	//1: v_y
	//2: w_z
	//3: v_z + w_y
	//4: u_z + w_x
	//5: u_y + v_x
	
	C.Print();
	
	for(Mesh::iteratorElement it = m->BeginElement(CELL|FACE); it != m->EndElement(); ++it)
	{
		it->Centroid(c);
		
		double x = c[0];
		double y = c[1];
		double z = c[2];
		
		//std::cout << "radius: " << get_value(u_r) << " angle: " << get_value(u_phi) << std::endl;
		
		rMatrix sol(3,1);
		rMatrix sigma(6,1);
		rMatrix N(3,6), flux(3,1);
		
		sol(0,0) = -beta*y*z;
		sol(1,0) = beta*x*z;
		sol(2,0) = beta*x*y;
		
		
		
		sigma(0,0) = sigma(1,0) = sigma(2,0) = sigma(5,0) = 0; //xx,yy,zz,xy
		sigma(3,0) = E*beta/(1+nu)*x;//yz = v_z + w_y
		sigma(4,0) = 0;//xz
		
		double minus1 = -1;
		for(int i = 1; i < 100; ++i)
		{
			int n = i;
			sol(2,0) += 32*beta/(pi*pi*pi)*minus1/((2*n-1)*(2*n-1)*(2*n-1)*cosh((2*n-1)*pi/2))*sin((2*n-1)*pi*x/2)*sinh((2*n-1)*pi*y/2);
			sigma(3,0) += 8*E*beta/(pi*pi*(1+nu))*minus1/((2*n-1)*(2*n-1)*cosh((2*n-1)*pi/2))*sin((2*n-1)*pi*x/2)*cosh((2*n-1)*pi*y/2);
			sigma(4,0) += 8*E*beta/(pi*pi*(1+nu))*minus1/((2*n-1)*(2*n-1)*cosh((2*n-1)*pi/2))*cos((2*n-1)*pi*x/2)*sinh((2*n-1)*pi*y/2);
			minus1 *= -1;
		}
		
		if( it->GetElementType() == CELL )
		{
			//it->RealArray(force)[0] = beta;
			//it->RealArray(force)[1] = -beta;
			//tensor layout
			//C_1  C_7  C_8  C_9  C_10 C_11
			//     C_2  C_12 C_13 C_14 C_15
			//          C_3  C_16 C_17 C_18
			//               C_4  C_19 C_20
			//                    C_5  C_21
			//                         C_6
			Storage::real_array perm = it->RealArray(tensor);
			perm[ 0] = C(0,0);
			perm[ 1] = C(0,1);
			perm[ 2] = C(0,2);
			perm[ 3] = C(0,3);
			perm[ 4] = C(0,4);
			perm[ 5] = C(0,5);
			perm[ 6] = C(1,1);
			perm[ 7] = C(1,2);
			perm[ 8] = C(1,3);
			perm[ 9] = C(1,4);
			perm[10] = C(1,5);
			perm[11] = C(2,2);
			perm[12] = C(2,3);
			perm[13] = C(2,4);
			perm[14] = C(2,5);
			perm[15] = C(3,3);
			perm[16] = C(3,4);
			perm[17] = C(3,5);
			perm[18] = C(4,4);
			perm[19] = C(4,5);
			perm[20] = C(5,5);
			
	 
			Storage::real_array rsol = it->RealArray(solution_val);
			rsol[0] = sol(0,0);
			rsol[1] = sol(1,0);
			rsol[2] = sol(2,0);
			
			
			Storage::real_array s = it->RealArray(stress_val);
			s[0] = sigma(0,0);
			s[1] = sigma(1,0);
			s[2] = sigma(2,0);
			s[3] = sigma(3,0);
			s[4] = sigma(4,0);
			s[5] = sigma(5,0);
			
		}
		
		
		
		else if( it->GetElementType() == FACE )
		{
			it->getAsFace()->UnitNormal(nrm);
			
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
			
			flux = N*sigma;
			
			
			Storage::real_array rflux = it->RealArray(solution_flux);
			
			rflux[0] = flux(0,0);
			rflux[1] = flux(1,0);
			rflux[2] = flux(2,0);
				
			if( it->getAsFace()->Boundary() )
			{
				
				// alpha*u + beta*(I-p*n*n^T)*N*sigma = gamma
				//if(  c[2] > 10 - eps ) //dirichlet
				if( c[2] < eps )
				{
					Storage::real_array bc = it->RealArray(bndcond);
					bc[0] = 1.0; //alpha
					bc[1] = 0.0; //beta
					bc[2] = 0.0; //projection, gamma
					//right hand side
					bc[3] = sol(0,0);
					bc[4] = sol(1,0);
					bc[5] = sol(2,0);
				}
				else //if( c[2] < eps )
					if(  c[2] > 10 - eps )
				{
					Storage::real_array bc = it->RealArray(bndcond);
					bc[0] = 0.0; //alpha
					bc[1] = 1.0; //beta
					bc[2] = 0.0; //projection, gamma
					//right hand side
					bc[3] = flux(0,0);
					bc[4] = flux(1,0);
					bc[5] = flux(2,0);
				}
			}
			
		}
	}
	
	
	std::cout << "Saving output to " << argv[2] << std::endl;
	
	m->Save(argv[2]);
	
	delete m;
}


