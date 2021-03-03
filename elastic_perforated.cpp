#include "inmost.h"

using namespace INMOST;

bool print = false;
const Storage::real pi = 3.1415926535897932384626433832795;
const Storage::real eps = 1.0e-6;
std::string problem_name = "elastic_perforated_strip";
//defined from
// https://www.math.hu-berlin.de/~cc/cc_homepage/download/2002-AJ_CC_FS_KR-Matlab_Implementation_FEM_Elasticity.pdf


template<int k>
double fsigma(double * x, double)
{
	const double fx = 10000;
	const double a = 0.5;
	double r = sqrt(x[0]*x[0]+x[1]*x[1]);
	double phi = atan2(x[1],x[0]);
	if( k == 0 )
		return fx*(1 - pow(a/r,2)*(1.5*cos(2*phi)+cos(4*phi))+1.5*pow(a/r,4)*cos(4*phi));
	else if( k == 1 )
		return fx*( - pow(a/r,2)*(0.5*cos(2*phi)-cos(4*phi))-1.5*pow(a/r,4)*cos(4*phi));
	else if (k == 2 )
		return fx*( - pow(a/r,2)*(0.5*sin(2*phi)+sin(4*phi))+1.5*pow(a/r,4)*sin(4*phi)); //xy
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
	
	//rescale mesh by 2 in x and y
	for(Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); ++it)
	{
		Storage::real_array c = it->Coords();
		c[0] *= 2;
		c[1] *= 2;
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
	Tag stress_val = m->CreateTag("REFERENCE_STRESS",DATA_REAL,CELL,NONE,6);
	Tag solution_flux = m->CreateTag("REFERENCE_FLUX",DATA_REAL,FACE,NONE,3);
	Tag bndcond = m->CreateTag("BOUNDARY_CONDITION",DATA_REAL,FACE,FACE,6);
	Tag force = m->CreateTag("FORCE",DATA_REAL,CELL,NONE,3);
	

	{
		Storage::bulk_array name = m->self()->BulkArray(m->CreateTag("PROBLEMNAME",DATA_BULK,MESH,NONE));
		name.replace(name.begin(),name.end(),problem_name.begin(),problem_name.end());
	}
	
	rMatrix C(6,6); //define elastic tensor with lame coefficients
	
	//Lame constants, young modulus E = 100000
	double E = 1E+7;
	double nu = 0.3;
	double lambda = E*nu/((1+nu)*(1-2*nu));
	double mu = E/(2*(1+nu));
	double fx = 10000;
	double a = 0.5;
	double b = 2;
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
	
	int bc1 = 0, bc2 = 0, bc3 = 0;
	
	for(Mesh::iteratorElement it = m->BeginElement(CELL|FACE); it != m->EndElement(); ++it)
	{
		it->Centroid(c);
		
		double x = c[0];
		double y = c[1];
		double z = c[2];
		double r = sqrt(x*x+y*y);
		double phi = atan2(y,x);
		
		rMatrix sigma(6,1);
		rMatrix N(3,6), flux(3,1);
		
		/*
		sigma(0,0) = it->Mean(fsigma<0>,0);//fx*(1 - pow(a/r,2)*(1.5*cos(2*phi)+cos(4*phi))+1.5*pow(a/r,4)*cos(4*phi)); //xx
		sigma(1,0) = it->Mean(fsigma<1>,0);//fx*( - pow(a/r,2)*(0.5*cos(2*phi)-cos(4*phi))-1.5*pow(a/r,4)*cos(4*phi)); //yy
		sigma(2,0) = nu*(sigma(0,0)+sigma(1,0)); //zz is not zero
		sigma(3,0) = 0; //yz
		sigma(4,0) = 0; //xz
		sigma(5,0) = it->Mean(fsigma<2>,0); //fx*( - pow(a/r,2)*(0.5*sin(2*phi)+sin(4*phi))+1.5*pow(a/r,4)*sin(4*phi)); //xy
		*/
		
		sigma(0,0) = fx*(1 - pow(a/r,2)*(1.5*cos(2*phi)+cos(4*phi))+1.5*pow(a/r,4)*cos(4*phi)); //xx
		sigma(1,0) = fx*( - pow(a/r,2)*(0.5*cos(2*phi)-cos(4*phi))-1.5*pow(a/r,4)*cos(4*phi)); //yy
		sigma(2,0) = nu*(sigma(0,0)+sigma(1,0)); //zz is not zero
		sigma(3,0) = 0; //yz
		sigma(4,0) = 0; //xz
		sigma(5,0) = fx*( - pow(a/r,2)*(0.5*sin(2*phi)+sin(4*phi))+1.5*pow(a/r,4)*sin(4*phi)); //xy
		
		if( it->GetElementType() == CELL )
		{
			//it->RealArray(force)[0] = fx;
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
				Storage::real_array bc = it->RealArray(bndcond);
				// alpha*u + beta*(I-p*n*n^T)*N*sigma = gamma
				if( c[2] > 1.0-eps || c[2] < eps  ) //pure roller on z
				{
					bc[0] = 0.0; //alpha
					bc[1] = 1.0; //beta
					bc[2] = 1.0; //projection, gamma
					//right hand side
					bc[3] = 0.0;
					bc[4] = 0.0;
					bc[5] = 0.0;
					bc1++;
				}
				else if( c[0] < eps || c[1] < eps ) //symmetry on left and bottom
				{
					bc[0] = 0.0; //alpha
					bc[1] = 1.0; //beta
					bc[2] = 1.0; //projection, gamma
					//right hand side
					bc[3] = 0.0;
					bc[4] = 0.0;
					bc[5] = 0.0;
					bc1++;
				}
				else if( c[0] > 2-eps || c[1] > 2-eps ) //neumann on right and top
				{
					bc[0] = 0.0; //alpha
					bc[1] = 1.0; //beta
					bc[2] = 0.0; //projection, gamma
					//right hand side
					bc[3] = flux(0,0);
					bc[4] = flux(1,0);
					bc[5] = flux(2,0);
					bc2++;
				}
				else //zero traction on perforation
				{
					bc[0] = 0.0; //alpha
					bc[1] = 1.0; //beta
					bc[2] = 0.0; //projection, gamma
					//right hand side
					bc[3] = 0.0;
					bc[4] = 0.0;
					bc[5] = 0.0;
					bc3++;
				}
			}
			
		}
	}
	
	std::cout << "BC, pure roller " << bc1 << " traction control " << bc2 << " traction-free " << bc3 << std::endl;
	
	
	std::cout << "Saving output to " << argv[2] << std::endl;
	
	m->Save(argv[2]);
	
	delete m;
}


