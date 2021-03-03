#include "inmost.h"

using namespace INMOST;

bool print = false;
const Storage::real pi = 3.1415926535897932384626433832795;
const Storage::real eps = 1.0e-6;
std::string problem_name = "elastic_membrane";
//defined from
// http://www.isima.fr/~jkoko/Codes/Elas2d.pdf


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
	
	//rescale mesh to the unit cube
	for(Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); ++it)
	{
		Storage::real_array c = it->Coords();
		c[0] = (c[0]-min[0])/(max[0]-min[0]);
		c[1] = (c[1]-min[1])/(max[1]-min[1]);
		c[2] = (c[2]-min[2])/(max[2]-min[2]);
	}
	
	Mesh::GeomParam t;
	//t[MEASURE] = FACE |CELL;
	t[ORIENTATION] = FACE;
	t[NORMAL] = FACE;
	t[CENTROID] = CELL | FACE;
	m->RemoveGeometricData(t);
	
	
	if( m->HaveTag("ELASTIC_TENSOR") ) m->DeleteTag(m->GetTag("ELASTIC_TENSOR"));
	//if( m->HaveTag("MATERIAL") ) m->DeleteTag(m->GetTag("MATERIAL"));
	
	
	
	Tag tensor = m->CreateTag("ELASTIC_TENSOR",DATA_REAL,CELL,NONE,21);
	Tag bndcond = m->CreateTag("BOUNDARY_CONDITION",DATA_REAL,FACE,FACE,6);
	Tag force = m->CreateTag("FORCE",DATA_REAL,CELL,NONE,3);
	Tag contact = m->CreateTag("CONTACT_TYPE",DATA_INTEGER,FACE,FACE,1);

	{
		Storage::bulk_array name = m->self()->BulkArray(m->CreateTag("PROBLEMNAME",DATA_BULK,MESH,NONE));
		name.replace(name.begin(),name.end(),problem_name.begin(),problem_name.end());
	}
	
	rMatrix C0(6,6), C1(6,6), C(6,6); //define elastic tensor with lame coefficients
	
	//Lame constants, young modulus E = 100000
	double E = 3000;
	double nu = 0.4;
	double lambda = E*nu/((1+nu)*(1-2*nu));
	double mu = E/(2*(1+nu));
	
	C0.Zero();
	C0(0,0) = C0(1,1) = C0(2,2) = lambda + 2*mu;
	C0(0,1) = C0(1,0) = lambda;
	C0(0,2) = C0(2,0) = lambda;
	C0(1,2) = C0(2,1) = lambda;
	C0(3,3) = C0(4,4) = C0(5,5) = mu;
	
	//make it pseudo-2d
	//C(0,2) = C(2,0) = C(1,2) = C(2,1) = 0;
	//C(2,2) = C(3,3) = C(4,4) = 1;

	//0: u_x
	//1: v_y
	//2: w_z
	//3: v_z + w_y
	//4: u_z + w_x
	//5: u_y + v_x
	
	double Cv[21] =
	{
		1.323 , 0.0726 , 0.263 , 0.108 , -0.08 , -0.239,
		1.276 , -0.318 , 0.383 , 0.108 , 0.501 ,
		0.943 , -0.183 , 0.146 , 0.182 ,
		1.517 , -0.0127 , -0.304 ,
		1.209 , -0.326 ,
		1.373
	};
	C1 = rMatrix::FromTensor(Cv,21,6);
	C = C0 + C1*100.0;
	
	C.Print();
	
	for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
	{
		it->Centroid(c);
		
		C = C0 + C1*(100.0*c[0]);
		
		Storage::real_array f = it->RealArray(force);
		f[0] = 0;
		f[1] = 0;//-1.5;
		f[2] = 0;//-1.5;
		
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
	}
	int fcontact = 0;
	
	for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it)
	{
		it->Centroid(c);
		
				
		if( it->Boundary() )
		{
			Storage::real_array bc = it->RealArray(bndcond);
			// alpha*u + beta*(I-p*n*n^T)*N*sigma = gamma
			//if( false )
			if( false ) //c[2] > 1.0-eps || c[2] < 0.0+eps ) //pure roller on z
			{
				bc[0] = 0.0; //alpha
				bc[1] = 1.0; //beta
				bc[2] = 1.0; //projection, gamma
				//right hand side
				bc[3] = 0.0;
				bc[4] = 0.0;
				bc[5] = 0.0;
			}
			else if( c[0] < 0.0+eps ) //dirichlet on x = 0
			{
				bc[0] = 1.0; //alpha
				bc[1] = 0.0; //beta
				bc[2] = 0.0; //projection, gamma
				//right hand side
				bc[3] = 0;
				bc[4] = 0;
				bc[5] = 0;
			}
			else if( c[0] > 1.0-eps ) //shear load on x = 35
			{
				bc[0] = 0.0; //alpha
				bc[1] = 1.0; //beta
				bc[2] = 0.0; //projection, gamma
				//right hand side
				bc[3] = 0;
				bc[4] = 20;
				bc[5] = 20;
			}
			else //pure Neumann
			{
				bc[0] = 0.0; //alpha
				bc[1] = 1.0; //beta
				bc[2] = 0.0; //projection, gamma
				//right hand side
				bc[3] = 0;
				bc[4] = 0;
				bc[5] = 0;
			}
		}
		else
		{
			if( c[1] > 0.5-eps && c[1] < 0.5+eps &&
			    c[0] > 0.2 && c[0] < 0.8 &&
			    c[2] > 0.2 && c[2] < 0.8)
			{
				it->Integer(contact) = 2;
				fcontact++;
			}
		}
	}
	
	std::cout << "faces with contact: " << fcontact << std::endl;
	
	
	//rescale mesh to fit quad (0, −5), (35, −2), (35, 2) (0, 5) in Oxy
	for(Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); ++it)
	{
		Storage::real_array c = it->Coords();
		double a = c[0], b = c[1], q = c[2];
		c[0] = c[0]*35; // x is from 0 to 35 now
		c[1] = c[1]*10 - 5; // y is from -5 to 5 now
		c[2] = c[2]*10 - 5; // z is from -5 to 5 now
		c[1] += a*3*(1-2*b); //now y changes from -5 to -2 on bottom and from 5 to 2 on top as x moves from 0 to 35
		c[2] += a*3*(1-2*q);
	}
	
	m->PrepareGeometricData(t);
	
	
	std::cout << "Saving output to " << argv[2] << std::endl;
	
	m->Save(argv[2]);
	
	delete m;
}


