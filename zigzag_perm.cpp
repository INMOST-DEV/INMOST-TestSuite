

#include "inmost.h"

using namespace INMOST;
std::string problem_name = "zigzag_permeability";
const double pi = 3.1415926535897932384626433832795;

int main(int argc, char ** argv)
{
	if( argc < 3 )
	{
		std::cout << "Usage: " << argv[0] << " mesh mesh_out [Kx=1000] [Ky=10]" << std::endl;
		return 0;
	}
	
	Mesh * m = new Mesh;
	m->SetFileOption("VERBOSITY","2");
	try{m->Load(argv[1]);} catch(...) { std::cout << "Cannot load the mesh " << argv[1] << std::endl; return -1;}


	double max[3] = {-1.0e20, -1.0e20, -1.0e20}, min[3] = {1.0e20,1.0e20,1.0e20};
	Storage::real c[3] = {0,0,0};
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
	if( max[2] <= min[2] ) min[2] = max[2] = 0.0; //2d mesh

	Tag material;
	if( m->HaveTag("MATERIAL") ) material = m->GetTag("MATERIAL");
	if( m->HaveTag("PERM") ) m->DeleteTag(m->GetTag("PERM"));

	TagRealArray tensor = m->CreateTag("PERM",DATA_REAL,CELL,NONE,6);

	Storage::real Kx = 1000.0, Ky = 10.0, Kz = 10.0;

	if( argc > 3 ) Kx = atof(argv[3]);
	if( argc > 4 ) Ky = atof(argv[4]);
    
    {
        Storage::bulk_array name = m->self()->BulkArray(m->CreateTag("PROBLEMNAME",DATA_BULK,MESH,NONE));
        name.replace(name.begin(),name.end(),problem_name.begin(),problem_name.end());
    }


	std::cout << "Setting zigzag permeability" << std::endl;

	for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
	{
		it->Centroid(c);
		Storage::real alpha = (c[0]-min[0])/(max[0]-min[0]);
		Storage::real beta  = (c[1]-min[1])/(max[1]-min[1]);
		Storage::real angle = 0;
		Storage::real_array perm = it->RealArrayDF(tensor);
		std::fill(perm.begin(),perm.end(),0.0);

		if (alpha + beta < 0.5)
			angle = 45;
		else if (alpha + beta < 1.0)
			angle = 0;
		else if (alpha + beta < 1.5)
			angle = 90;
		else angle = 45;
		
		Storage::real acos = cos(angle / 180.0 * pi);
		Storage::real asin = sin(angle / 180.0 * pi);
		perm[0] = Kx * acos * acos + Ky * asin * asin;
		perm[1] = (Kx - Ky) * acos * asin;
		perm[2] = 0.0;
		perm[3] = Ky * acos * acos + Kx * asin * asin;
		perm[4] = 0.0;
		perm[5] = Kz;
	}
	std::cout << "Saving output to " << argv[2] << std::endl;
	m->Save(argv[2]);
	delete m;
}

