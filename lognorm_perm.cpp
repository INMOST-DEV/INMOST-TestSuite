

#include "inmost.h"
#include <random>
using namespace INMOST;
std::string problem_name = "lognormal_scalar_permeability";

int main(int argc, char ** argv)
{
	if( argc < 3 )
	{
		std::cout << "Usage: " << argv[0] << " mesh mesh_out [mean=1] [variance=0.01]" << std::endl;
		return 0;
	}
	
	Mesh * m = new Mesh;
	m->SetFileOption("VERBOSITY","2");
	try{m->Load(argv[1]);} catch(...) { std::cout << "Cannot load the mesh " << argv[1] << std::endl; return -1;}


	
	if( m->HaveTag("PERM") ) m->DeleteTag(m->GetTag("PERM"));

	TagReal perm = m->CreateTag("PERM",DATA_REAL,CELL,NONE,1);
	double mean = 1.0, variance = 0.01;
	if( argc > 3 ) mean = atof(argv[3]);
	if( argc > 4 ) variance = atof(argv[4]);
    
    {
        Storage::bulk_array name = m->self()->BulkArray(m->CreateTag("PROBLEMNAME",DATA_BULK,MESH,NONE));
        name.replace(name.begin(),name.end(),problem_name.begin(),problem_name.end());
    }


	std::cout << "Setting lognormal scalar permeability with mean " << mean << " and variance " << variance << std::endl;

	std::random_device rd;
	std::mt19937 mt(rd());
	std::lognormal_distribution<double> dist(mean, variance);

	for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
		perm[*it] = dist(mt);
	std::cout << "Saving output to " << argv[2] << std::endl;
	m->Save(argv[2]);
	delete m;
}

