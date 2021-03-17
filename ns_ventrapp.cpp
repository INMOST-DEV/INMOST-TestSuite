#include "inmost.h"

using namespace INMOST;

//setup BC for ventricular appendage emulation, generate mesh using gmsh with ns_ventrapp.geo file
//
// requires GMSH_TAGS on mesh, physical surface with tag 2 - inlet, with tag 3 - outlet
//
//setup on mesh:
//
//
// BOUNDARY_CONDITION_VELOCITY - 7 entries
// r = (a5,a6,a7)
// n.(a1*u + a2*t) = n.r
// (I - nn.)(a3*u + a4*t) = (I-nn.)r
//
// BOUNDARY_CONDITION_PRESSURE - 1 entry
//
// BOUNDARY_BLOOD - 1 entry
// specify inlet position



typedef Storage::real real;
const double eps = 1.0e-5;


int main(int argc, char ** argv)
{
	
	if( argc < 2 )
	{
		std::cout << "Usage: " << argv[0] << " mesh [p_in=1] [p_out=0] [mesh_out=grid_out.pmf]" << std::endl;
		return 0;
	}
	
	Mesh * m = new Mesh;
	m->SetFileOption("VERBOSITY","2");
	try
	{
		m->Load(argv[1]);
	}
	catch(...)
	{
		std::cout << "Cannot load the mesh " << argv[1] << std::endl;
		return -1;
	}
	
	//read params
	std::string fout = "grid_out.pmf";
	double p_in = 1;
	double p_out = 0;
	if( argc > 2 ) p_in  = atof(argv[2]);
	if( argc > 3 ) p_out = atof(argv[3]);
	if( argc > 4 ) fout  = std::string(argv[4]);
	
	
	if( !m->HaveTag("GMSH_TAGS") )
	{
		std::cout << "Mesh does not have GMSH_TAGS data" << std::endl;
		return -1;
	}
	
	//just for infos
	double cmax[3] = {-1.0e20,-1.0e20,-1.0e20}, cmin[3] = {1.0e20,1.0e20,1.0e20};
	for(Mesh::iteratorNode n = m->BeginNode(); n != m->EndNode(); ++n)
	{
		Storage::real_array c = n->Coords();
		for(int k = 0; k < 3; ++k)
		{
			if( cmax[k] < c[k] ) cmax[k] = c[k];
			if( cmin[k] > c[k] ) cmin[k] = c[k];
		}
	}
	
	char XYZ[3] = {'x','y','z'};
	for(int d = 0; d < 3; ++d)
		std::cout << XYZ[d] << " " << cmin[d] << ":" << cmax[d] << std::endl;
			
	
	
	//TagRealArray force = m->CreateTag("FORCE",DATA_REAL,CELL,NONE,3);
	TagRealArray  bc      = m->CreateTag("BOUNDARY_CONDITION_VELOCITY",DATA_REAL,FACE,FACE,7);
	TagReal       bcp     = m->CreateTag("BOUNDARY_CONDITION_PRESSURE",DATA_REAL,FACE,FACE,1);
	TagInteger    bb      = m->CreateTag("BOUNDARY_BLOOD",DATA_INTEGER,FACE,FACE,1);
	TagRealArray  uvw     = m->CreateTag("UVW",DATA_REAL,CELL,NONE,3); //keep zero
	TagReal       p       = m->CreateTag("P",DATA_REAL,CELL,NONE,1); //keep zero
	TagIntegerArray gmsh  = m->GetTag("GMSH_TAGS");
	
	if( !gmsh.isDefined(FACE) )
	{
		std::cout << "gmsh tag was not defined on face" << std::endl;
		return -1;
	}
	
	
	//~ bool onep = false;
	for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it) if( it->Boundary() )
	{
		Storage::real  n[3], c[3];
		it->UnitNormal(n);
		it->Centroid(c);
		bb[*it] = 0;
		if( gmsh[*it][1] == 3 ) // outflow
		{
			bcp[*it] = p_out;
			//~ if( !onep )
			/*
			{
				//~ onep = true;
				bc[*it][0] = 0;
				bc[*it][1] = 1;
				bc[*it][2] = 0;
				bc[*it][3] = 1;
				bc[*it][4] = 0;
				bc[*it][5] = 0;
				bc[*it][6] = 0;
			}
			*/
		}
		else if( gmsh[*it][1] == 2 ) //inflow
		{
			bb[*it] = 1;
			bcp[*it] = p_in;
			/*
			bc[*it][0] = 1;
			bc[*it][1] = 0;
			bc[*it][2] = 1;
			bc[*it][3] = 0;
			bc[*it][4] = 17.3 * c[1]*(3.0-c[1])/9.0 * 4.0 * (3.0 / 2.0); //parabolic profile with average of 17.3
			bc[*it][5] = 0;
			bc[*it][6] = 0;
			*/
		}
		/*
		else //slip wall
		{
			bc[*it][0] = 1;
			bc[*it][1] = 0;
			bc[*it][2] = 0;
			bc[*it][3] = 1;
			bc[*it][4] = 0;
			bc[*it][5] = 0;
			bc[*it][6] = 0;
		}
		*/
		else //no-slip walls
		{
			bc[*it][0] = 1;
			bc[*it][1] = 0;
			bc[*it][2] = 1;
			bc[*it][3] = 0;
			bc[*it][4] = 0;
			bc[*it][5] = 0;
			bc[*it][6] = 0;
		}
	}
	
	std::cout << "Saving output to " << fout << std::endl;
	
	m->Save(fout);
	
	return 0;
}
