#include "inmost.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


using namespace INMOST;
std::string problem_name = "single_well";

#define V_ID(x, y, z) ((x)*n*n + (y)*n + (z))


void prepare_tensor1(Storage::real t[21])
{
	//diagonal
	t[0] = 13;
	t[6] = 12.5;
	t[11] = 48;
	t[15] = 20;
	t[18] = 20;
	t[20] = 36;
	//off-diagonal top line
	t[1] = -6;
	t[2] = -12;
	t[3] = -7;
	t[4] = 10;
	t[5] = 0;
	//off-diagonal second line
	t[7] = -12;
	t[8] = 7;
	t[9] = -10;
	t[10] = 0;
	//off-diagonal third line
	t[12] = 0;
	t[13] = 0;
	t[14] = 0;
	//off-diagonal forth line
	t[16] = 0;
	t[17] = -20;
	//last
	t[19] = -14;
}

void prepare_tensor2(Storage::real t[21])
{
	//diagonal
	t[0] = 32;
	t[6] = 29;
	t[11] = 18;
	t[15] = 21;
	t[18] = 25;
	t[20] = 34;
	//off-diagonal top line
	t[1] = -8;
	t[2] = -14;
	t[3] = -9;
	t[4] = 11;
	t[5] = 0.9;
	//off-diagonal second line
	t[7] = -11;
	t[8] = 8;
	t[9] = -9;
	t[10] = 0.1;
	//off-diagonal third line
	t[12] = 0.1;
	t[13] = 1;
	t[14] = 2;
	//off-diagonal forth line
	t[16] = 2;
	t[17] = -11;
	//last
	t[19] = -12;
}

int main(int argc, char *argv[])
{
	
	if (argc < 2)
	{
  std::cout << "Usage: " << argv[0] << " nx|grid_name [alpha=0.4] [outer_boundary_pressure=0.0 0.0 0.0] [inner_boundary_pressure=1.0 1.0 1.0] [cut_grid=1]" << std::endl;
  return -1;
	}
	
	Mesh * mesh;
	double alpha=0.05;
	double h;
	double outer_boundary_displacement[3] = {0,0,0};
	double inner_boundary_displacement[3] = {1,1,1};
	int cut_grid = 1;
	int n = 20 + 1;
	
	mesh = new Mesh();
	
	if (argc > 1)
	{
		if( atoi(argv[1]) )
		{
			n = atoi(argv[1])+1;
			h = 1.0 / (n-1);
			std::cout << "Mesh: cube " << n << std::endl;
		}
		else
		{
			mesh->Load(argv[1]);
			std::cout << "Mesh: " << argv[1] << std::endl;
			n = 0;
			h = 0;
			for(Mesh::iteratorCell e = mesh->BeginCell(); e != mesh->EndCell(); ++e)
			{
				Storage::real maxmin[6];
				maxmin[0] = -1e20;
				maxmin[1] = 1e20;
				maxmin[2] = -1e20;
				maxmin[3] = 1e20;
				maxmin[4] = -1e20;
				maxmin[5] = 1e20;
				ElementArray<Node> nodes = e->getNodes();
				for (ElementArray<Node>::iterator it = nodes.begin(); it != nodes.end(); it++)
				{
					Storage::real_array c = it->Coords();
					for (int i = 0; i < (int)c.size(); i++) {
						if (maxmin[2 * i + 0] < c[i]) maxmin[2 * i + 0] = c[i]; //max
						if (maxmin[2 * i + 1] > c[i]) maxmin[2 * i + 1] = c[i]; //min
					}
					if( c.size() < 3 )
					{
						for(int i = c.size(); i < 3; i++)
						{
							maxmin[2*i+0] = maxmin[2*i+1] = 0;
						}
					}
				}
				h = std::max(h,maxmin[0]-maxmin[1]);
				h = std::max(h,maxmin[2]-maxmin[3]);
				//h = std::max(h,maxmin[4]-maxmin[5]);
			}
		}
		
	}
	
	std::cout << "Mesh radius: " << h << std::endl;
	if( argc > 2 )
		alpha = atof(argv[2]);
	
	std::cout << "Alpha: " << alpha << std::endl;
	
	if( argc > 3 ) outer_boundary_displacement[0] = atof(argv[3]);
	if( argc > 4 ) outer_boundary_displacement[1] = atof(argv[4]);
	if( argc > 5 ) outer_boundary_displacement[2] = atof(argv[5]);
	
	if( argc > 6 ) inner_boundary_displacement[0] = atof(argv[6]);
	if( argc > 7 ) inner_boundary_displacement[1] = atof(argv[7]);
	if( argc > 8 ) inner_boundary_displacement[2] = atof(argv[8]);
	
	std::cout << "outer displacement: ";
	std::cout << outer_boundary_displacement[0] << ",";
	std::cout << outer_boundary_displacement[1] << ",";
	std::cout << outer_boundary_displacement[2] << ".";
	std::cout << std::endl;
	
	std::cout << "inner displacement: ";
	std::cout << inner_boundary_displacement[0] << ",";
	std::cout << inner_boundary_displacement[1] << ",";
	std::cout << inner_boundary_displacement[2] << ".";
	std::cout << std::endl;
	
	if( argc > 9 )
		cut_grid = atoi(argv[9]);
	
	if( cut_grid )
		std::cout << "Cutting center of the grid." << std::endl;
	
	
	
	
	srand(0);//(unsigned)time(NULL));
	
	if( n )
	{
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				for (int k = 0; k < n; k++)
				{
					Storage::real xyz[3];
					bool mark = false;
					xyz[0] = i * 1.0 / (n - 1);
					xyz[1] = j * 1.0 / (n - 1);
					xyz[2] = k * 1.0 / (n - 1);
					Node c = mesh->CreateNode(xyz);
					if (c->LocalID() != V_ID(i, j, k)) printf("v_id = %d, [i,j,k] = %d\n", c->LocalID(), V_ID(i, j, k));
				}
			}
		}
		
		
		const INMOST_DATA_INTEGER_TYPE nvf[24] = { 2, 3, 1, 0, 4, 5, 7, 6, 0, 1, 5, 4, 3, 2, 6, 7, 2, 0, 4, 6, 1, 3, 7, 5 };
		const INMOST_DATA_INTEGER_TYPE numnodes[6] = { 4, 4, 4, 4, 4, 4 };
		for (int i = 1; i < n; i++)
		{
			for (int j = 1; j < n; j++)
			{
				for (int k = 1; k < n; k++)
				{
					
					ElementArray<Node> verts(mesh,8);
					verts[0] = mesh->NodeByLocalID(V_ID(i - 1, j - 1, k - 1));
					verts[1] = mesh->NodeByLocalID(V_ID(i - 0, j - 1, k - 1));
					verts[2] = mesh->NodeByLocalID(V_ID(i - 1, j - 0, k - 1));
					verts[3] = mesh->NodeByLocalID(V_ID(i - 0, j - 0, k - 1));
					verts[4] = mesh->NodeByLocalID(V_ID(i - 1, j - 1, k - 0));
					verts[5] = mesh->NodeByLocalID(V_ID(i - 0, j - 1, k - 0));
					verts[6] = mesh->NodeByLocalID(V_ID(i - 1, j - 0, k - 0));
					verts[7] = mesh->NodeByLocalID(V_ID(i - 0, j - 0, k - 0));
					
					mesh->CreateCell(verts,nvf,numnodes,6);
				}
			}
		}
	}
	
	
	if( cut_grid )
	{
		for (Mesh::iteratorCell it = mesh->BeginCell(); it != mesh->EndCell(); ++it)
		{
			Storage::real cnt[3];
			it->Centroid(cnt);
			if( cnt[0] > 0.4 && cnt[0] < 0.6 && cnt[1] > 0.4 && cnt[1] < 0.6 && cnt[2] > 0.4 && cnt[2] < 0.6 )
				it->Delete();
		}
		
		for (Mesh::iteratorElement it = mesh->BeginElement(FACE|EDGE|NODE); it != mesh->EndElement(); ++it)
			if( it->nbAdjElements(CELL) == 0 ) it->Delete();
	}
	
	if( alpha > 0.0 )
	{
		for( Mesh::iteratorNode it = mesh->BeginNode(); it != mesh->EndNode(); ++it) if( !it->Boundary() )
		{
			it->Coords()[0] += alpha*h*0.5*(2.0*rand()/RAND_MAX-1.0);
			it->Coords()[1] += alpha*h*0.5*(2.0*rand()/RAND_MAX-1.0);
			it->Coords()[2] += alpha*h*0.5*(2.0*rand()/RAND_MAX-1.0);
		}
		{
			std::stringstream str;
			str << problem_name << "_disturb_" << alpha;
			problem_name = str.str();
		}
	}
	
	mesh->ReorderEmpty(CELL|FACE|EDGE|NODE);
	
	printf("nodes: %d edges: %d faces: %d cells: %d\n", mesh->NumberOfNodes(), mesh->NumberOfEdges(), mesh->NumberOfFaces(), mesh->NumberOfCells());
	
	
	{
		Storage::bulk_array name = mesh->self()->BulkArray(mesh->CreateTag("PROBLEMNAME",DATA_BULK,MESH,NONE));
		name.replace(name.begin(),name.end(),problem_name.begin(),problem_name.end());
	}
	
	Tag bndcond = mesh->CreateTag("BOUNDARY_CONDITION",DATA_REAL,FACE|NODE,FACE|NODE,6);
	
	int numinner = 0, numouter = 0;
	const double eps = 1.0e-6;
	for(Mesh::iteratorElement it = mesh->BeginElement(FACE); it != mesh->EndElement(); ++it) if( it->Boundary() )
	{
		Storage::real cnt[3];
		it->Centroid(cnt);
		
		
		//internal part bounary
		//internal part bounary
		//cnt[0] > 0.1 && cnt[0] < 0.9 && cnt[1] > 0.1 && cnt[1] < 0.9 && cnt[2] > 0.1 && cnt[2] < 0.9 )
		if( (cnt[0] > 0.4-eps && cnt[0] < 0.6+eps && cnt[1] > 0.4-eps && cnt[1] < 0.6+eps && cnt[2] > 0.4-eps && cnt[2] < 0.6+eps))
	{
			Storage::real_array bnd = it->RealArray(bndcond);
			bnd[0] = 1.0; //dirichlet
			bnd[1] = 0.0;
			bnd[2] = 0.0;
			bnd[3] = inner_boundary_displacement[0];
			bnd[4] = inner_boundary_displacement[1];
			bnd[5] = inner_boundary_displacement[2];
			if( it->GetElementType() == FACE ) numinner++;

		}
		else if( (cnt[0] < eps || cnt[0] > 1.0-eps || cnt[1] < eps || cnt[1] > 1 - eps) || (cnt[2] < 1-eps || cnt[2] > 1 - eps) )
		{
			//std::cout << ElementTypeName(it->GetElementType()) << " " << cnt[0] << " " << cnt[1] << " " << cnt[2] << " " << (cnt[0] < eps || cnt[0] > 1.0-eps || cnt[1] < eps || cnt[1] > eps) << std::endl;
			Storage::real_array bnd = it->RealArray(bndcond);
			bnd[0] = 1.0; //dirichlet
			bnd[1] = 0.0;
			bnd[2] = 0.0;
			bnd[3] = outer_boundary_displacement[0];
			bnd[4] = outer_boundary_displacement[1];
			bnd[5] = outer_boundary_displacement[2];
			if( it->GetElementType() == FACE )  numouter++;
		}
	}
	std::cout << "Outer faces " << numouter << " inner faces " << numinner << std::endl;
	
	
	Storage::real mat1[21], mat2[21], cnt[3];
	prepare_tensor1(mat1);
	prepare_tensor2(mat2);
	
	
	
	Tag tensor = mesh->CreateTag("ELASTIC_TENSOR", DATA_REAL, CELL, NONE, 21);
	for (Mesh::iteratorCell it = mesh->BeginCell(); it != mesh->EndCell(); ++it)
	{
		it->Centroid(cnt);
		if( cnt[0] < 0.5 )
			memcpy(&it->RealArray(tensor)[0],mat1,sizeof(Storage::real)*21);
		else
			memcpy(&it->RealArray(tensor)[0],mat2,sizeof(Storage::real)*21);
	}
	
	printf("I'm ready!\n");
	
	//mesh->Save("grid.vtk");
	mesh->Save("grid_out.pmf");
	//mesh->Save("grid.gmv");
	
	printf("File written!\n");
	
	delete mesh;
	return 0;
}
