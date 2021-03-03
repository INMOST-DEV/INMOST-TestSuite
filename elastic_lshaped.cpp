#include "inmost.h"

using namespace INMOST;

bool print = false;
const Storage::real pi = 3.1415926535897932384626433832795;
const Storage::real eps = 1.0e-6;
std::string problem_name = "elastic_lshaped";
//defined from
// https://www.math.hu-berlin.de/~cc/cc_homepage/download/2002-AJ_CC_FS_KR-Matlab_Implementation_FEM_Elasticity.pdf


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
	
	//erase lower-left corner of the mesh
	for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
	{
		it->Centroid(c);
		if( c[0] < 0.5 && c[1] < 0.5 )
			it->Delete();
	}
	for(ElementType etype = FACE; etype >= NODE; etype = PrevElementType(etype) )
		for(Mesh::iteratorElement it = m->BeginElement(etype); it != m->EndElement(); ++it)
			if( it->nbAdjElements(NextElementType(etype)) == 0 ) it->Delete();
	//unite sliced parts
	if( m->HaveTag("SLICED") )
	{
		//Tag sliced = m.CreateTag("SLICED",DATA_BULK,FACE|EDGE|NODE,FACE|EDGE|NODE,1);
		Tag sliced = m->GetTag("SLICED");
		for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it) if( it->HaveData(sliced) )
		{
			ElementArray<Cell> cells = it->getCells();
			if( cells.size() > 1 ) Cell::UniteCells(cells,0);
		}
		for(Mesh::iteratorEdge it = m->BeginEdge(); it != m->EndEdge(); ++it) if( it->HaveData(sliced) )
		{
			ElementArray<Face> faces = it->getFaces();
			//check all faces are on one plane
			bool dosplit = faces.size() > 1;
			double nrm0[3], nrm1[3], dot;
			faces[0].UnitNormal(nrm0);
			for(int k = 1; k < faces.size() && dosplit; ++k)
			{
				faces[k].UnitNormal(nrm1);
				dot = nrm0[0]*nrm1[0]+nrm0[1]*nrm1[1]+nrm0[2]*nrm1[2];
				if( fabs(fabs(dot)-1.0) > 1.0e-6 ) dosplit = false;
			}
			if( dosplit ) Face::UniteFaces(faces,0);
		}
		for(Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); ++it) if( it->HaveData(sliced)  )
		{
			ElementArray<Edge> edges = it->getEdges();
			//check all edges are on one line
			bool dosplit = edges.size() > 1;
			double v0[3], v1[3], dot;
			v0[0] = edges[0].getEnd().Coords()[0]-edges[0].getBeg().Coords()[0];
			v0[1] = edges[0].getEnd().Coords()[1]-edges[0].getBeg().Coords()[1];
			v0[2] = edges[0].getEnd().Coords()[2]-edges[0].getBeg().Coords()[2];
			dot = sqrt(v0[0]*v0[0]+v0[1]*v0[1]+v0[2]*v0[2]);
			v0[0] /= dot;
			v0[1] /= dot;
			v0[2] /= dot;
			for(int k = 1; k < edges.size() && dosplit; ++k)
			{
				v1[0] = edges[k].getEnd().Coords()[0]-edges[k].getBeg().Coords()[0];
				v1[1] = edges[k].getEnd().Coords()[1]-edges[k].getBeg().Coords()[1];
				v1[2] = edges[k].getEnd().Coords()[2]-edges[k].getBeg().Coords()[2];
				dot = sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]);
				v1[0] /= dot;
				v1[1] /= dot;
				v1[2] /= dot;
				dot = v0[0]*v1[0] + v0[1]*v1[1] + v0[2]*v1[2];
				if( fabs(fabs(dot)-1.0) > 1.0e-6 ) dosplit = false;
			}
			if( dosplit ) Edge::UniteEdges(edges,0);
		}
	}
	//shift mesh be translation from (0.5,0.5,0.5) to (0,0,0)
	for(Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); ++it)
	{
		Storage::real_array c = it->Coords();
		c[0] -= 0.5;
		c[1] -= 0.5;
		c[2] -= 0.5;
	}
	//rotate mesh in Oxy plane by -45 degrees
	rMatrix R(3,3);
	R(0,0) = R(1,1) = cos(-45.0/180.0*pi);
	R(0,1) = -sin(-45.0/180.0*pi);
	R(1,0) = -R(0,1);
	R(2,2) = 1;
	for(Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); ++it)
	{
		Storage::real_array c = it->Coords();
		rMatrix xyz(c.data(),3,1);
		xyz = R*xyz;
		c[0] = xyz(0,0);
		c[1] = xyz(1,0);
		c[2] = xyz(2,0);
	}
	//rescale mesh in all directions by 2/sqrt(2)
	for(Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); ++it)
	{
		Storage::real_array c = it->Coords();
		c[0] *= 2.0/sqrt(2.0);
		c[1] *= 2.0/sqrt(2.0);
		c[2] *= 2.0/sqrt(2.0);
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
	Tag grad_val = m->CreateTag("REFERENCE_GRADIENT",DATA_REAL,CELL,NONE,9);

	{
		Storage::bulk_array name = m->self()->BulkArray(m->CreateTag("PROBLEMNAME",DATA_BULK,MESH,NONE));
		name.replace(name.begin(),name.end(),problem_name.begin(),problem_name.end());
	}
	
	rMatrix C(6,6); //define elastic tensor with lame coefficients
	
	//Lame constants, young modulus E = 100000
	double E = 100000;
	double nu = 0.3;
	double lambda = E*nu/((1+nu)*(1-2*nu));
	double mu = E/(2*(1+nu));
	
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
	
	//constants defining problem
	double alpha = 0.54448373678246386742074491849052719771862030029297; //solve alpha*sin(2*omega) + sin(2*omega*alpha) = 0
	double omega = 3.0*pi/4.0;
	double C1 = -cos((alpha+1)*omega)/cos((alpha-1)*omega);
	double C2 = 2*(lambda+2*mu)/(lambda+mu);
	
	for(Mesh::iteratorElement it = m->BeginElement(CELL|FACE); it != m->EndElement(); ++it)
	{
		it->Centroid(c);
		
		unknown x(c[0],0);
		unknown y(c[1],1);
		unknown z(c[2],2);
		variable r = sqrt(x*x+y*y);
		variable phi = atan2(y,x);
		//variable u_r = pow(r,alpha)/(2.0*mu)*(C1*(C2-alpha-1)*cos((alpha-1)*phi) - (alpha+1)*cos((alpha+1)*phi));
		//variable u_phi = pow(r,alpha)/(2.0*mu)*((alpha+1)*sin((alpha+1)*phi) + C1*(C2+alpha-1)*sin((alpha-1)*phi));
		if( phi < -pi || phi > pi ) std::cout << "phi is " << get_value(phi) << std::endl;
		variable u_r = pow(r,alpha)/(2.0*mu)*((C2-alpha-1)*C1*cos((alpha-1)*phi) - (alpha+1)*cos((alpha+1)*phi));
		variable u_phi = pow(r,alpha)/(2.0*mu)*((alpha+1)*sin((alpha+1)*phi) + (C2+alpha-1)*C1*sin((alpha-1)*phi));
	
		
		//std::cout << "radius: " << get_value(u_r) << " angle: " << get_value(u_phi) << std::endl;
		
		vMatrix sol(3,1);
		rMatrix epsilon(6,1), sigma(6,1);
		rMatrix N(3,6), sigma_value(6,1), flux(3,1);
		
		sol(0,0) = u_r*cos(phi) - u_phi*sin(phi);
		sol(1,0) = u_r*sin(phi) + u_phi*cos(phi);
		sol(2,0) = 0;
		
		epsilon(0,0) = sol(0,0).GetRow()[0]; //u_x
		epsilon(1,0) = sol(1,0).GetRow()[1]; //v_y
		epsilon(2,0) = sol(2,0).GetRow()[2]; //w_z
		epsilon(3,0) = sol(1,0).GetRow()[2] + sol(2,0).GetRow()[1]; //v_z + w_y
		epsilon(4,0) = sol(0,0).GetRow()[2] + sol(2,0).GetRow()[0]; //u_z + w_x
		epsilon(5,0) = sol(0,0).GetRow()[1] + sol(1,0).GetRow()[0]; //u_y + v_x
		
		
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
			rsol[0] = sol(0,0).GetValue();
			rsol[1] = sol(1,0).GetValue();
			rsol[2] = sol(2,0).GetValue();
			
			
			Storage::real_array s = it->RealArray(stress_val);
			s[0] = sigma(0,0);
			s[1] = sigma(1,0);
			s[2] = sigma(2,0);
			s[3] = sigma(3,0);
			s[4] = sigma(4,0);
			s[5] = sigma(5,0);
			
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
				if( c[2] > 1.0/sqrt(2.0)-eps || c[2] < -1.0/sqrt(2.0)+eps ) //pure roller on z
				{
					bc[0] = 0.0; //alpha
					bc[1] = 1.0; //beta
					bc[2] = 1.0; //projection, gamma
					//right hand side
					bc[3] = 0.0;
					bc[4] = 0.0;
					bc[5] = 0.0;
				}
				else //dirichlet
				{
					bc[0] = 1.0; //alpha
					bc[1] = 0.0; //beta
					bc[2] = 0.0; //projection, gamma
					//right hand side
					bc[3] = sol(0,0).GetValue();
					bc[4] = sol(1,0).GetValue();
					bc[5] = sol(2,0).GetValue();
				}
			}
			
		}
	}
	
	/*
	for(Mesh::iteratorElement it = m->BeginElement(CELL|FACE); it != m->EndElement(); ++it)
	{
		it->Centroid(c);
		
		double x = c[0];
		double y = c[1];
		double z = c[2];
		double r = sqrt(x*x+y*y);
		double phi = atan2(y,x);
		if( phi < -pi || phi > pi ) std::cout << "phi is " << phi << std::endl;
		double u_r = pow(r,alpha)/(2.0*mu)*((C2-alpha-1)*C1*cos((alpha-1)*phi) - (alpha+1)*cos((alpha+1)*phi));
		double u_phi = pow(r,alpha)/(2.0*mu)*((alpha+1)*sin((alpha+1)*phi) + (C2+alpha-1)*C1*sin((alpha-1)*phi));
		
		//std::cout << "radius: " << get_value(u_r) << " angle: " << get_value(u_phi) << std::endl;
		
		rMatrix sol(3,1);
		
		sol(0,0) = u_r*cos(u_phi);
		sol(1,0) = u_r*sin(u_phi);
		sol(2,0) = 0;
		
		
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
		}
		
		
		
		else if( it->GetElementType() == FACE )
		{
			if( it->getAsFace()->Boundary() )
			{
				Storage::real_array bc = it->RealArray(bndcond);
				// alpha*u + beta*(I-p*n*n^T)*N*sigma = gamma
				if( c[2] > 1.0/sqrt(2.0)-eps || c[2] < -1.0/sqrt(2.0)+eps ) //pure roller on z
				{
					bc[0] = 0.0; //alpha
					bc[1] = 1.0; //beta
					bc[2] = 1.0; //projection, gamma
					//right hand side
					bc[3] = 0.0;
					bc[4] = 0.0;
					bc[5] = 0.0;
				}
				else //dirichlet
				{
					bc[0] = 1.0; //alpha
					bc[1] = 0.0; //beta
					bc[2] = 0.0; //projection, gamma
					//right hand side
					bc[3] = sol(0,0);
					bc[4] = sol(1,0);
					bc[5] = sol(2,0);
				}
			}
			
		}
	}
	*/
	
	std::cout << "Saving output to " << argv[2] << std::endl;
	
	m->Save(argv[2]);
	
	delete m;
}


