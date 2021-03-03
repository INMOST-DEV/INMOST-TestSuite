#include "inmost.h"

using namespace INMOST;

//#define TEST_INTRP

bool print = false;
//const Storage::real pi = 3.1415926535897932384626433832795;
const Storage::real eps = 1.0e-6;
std::string problem_name = "elastic_analytic_test_3";

//discontinuity plane equation ax+by+cz=d
double pnrm[3] = {2.0/7.0,6.0/7.0,3.0/7.0};
//double pnrm[3] = {1,0,0};
double pcnt[3] = {0.5,0.5,0.5};

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
}


void GetNormalProjection(const rMatrix & C, const rMatrix & n, rMatrix & L)
{
	
	rMatrix K(3,3); //co-normal matrix
	// extract \mathbb{K}_1
	K(0,0) = C(0,0);
	K(1,1) = C(5,5);
	K(2,2) = C(4,4);
	K(0,1) = K(1,0) = C(0,5);
	K(0,2) = K(2,0) = C(0,4);
	K(1,2) = K(2,1) = C(4,5);
	L(0,0) = n.DotProduct(K*n);
	
	// extract \mathbb{K}_2
	K(0,0) = C(5,5);
	K(1,1) = C(1,1);
	K(2,2) = C(3,3);
	K(0,1) = K(1,0) = C(1,5);
	K(0,2) = K(2,0) = C(3,5);
	K(1,2) = K(2,1) = C(1,3);
	L(1,1) = n.DotProduct(K*n);
	
	// extract \mathbb{K}_3
	K(0,0) = C(4,4);
	K(1,1) = C(3,3);
	K(2,2) = C(2,2);
	K(0,1) = K(1,0) = C(3,4);
	K(0,2) = K(2,0) = C(2,4);
	K(1,2) = K(2,1) = C(2,3);
	L(2,2) = n.DotProduct(K*n);
	
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
}


void GetTensors(const rMatrix & C, rMatrix K[6])
{
	// extract \mathbb{K}_1
	K[0].Resize(3,3);
	K[0](0,0) = C(0,0);
	K[0](1,1) = C(5,5);
	K[0](2,2) = C(4,4);
	K[0](0,1) = K[0](1,0) = C(0,5);
	K[0](0,2) = K[0](2,0) = C(0,4);
	K[0](1,2) = K[0](2,1) = C(4,5);
	
	// extract \mathbb{K}_2
	K[1].Resize(3,3);
	K[1](0,0) = C(5,5);
	K[1](1,1) = C(1,1);
	K[1](2,2) = C(3,3);
	K[1](0,1) = K[1](1,0) = C(1,5);
	K[1](0,2) = K[1](2,0) = C(3,5);
	K[1](1,2) = K[1](2,1) = C(1,3);
	
	// extract \mathbb{K}_3
	K[2].Resize(3,3);
	K[2](0,0) = C(4,4);
	K[2](1,1) = C(3,3);
	K[2](2,2) = C(2,2);
	K[2](0,1) = K[2](1,0) = C(3,4);
	K[2](0,2) = K[2](2,0) = C(2,4);
	K[2](1,2) = K[2](2,1) = C(2,3);
	
	
	// extract \mathbb{K}_4
	K[5].Resize(3,3);
	K[5](0,0) = C(0,5);
	K[5](0,1) = C(5,5);
	K[5](0,2) = C(4,5);
	K[5](1,0) = C(0,1);
	K[5](1,1) = C(1,5);
	K[5](1,2) = C(1,4);
	K[5](2,0) = C(0,3);
	K[5](2,1) = C(3,5);
	K[5](2,2) = C(3,4);
	
	
	
	// extract \mathbb{K}_5
	K[4].Resize(3,3);
	K[4](0,0) = C(0,4);
	K[4](0,1) = C(4,5);
	K[4](0,2) = C(4,4);
	K[4](1,0) = C(0,3);
	K[4](1,1) = C(3,5);
	K[4](1,2) = C(3,4);
	K[4](2,0) = C(0,2);
	K[4](2,1) = C(2,5);
	K[4](2,2) = C(2,4);
	
	
	// extract \mathbb{K}_6
	K[3].Resize(3,3);
	K[3](0,0) = C(4,5);
	K[3](0,1) = C(1,4);
	K[3](0,2) = C(3,4);
	K[3](1,0) = C(3,5);
	K[3](1,1) = C(1,3);
	K[3](1,2) = C(3,3);
	K[3](2,0) = C(2,5);
	K[3](2,1) = C(1,2);
	K[3](2,2) = C(2,3);
}

/*
 void GenTensor(rMatrix & C, int seed, bool print = false)
 {
	rMatrix U,S,V;
	srand(seed);
	//construct some symmetric C
	for(int i = 0; i < 6; ++i)
		for(int j = i; j < 6; ++j)
			C(i,j) = C(j,i) = ceil(rand()/(1.*RAND_MAX)*100.0);
	if( print )
	{
		std::cout << "Original matrix: " << std::endl;
		C.Print();
	}
	//perform SVD
	C.SVD(U,S,V);
	//check column of U is the same sign as row of V^T
	SVD2Eigen(U,S,V);
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
	for(int i = 0; i < 6; ++i)
		if( S(i,i) < 0.0 )
			S(i,i) = std::abs(S(i,i));
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
		std::cout << "Original U-V: " << std::endl;
		(U-V).Print();
	}
	
}
 */

bool CheckEigen(const rMatrix & U, const rMatrix & S, const rMatrix & V)
{
	bool negative = false;
	for(int i = 0; i < V.Cols() && !negative; ++i)
	{
		double dot = 0.0;
		for(int j = 0; j < V.Rows(); ++j)
			dot += U(j,i)*V(j,i);
		if( S(i,i) * (dot < 0.0 ? -1.0 : 1.0) < 0.0 )
			negative = true;
	}
	return negative;
}


int GenTensor(rMatrix & C, int seed, bool print = false)
{
	rMatrix U(6,6),S(6,6),V(6,6);
	for(int k = seed; k < INT_MAX; ++k)
	{
		std::cout << k << "\r";
		std::cout.flush();
		srand(k);
		//construct some symmetric C
		for(int i = 0; i < 6; ++i)
			for(int j = i; j < 6; ++j)
				C(i,j) = C(j,i) = ceil(rand()/(1.*RAND_MAX)*100.0);
		//perform SVD
		C.SVD(U,S,V);
		//check column of U is the same sign as row of V^T
		if( CheckEigen(U,S,V) )
			continue;
		else
		{
			SVD2Eigen(U,S,V);
			std::cout << "eigenvalues: ";
			for(int q = 0; q < 6; ++q)
				std::cout << S(q,q) << " ";
			std::cout << std::endl;
			return k;
		}
	}
	return -1;
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
	//if( m->HaveTag("MATERIAL") ) m->DeleteTag(m->GetTag("MATERIAL"));
	
	
	
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
	
	
	rMatrix C1(6,6), C2(6,6), G1(3,3), G2(3,3), d1(3,1), d2(3,1), M(9,6);
	rMatrix sigma(6,1), np(pnrm,3,1);
	rMatrix NP(3,6), B(3,3);
	rMatrix Q(3,3), xp(pcnt,3,1);
	const rMatrix I = rMatrix::Unit(3);
	
	//Voigt notation for normal matrix
	// NP*sigma = E : 0.5(G+G^T)np
	NP(0,0) = np(0,0);
	NP(1,1) = np(1,0);
	NP(2,2) = np(2,0);
	NP(0,4) = np(2,0);
	NP(0,5) = np(1,0);
	NP(1,3) = np(2,0);
	NP(1,5) = np(0,0);
	NP(2,3) = np(1,0);
	NP(2,4) = np(0,0);
	
	double vB[] =
	{
		1,0,0,0,0,0,
		0,0,0,0,0,1,
		0,0,0,0,1,0,
		0,0,0,0,0,1,
		0,1,0,0,0,0,
		0,0,0,1,0,0,
		0,0,0,0,1,0,
		0,0,0,1,0,0,
		0,0,1,0,0,0
	};
	
	M = rMatrix(vB,9,6);
	
	/*
	93         46         22         13         72         35
	46         95         41         62         56         24
	22         41         89         25         33         21
	13         62         25         87         13         25
	72         56         33         13         99         57
	35         24         21         25         57         78
	 eigenvalues: 280.365 107.553 65.4832 56.2158 24.6995 6.68368
	*/
	double C1v[] =
	{
		93,46,22,13,72,35,
		   95,41,64,56,24,
		      89,25,33,21,
		         87,13,25,
		            99,57,
		               78
	};
	/*
	 81         33          5         14         37         32
	 33        100         26         58         58         73
	 5         26         90         23         60         28
	 14         58         23         77         29          3
	 37         58         60         29         84         60
	 32         73         28          3         60         99
	 eigenvalues: 283.182 89.3622 84.024 60.4678 12.1471 1.81666
	 */
	double C2v[] =
	{
		81,33,5,14,37,32,
		   100,26,58,58,73,
		       90,23,60,28,
		          77,29,3,
		             84,60,
		                99
	};
	
	//int q = 0;
	//q = GenTensor(C1,q)+1;
	//q = GenTensor(C2,q)+1;
	C1 = rMatrix::FromTensor(C1v,21,6);
	C2 = rMatrix::FromTensor(C2v,21,6);
	
	GetNormalProjection(C2,np,B);

	double G1v[9] = 
	{
		  1    ,     14    ,     76  ,
        54    ,     22   ,       5 ,
        68    ,     94    ,     39
	};

	double d1v[3] =
	{
		 46,
        68, 
        52
	};
	/*
	for(int i = 0; i < 3; ++i)
	{
		for(int j = 0; j < 3; ++j)
			G1(i,j) = ceil(rand()/(1.*RAND_MAX)*100.0);
		d1(i,0) = ceil(rand()/(1.*RAND_MAX)*100.0);
	}
	*/
	G1 = rMatrix(G1v,3,3);
	d1 = rMatrix(d1v,3,1);
	sigma(0,0) = G1(0,0);
	sigma(1,0) = G1(1,1);
	sigma(2,0) = G1(2,2);
	sigma(3,0) = G1(1,2) + G1(2,1);
	sigma(4,0) = G1(0,2) + G1(2,0);
	sigma(5,0) = G1(0,1) + G1(1,0);
	Q = B.Invert()*NP*(C1-C2)*sigma*np.Transpose();
	G2 = G1 + Q;
	d2 = d1 - Q*xp;
	// G1*xp + d1 = G2*xp+d2 = G1*xp + Q*xp + d1 - Q*xp
	
	std::cout << "C1:" << std::endl;
	C1.Print();
	std::cout << "G1:" << std::endl;
	G1.Print();
	std::cout << "d1:" << std::endl;
	d1.Print();
	
	//std::cout << std::scientific;
	std::cout << "C2:" << std::endl;
	C2.Print();
	std::cout << "G2:" << std::endl;
	G2.Print();
	std::cout << "d2:" << std::endl;
	d2.Print();
	
	std::cout << "B*49:" << std::endl;
	(B*49.).Print();
	std::cout << "B^{-1}:" << std::endl;
	((B*49.).Invert()).Print();
	std::cout << "np*7:" << std::endl;
	(np*7.).Print();
	std::cout << "np*np^T*49:" << std::endl;
	(np*np.Transpose()*49.).Print();
	std::cout << "Q*49:" << std::endl;
	(Q*49.).Print();
	std::cout << "Q*xp*49*2:" << std::endl;
	(Q*xp*49.*2.).Print();
	std::cout << "M" <<std::endl;
	M.Print();
	std::cout << "I" << std::endl;
	I.Print();
	std::cout << "(I \\otimes n^T)*M*C1*M^T*G1" << std::endl;
	(I.Kronecker(np.Transpose())*M*C1*M.Transpose()*rMatrix(G1.data(),9,1)).Print();
	std::cout << "(I \\otimes n^T)*M*C2*M^T*G2" << std::endl;
	(I.Kronecker(np.Transpose())*M*C2*M.Transpose()*rMatrix(G2.data(),9,1)).Print();
	std::cout << "UVW1:" <<std::endl;
	(d1 + G1*xp).Print();
	std::cout << "UVW2:" <<std::endl;
	(d2 + G2*xp).Print();
	std::cout << "UVW1':" <<std::endl;
	(d1 + I.Kronecker(xp.Transpose())*M*(M.Transpose()*M).Invert()*M.Transpose()*rMatrix(G1.data(),9,1)).Print();
	std::cout << "UVW2':" <<std::endl;
	(d2 + I.Kronecker(xp.Transpose())*M*(M.Transpose()*M).Invert()*M.Transpose()*rMatrix(G2.data(),9,1)).Print();
	std::cout << "P*G1:" << std::endl;
	(M*(M.Transpose()*M).Invert()*M.Transpose()*rMatrix(G1.data(),9,1)).Transpose().Print();
	std::cout << "P*G2:" << std::endl;
	(M*(M.Transpose()*M).Invert()*M.Transpose()*rMatrix(G2.data(),9,1)).Transpose().Print();
	//d2 = d1 + I.Kronecker(xp.Transpose())*M*(M.Transpose()*M).Invert()*M.Transpose()*rMatrix((G1-G2).data(),9,1);
	//std::cout << "UVW2' with new d2:" <<std::endl;
	//(d2 + I.Kronecker(xp.Transpose())*M*(M.Transpose()*M).Invert()*M.Transpose()*rMatrix(G2.data(),9,1)).Print();
	//Test for staggered scheme interpolation
	std::cout << " B (B^T B)^{-1} B^T " << std::endl;
	(M*(M.Transpose()*M).Invert()*M.Transpose()).Print();
	std::cout << " I - B (B^T B)^{-1} B^T " << std::endl;
	(rMatrix::Unit(9) - M*(M.Transpose()*M).Invert()*M.Transpose()).Print();
	std::cout << " (B (B^T B)^{-1} B^T)^{-1} " << std::endl;
	(M*(M.Transpose()*M).Invert()*M.Transpose()).PseudoInvert().Print();
	std::cout << " N^T B (B^T B)^{-1} B^T N * 49" << std::endl;
	(I.Kronecker(np.Transpose())*M*(M.Transpose()*M).Invert()*M.Transpose()*I.Kronecker(np)*49).Print();
#if defined(TEST_INTRP)
	C1 = C1.Invert();
	C2 = C2.Invert();
	rMatrix U1,S1,V1, U2,S2,V2, w1(6,1), w2(6,1), B1(3,3), B2(3,3), C1r, C2r, E1r, e1(6,1), E2r, e2(6,1);
	double p1,p2, K1, K2;
	C1.SVD(U1,S1,V1);
	C2.SVD(U2,S2,V2);
	for(int k = 0; k < 6; ++k)
	{
		w1(k,0) = U1(k,5)*sqrt(3.0);
		w2(k,0) = U2(k,5)*sqrt(3.0);
	}
	C1r = (rMatrix::Unit(6)-w1*w1.Transpose()/3.0)*C1;
	C2r = (rMatrix::Unit(6)-w2*w2.Transpose()/3.0)*C2;
	E1r = (rMatrix::Unit(6)-w1*w1.Transpose()/3.0)*(C1r+w1*w1.Transpose()/3.0).Invert();
	E2r = (rMatrix::Unit(6)-w2*w2.Transpose()/3.0)*(C2r+w2*w2.Transpose()/3.0).Invert();
	B1(0,0) = w1(0,0);
	B1(1,1) = w1(1,0);
	B1(2,2) = w1(2,0);
	B1(0,1) = B1(1,0) = w1(5,0);
	B1(0,2) = B1(2,0) = w1(4,0);
	B1(1,2) = B1(2,1) = w1(3,0);
	B2(0,0) = w2(0,0);
	B2(1,1) = w2(1,0);
	B2(2,2) = w2(2,0);
	B2(0,1) = B2(1,0) = w2(5,0);
	B2(0,2) = B2(2,0) = w2(4,0);
	B2(1,2) = B2(2,1) = w2(3,0);
	std::cout << "U1:" << std::endl;
	U1.Print();
	std::cout << "w1:" << std::endl;
	w1.Print();
	std::cout << "|w1|:" << w1.DotProduct(w1) << std::endl;
	std::cout << "C1r:" << std::endl;
	C1r.Print();
	std::cout << "E1r:" << std::endl;
	E1r.Print();
	std::cout << "B1:" << std::endl;
	B1.Print();
	
	
	C1r.SVD(U1,S1,V1);
	std::cout << "SC1r:" << std::endl;
	S1.Print();
	
	E1r.SVD(U1,S1,V1);
	std::cout << "SE1r:" << std::endl;
	S1.Print();
	
	B1.SVD(U1,S1,V1);
	SVD2Eigen(U1,S1,V1);
	std::cout << "B_S1" << std::endl;
	S1.Print();
	
	
	e1(0,0) = G1(0,0);
	e1(1,0) = G1(1,1);
	e1(2,0) = G1(2,2);
	e1(3,0) = G1(1,2) + G1(2,1);
	e1(4,0) = G1(0,2) + G1(2,0);
	e1(5,0) = G1(0,1) + G1(1,0);
	
	std::cout << "e1:" << std::endl;
	e1.Print();
	
	K1 = (w1.DotProduct(C1*w1));
	p1 = w1.DotProduct(e1);
	
	std::cout << "K1: " << K1 << std::endl;
	std::cout << "p1: " << p1 << std::endl;
	
	
	std::cout << "U2:" << std::endl;
	U2.Print();
	std::cout << "w2:" << std::endl;
	w2.Print();
	std::cout << "|w2|:" << w2.DotProduct(w2) << std::endl;
	std::cout << "C2r:" << std::endl;
	C2r.Print();
	std::cout << "E2r:" << std::endl;
	E2r.Print();
	std::cout << "B2:" << std::endl;
	B2.Print();
	
	
	
	C2r.SVD(U2,S2,V2);
	std::cout << "SC2r:" << std::endl;
	S2.Print();
	
	E2r.SVD(U2,S2,V2);
	std::cout << "SE2r:" << std::endl;
	S2.Print();
	
	
	B2.SVD(U2,S2,V2);
	SVD2Eigen(U2,S2,V2);
	std::cout << "B_S2" << std::endl;
	S2.Print();
	
	
	e2(0,0) = G2(0,0);
	e2(1,0) = G2(1,1);
	e2(2,0) = G2(2,2);
	e2(3,0) = G2(1,2) + G2(2,1);
	e2(4,0) = G2(0,2) + G2(2,0);
	e2(5,0) = G2(0,1) + G2(1,0);
	std::cout << "e2:" << std::endl;
	e2.Print();
	
	K2 = (w2.DotProduct(C2*w2));
	p2 = w2.DotProduct(e2);
	
	std::cout << "K2: " << K2 << std::endl;
	std::cout << "p2: " << p2 << std::endl;
	
	Q = B2.Invert()*B1;
	
	std::cout << "Q:" << std::endl;
	Q.Print();
	
	rMatrix KK1[6], KK2[6], KK1n[9], KK2n[9];
	GetTensors(E1r,KK1);
	GetTensors(E2r,KK2);
	
	rMatrix T1(3,3), T2(3,3);
	
	T1(0,0) = np.DotProduct(KK1[0]*np);
	T1(1,1) = np.DotProduct(KK1[1]*np);
	T1(2,2) = np.DotProduct(KK1[2]*np);
	T1(0,1) = T1(1,0) = np.DotProduct(KK1[5]*np);
	T1(0,2) = T1(2,0) = np.DotProduct(KK1[4]*np);
	T1(1,2) = T1(2,1) = np.DotProduct(KK1[3]*np);
	
	T2(0,0) = np.DotProduct(KK2[0]*np);
	T2(1,1) = np.DotProduct(KK2[1]*np);
	T2(2,2) = np.DotProduct(KK2[2]*np);
	T2(0,1) = T2(1,0) = np.DotProduct(KK2[5]*np);
	T2(0,2) = T2(2,0) = np.DotProduct(KK2[4]*np);
	T2(1,2) = T2(2,1) = np.DotProduct(KK2[3]*np);
	
	std::cout << "T1: " << std::endl;
	T1.Print();
	
	std::cout << "T2: " << std::endl;
	T2.Print();
	
	// 0  5  4
	// 5T 1  3
	// 4T 3T 2
	//
	// 0 1 2
	// 3 4 5
	// 6 7 8
	//
	// KK1n: [0 5 4 5T 1 3 4T 3T 2]
	KK1n[0] = KK1[0]*np;
	KK1n[1] = KK1[5]*np;
	KK1n[2] = KK1[4]*np;
	KK1n[3] = KK1[5].Transpose()*np;
	KK1n[4] = KK1[1]*np;
	KK1n[5] = KK1[3]*np;
	KK1n[6] = KK1[4].Transpose()*np;
	KK1n[7] = KK1[3].Transpose()*np;
	KK1n[8] = KK1[2]*np;
	
	KK2n[0] = KK2[0]*np;
	KK2n[1] = KK2[5]*np;
	KK2n[2] = KK2[4]*np;
	KK2n[3] = KK2[5].Transpose()*np;
	KK2n[4] = KK2[1]*np;
	KK2n[5] = KK2[3]*np;
	KK2n[6] = KK2[4].Transpose()*np;
	KK2n[7] = KK2[3].Transpose()*np;
	KK2n[8] = KK2[2]*np;
	
	rMatrix KK1U, KK1V, KK1W;
	rMatrix KK2U, KK2V, KK2W;
	rMatrix KKK;
	
	
	KK1U = KK1n[0].Transpose().ConcatRows(KK1n[3].Transpose()).ConcatRows(KK1n[6].Transpose());
	KK1V = KK1n[1].Transpose().ConcatRows(KK1n[4].Transpose()).ConcatRows(KK1n[7].Transpose());
	KK1W = KK1n[2].Transpose().ConcatRows(KK1n[5].Transpose()).ConcatRows(KK1n[8].Transpose());
	
	std::cout << "KK1U:" << std::endl;
	KK1U.Print();
	std::cout << "KK1V:" << std::endl;
	KK1V.Print();
	std::cout << "KK1W:" << std::endl;
	KK1W.Print();
	
	KK2U = KK2n[0].Transpose().ConcatRows(KK2n[3].Transpose()).ConcatRows(KK2n[6].Transpose());
	KK2V = KK2n[1].Transpose().ConcatRows(KK2n[4].Transpose()).ConcatRows(KK2n[7].Transpose());
	KK2W = KK2n[2].Transpose().ConcatRows(KK2n[5].Transpose()).ConcatRows(KK2n[8].Transpose());
	
	std::cout << "KK2U:" << std::endl;
	KK2U.Print();
	std::cout << "KK2V:" << std::endl;
	KK2V.Print();
	std::cout << "KK2W:" << std::endl;
	KK2W.Print();
	
	//Q = rMatrix::Unit(3);
	std::cout << "Q:" << std::endl;
	Q.Print();
	

	
	rMatrix KKU = KK1U - (Q(0,0)*KK2U + Q(1,0)*KK2V + Q(2,0)*KK2W);
	rMatrix KKV = KK1V - (Q(0,1)*KK2U + Q(1,1)*KK2V + Q(2,1)*KK2W);
	rMatrix KKW = KK1W - (Q(0,2)*KK2U + Q(1,2)*KK2V + Q(2,2)*KK2W);
	
	std::cout << "KK1U - KK2U:" << std::endl;
	(KK1U-KK2U).Print();
	
	std::cout << "KK1V - KK2V:" << std::endl;
	(KK1V-KK2V).Print();
	
	std::cout << "KK1W - KK2W:" << std::endl;
	(KK1W-KK2W).Print();
	
	
	std::cout << "KK1U - q_11*KK2U - q_21*KK2V - q_31*KK2W:" << std::endl;
	KKU.Print();
	std::cout << "KK2V - q_12*KK2U - q_22*KK2V - q_32*KK2W:" << std::endl;
	KKV.Print();
	std::cout << "KK2W - q_13*KK2U - q_23*KK2V - q_33*KK2W:" << std::endl;
	KKW.Print();
	
	
	KKK = KKU.ConcatCols(KKV).ConcatCols(KKW);
	std::cout << "KKK:" << std::endl;
	KKK.Print();
	
	std::cout << "G1:" << std::endl;
	G1.Print();
	
	std::cout << "vector(G1):" << std::endl;
	G1.Repack(9,1).Print();
	
	
	
	rMatrix uG1 = G1.SubMatrix(0,0,1,3).Transpose();
	rMatrix vG1 = G1.SubMatrix(1,0,2,3).Transpose();
	rMatrix wG1 = G1.SubMatrix(2,0,3,3).Transpose();
	
	std::cout << "uG1:" << std::endl;
	uG1.Print();
	std::cout << "vG1:" << std::endl;
	vG1.Print();
	std::cout << "wG1:" << std::endl;
	wG1.Print();
	
	std::cout << "NP*C1*sigma:" << std::endl;
	(NP*C1r*sigma).Print();
	
	std::cout << "KK1U*uG1+KK1V*vG1+KK1W*wG1" << std::endl;
	(KK1U*uG1+KK1V*vG1+KK1W*wG1).Print();
	
	std::cout << "NP*C2*sigma:" << std::endl;
	(NP*C2r*sigma).Print();
	
	std::cout << "KK2U*uG1+KK2V*vG1+KK2W*wG1" << std::endl;
	(KK2U*uG1+KK2V*vG1+KK2W*wG1).Print();
	
	
	std::cout << "KKK*vector(G1):" << std::endl;
	(KKK*G1.Repack(9,1)).Print();
	
	std::cout << "KKU*uG1+KKV*vG1+KKW*wG1:" << std::endl;
	(KKU*uG1+KKV*vG1+KKW*wG1).Print();
	
	std::cout << "NP*(C1-C2)*sigma" << std::endl;
	(NP*(C1r-C2r)*sigma).Print();
	
	std::cout << "KK1U*uG1+KK1V*vG1+KK1W*wG1-KK2U*uG1+KK2V*vG1+KK2W*wG1" << std::endl;
	(KK1U*uG1+KK1V*vG1+KK1W*wG1-(KK2U*uG1+KK2V*vG1+KK2W*wG1)).Print();
	
	std::cout << "(KK1U-KK2U)*uG1+(KK1V-KK2V)*vG1+(KK1W-KK2W)*wG1" << std::endl;
	((KK1U-KK2U)*uG1+(KK1V-KK2V)*vG1+(KK1W-KK2W)*wG1).Print();
	
	
	rMatrix x(3,1), u2(3,1), u1(3,1), uI(3,1);
	x(0,0) = 1;
	x(1,0) = 1;
	x(2,0) = 1;
	double r1 = np.DotProduct(xp);
	double r2 = np.DotProduct(x-xp);
	std::cout << "r1: " << r1 << std::endl;
	std::cout << "r2: " << r2 << std::endl;
	u1 = G1*x+d1;
	u2 = G2*x+d2;
	
	std::cout << "u1:" << std::endl;
	u1.Print();
	
	std::cout << "pressure contribution:" << std::endl;
	(r2*T2.Invert()*(B1*p1/K1 - B2*p2/K2)*np).Print();
	
	std::cout << "flux contribution:" << std::endl;
	(r2*T2.Invert()*(KKK*G1.Repack(9,1))).Print();
	
	std::cout << "pressure+flux contribution:" << std::endl;
	(r2*T2.Invert()*(KKK*G1.Repack(9,1) + (B1*p1/K1 - B2*p2/K2)*np)).Print();
	
	std::cout << "interpolation contribution:" << std::endl;
	(Q*u1).Print();
	
	
	
	std::cout << "TB^{-1}" << std::endl;
	B.Invert().Print();
	
	std::cout << "TBr^{-1}" << std::endl;
	T2.Invert().Print();
	
	uI = Q*u1 + r2*T2.Invert()*(KKK*G1.Repack(9,1) + (B1*p1/K1 - B2*p2/K2)*np);
	
	//uI = u1 + r2*T2.Invert()*(NP*(C1r-C2r)*sigma + NP*(K1*w1*p1 - K2*w2*p2));
	//uI = u1 + r2*B.Invert()*(NP*(C1-C2)*sigma);
	//uI = u1 + r2*B.Invert()*(KKK*G1.Repack(9,1));
	
	std::cout << "u2:" << std::endl;
	u2.Print();
	std::cout << "uI:" << std::endl;
	uI.Print();
	
	std::cout << "np^T*B2^{-1}*B1*(d1+G1*xp)" << std::endl;
	(np.Transpose()*Q*(d1+G1*xp)).Print();
	std::cout << "np^T*(d2+G2*xp)" << std::endl;
	(np.Transpose()*(d2+G2*xp)).Print();
	
	C1 = C1.Invert();
	C2 = C2.Invert();
	
#endif //TEST_INTRP
	
	
	// u1 + G1 (x2-x1)

	/*
	double Gv2[9] =
	{
		
		9.974114e+00, 4.092234e+01, 8.946117e+01,
		4.251629e+01, -1.245114e+01, -1.222557e+01,
		7.001355e+01, 1.000407e+02, 4.202033e+01
	};
	
	double dv2[3] =
	{
		2.132119e+01,
		9.958022e+01,
		4.646273e+01
	};
	
	G2 = rMatrix(Gv2,3,3);
	d2 = rMatrix(dv2,3,1);
	*/
	
	
	
	for(Mesh::iteratorElement it = m->BeginElement(CELL|FACE); it != m->EndElement(); ++it)
	{
		it->Barycenter(c);
		
		int zone = 0;
		rMatrix xv(c,3,1);
		if( np.DotProduct(xv-xp) > 0 )
			zone = 1;
		
		
		unknown x(c[0],0);
		unknown y(c[1],1);
		unknown z(c[2],2);
		vMatrix xyz(3,1);
		hMatrix sol(3,1);
		vMatrix epsilon(6,1), sigma(6,1);
		rMatrix N(3,6), sigma_value(6,1), flux(3,1);
		
		xyz(0,0) = x;
		xyz(1,0) = y;
		xyz(2,0) = z;
		
		if( zone == 0 )
			sol = G1*xyz + d1;
		else
			sol = G2*xyz + d2;
		
		epsilon(0,0) = sol(0,0).GetVariable(0); //u_x
		epsilon(1,0) = sol(1,0).GetVariable(1); //v_y
		epsilon(2,0) = sol(2,0).GetVariable(2); //w_z
		epsilon(3,0) = sol(1,0).GetVariable(2) + sol(2,0).GetVariable(1); //v_z + w_y
		epsilon(4,0) = sol(0,0).GetVariable(2) + sol(2,0).GetVariable(0); //u_z + w_x
		epsilon(5,0) = sol(0,0).GetVariable(1) + sol(1,0).GetVariable(0); //u_y + v_x
		
		
		if( zone == 0 )
			sigma = C1*epsilon;
		else
			sigma = C2*epsilon;
		
		
		
		
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
			if( zone == 0 )
			{
				perm[ 0] = C1(0,0);
				perm[ 1] = C1(0,1);
				perm[ 2] = C1(0,2);
				perm[ 3] = C1(0,3);
				perm[ 4] = C1(0,4);
				perm[ 5] = C1(0,5);
				perm[ 6] = C1(1,1);
				perm[ 7] = C1(1,2);
				perm[ 8] = C1(1,3);
				perm[ 9] = C1(1,4);
				perm[10] = C1(1,5);
				perm[11] = C1(2,2);
				perm[12] = C1(2,3);
				perm[13] = C1(2,4);
				perm[14] = C1(2,5);
				perm[15] = C1(3,3);
				perm[16] = C1(3,4);
				perm[17] = C1(3,5);
				perm[18] = C1(4,4);
				perm[19] = C1(4,5);
				perm[20] = C1(5,5);
			}
			else
			{
				perm[ 0] = C2(0,0);
				perm[ 1] = C2(0,1);
				perm[ 2] = C2(0,2);
				perm[ 3] = C2(0,3);
				perm[ 4] = C2(0,4);
				perm[ 5] = C2(0,5);
				perm[ 6] = C2(1,1);
				perm[ 7] = C2(1,2);
				perm[ 8] = C2(1,3);
				perm[ 9] = C2(1,4);
				perm[10] = C2(1,5);
				perm[11] = C2(2,2);
				perm[12] = C2(2,3);
				perm[13] = C2(2,4);
				perm[14] = C2(2,5);
				perm[15] = C2(3,3);
				perm[16] = C2(3,4);
				perm[17] = C2(3,5);
				perm[18] = C2(4,4);
				perm[19] = C2(4,5);
				perm[20] = C2(5,5);
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


