

#include "inmost.h"

using namespace INMOST;
std::string problem_name = "randmap3d_permeability";
const double pi = 3.1415926535897932384626433832795;

#define ind(r,c,b) (((r)*N + (c))*N+(b))
const double noind = std::numeric_limits<double>::max();

void rand3d(double* arr, int N, int Nl, int Nr, int Nb, int Nt, int Nn, int Nf, double t)
{
	//std::cout << "Nl:Nr " << Nl <<":" << Nr << " Nb:Nt " << Nb << ":" << Nt << std::endl;
	if (Nr - Nl < 2 && Nt - Nb < 2 && Nf - Nn < 2)
	{
		//std::cout << "exit" << std::endl;
		return;
	}
	//const double t = 0.15;
	int Nk = (Nb + Nt) / 2;
	int Nm = (Nl + Nr) / 2;
	int Nq = (Nf + Nn) / 2;
	//std::cout << "Nk " << Nk << " Nm " << Nm << std::endl;
	double lbn = arr[ind(Nb, Nl, Nn)];
	double rbn = arr[ind(Nb, Nr, Nn)];
	double ltn = arr[ind(Nt, Nl, Nn)];
	double rtn = arr[ind(Nt, Nr, Nn)];
	double lbf = arr[ind(Nb, Nl, Nf)];
	double rbf = arr[ind(Nb, Nr, Nf)];
	double ltf = arr[ind(Nt, Nl, Nf)];
	double rtf = arr[ind(Nt, Nr, Nf)];
	if (lbn != lbn || rbn != rbn || ltn != ltn || rtn != rtn) throw - 1;
	if (lbf != lbf || rbf != rbf || ltf != ltf || rtf != rtf) throw - 1;
	if (arr[ind(Nk, Nl, Nn)] == noind) arr[ind(Nk, Nl, Nn)] = 0.5 * (lbn + ltn) + (2 * (rand() * 1.0 / RAND_MAX) - 1) * t;
	if (arr[ind(Nk, Nl, Nf)] == noind) arr[ind(Nk, Nl, Nf)] = 0.5 * (lbf + ltf) + (2 * (rand() * 1.0 / RAND_MAX) - 1) * t;
	if (arr[ind(Nk, Nr, Nn)] == noind) arr[ind(Nk, Nr, Nn)] = 0.5 * (rbn + rtn) + (2 * (rand() * 1.0 / RAND_MAX) - 1) * t;
	if (arr[ind(Nk, Nr, Nf)] == noind) arr[ind(Nk, Nr, Nf)] = 0.5 * (rbf + rtf) + (2 * (rand() * 1.0 / RAND_MAX) - 1) * t;
	if (arr[ind(Nb, Nm, Nn)] == noind) arr[ind(Nb, Nm, Nn)] = 0.5 * (lbn + rbn) + (2 * (rand() * 1.0 / RAND_MAX) - 1) * t;
	if (arr[ind(Nb, Nm, Nf)] == noind) arr[ind(Nb, Nm, Nf)] = 0.5 * (lbf + rbf) + (2 * (rand() * 1.0 / RAND_MAX) - 1) * t;
	if (arr[ind(Nt, Nm, Nn)] == noind) arr[ind(Nt, Nm, Nn)] = 0.5 * (ltn + rtn) + (2 * (rand() * 1.0 / RAND_MAX) - 1) * t;
	if (arr[ind(Nt, Nm, Nf)] == noind) arr[ind(Nt, Nm, Nf)] = 0.5 * (ltf + rtf) + (2 * (rand() * 1.0 / RAND_MAX) - 1) * t;
	if (arr[ind(Nb, Nl, Nq)] == noind) arr[ind(Nb, Nl, Nq)] = 0.5 * (lbn + lbf) + (2 * (rand() * 1.0 / RAND_MAX) - 1) * t;
	if (arr[ind(Nb, Nr, Nq)] == noind) arr[ind(Nb, Nr, Nq)] = 0.5 * (rbn + rbf) + (2 * (rand() * 1.0 / RAND_MAX) - 1) * t;
	if (arr[ind(Nt, Nl, Nq)] == noind) arr[ind(Nt, Nl, Nq)] = 0.5 * (ltn + ltf) + (2 * (rand() * 1.0 / RAND_MAX) - 1) * t;
	if (arr[ind(Nt, Nr, Nq)] == noind) arr[ind(Nt, Nr, Nq)] = 0.5 * (rtn + rtf) + (2 * (rand() * 1.0 / RAND_MAX) - 1) * t;
	if (arr[ind(Nb, Nm, Nq)] == noind) arr[ind(Nb, Nm, Nq)] = 0.25 * (lbn + lbf + rbn + rbf) + (2 * (rand() * 1.0 / RAND_MAX) - 1) * t;
	if (arr[ind(Nt, Nm, Nq)] == noind) arr[ind(Nt, Nm, Nq)] = 0.25 * (ltn + ltf + rtn + rtf) + (2 * (rand() * 1.0 / RAND_MAX) - 1) * t;
	if (arr[ind(Nk, Nl, Nq)] == noind) arr[ind(Nk, Nl, Nq)] = 0.25 * (lbn + lbf + ltn + ltf) + (2 * (rand() * 1.0 / RAND_MAX) - 1) * t;
	if (arr[ind(Nk, Nr, Nq)] == noind) arr[ind(Nk, Nr, Nq)] = 0.25 * (rtn + rtf + rbn + rbf) + (2 * (rand() * 1.0 / RAND_MAX) - 1) * t;
	if (arr[ind(Nk, Nm, Nn)] == noind) arr[ind(Nk, Nm, Nn)] = 0.25 * (lbn + rbn + ltn + rtn) + (2 * (rand() * 1.0 / RAND_MAX) - 1) * t;
	if (arr[ind(Nk, Nm, Nf)] == noind) arr[ind(Nk, Nm, Nf)] = 0.25 * (lbf + rbf + ltf + rtf) + (2 * (rand() * 1.0 / RAND_MAX) - 1) * t;
	arr[ind(Nk, Nm, Nq)] = 0.125 * (lbn + rbn + ltn + rtn + lbf + rbf + ltf + rtf) + (2 * (rand() * 1.0 / RAND_MAX) - 1) * t;
	rand3d(arr, N, Nl, Nm, Nb, Nk, Nn, Nq, t * 0.5); 
	rand3d(arr, N, Nm, Nr, Nb, Nk, Nn, Nq, t * 0.5);
	rand3d(arr, N, Nl, Nm, Nk, Nt, Nn, Nq, t * 0.5);
	rand3d(arr, N, Nm, Nr, Nk, Nt, Nn, Nq, t * 0.5);
	rand3d(arr, N, Nl, Nm, Nb, Nk, Nq, Nf, t * 0.5);
	rand3d(arr, N, Nm, Nr, Nb, Nk, Nq, Nf, t * 0.5);
	rand3d(arr, N, Nl, Nm, Nk, Nt, Nq, Nf, t * 0.5);
	rand3d(arr, N, Nm, Nr, Nk, Nt, Nq, Nf, t * 0.5);
}

void init3d(double* arr, int N, double mint, double maxt)
{
	for (int k = 0; k < N * N * N; ++k) arr[k] = noind;
	int L = 0, R = N - 1;
	arr[ind(L, L, L)] = mint + (rand() * 1.0 / RAND_MAX) * (maxt - mint);
	arr[ind(L, R, L)] = mint + (rand() * 1.0 / RAND_MAX) * (maxt - mint);
	arr[ind(R, L, L)] = mint + (rand() * 1.0 / RAND_MAX) * (maxt - mint);
	arr[ind(R, R, L)] = mint + (rand() * 1.0 / RAND_MAX) * (maxt - mint);
	arr[ind(L, L, R)] = mint + (rand() * 1.0 / RAND_MAX) * (maxt - mint);
	arr[ind(L, R, R)] = mint + (rand() * 1.0 / RAND_MAX) * (maxt - mint);
	arr[ind(R, L, R)] = mint + (rand() * 1.0 / RAND_MAX) * (maxt - mint);
	arr[ind(R, R, R)] = mint + (rand() * 1.0 / RAND_MAX) * (maxt - mint);
}

double intrp3d(const double* arr, int N, double x, double y, double z)
{
	int n = (int)ceil(x * (N - 1));
	int m = (int)ceil(y * (N - 1));
	int k = (int)ceil(z * (N - 1));
	if (n == 0) n = 1;
	if (m == 0) m = 1;
	if (k == 0) k = 1;
	double dh = 1.0 / (double)(N - 1);
	double kx = (x - (n - 1) * dh) / dh;
	double ky = (y - (m - 1) * dh) / dh;
	double kz = (z - (k - 1) * dh) / dh;
	//if( kx < 0 || kx > 1 ) std::cout << "bad kx: " << kx << " x is " << x << " n is " << n << " dh is " << dh << " N is " << N << std::endl;
	//if( ky < 0 || ky > 1 ) std::cout << "bad ky: " << ky << " y is " << y << " m is " << m << " dh is " << dh << " N is " << N << std::endl;
	double lbn = arr[ind(m - 1, n - 1, k - 1)];
	double rbn = arr[ind(m - 1, n + 0, k - 1)];
	double ltn = arr[ind(m + 0, n - 1, k - 1)];
	double rtn = arr[ind(m + 0, n + 0, k - 1)];
	double lbf = arr[ind(m - 1, n - 1, k + 0)];
	double rbf = arr[ind(m - 1, n + 0, k + 0)];
	double ltf = arr[ind(m + 0, n - 1, k + 0)];
	double rtf = arr[ind(m + 0, n + 0, k + 0)];
	if (lbn != lbn || rbn != rbn || ltn != ltn || rtn != rtn) throw - 1;
	if (lbf != lbf || rbf != rbf || ltf != ltf || rtf != rtf) throw - 1;
	return (1 - kz) * ((1 - ky) * (lbn * (1 - kx) + rbn * kx) + ky * (ltn * (1 - kx) + rtn * kx))
		+        kz * ((1 - ky) * (lbf * (1 - kx) + rbf * kx) + ky * (ltf * (1 - kx) + rtf * kx));
}


int main(int argc, char ** argv)
{
	if( argc < 3 )
	{
		std::cout << "Usage: " << argv[0] << " mesh mesh_out [K=1000] [dK=10] [L=4] [N=100]" << std::endl;
		return 0;
	}
	
	Mesh * m = new Mesh;
	m->SetFileOption("VERBOSITY","2");
	try{m->Load(argv[1]);} catch(...) { std::cout << "Cannot load the mesh " << argv[1] << std::endl; return -1;}


	double max[3] = {-1.0e20, -1.0e20, -1.0e20}, min[3] = {1.0e20,1.0e20,1.0e20};
	for(Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); ++it)
	{
		Storage::real_array c = it->Coords();
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

	TagReal perm = m->CreateTag("PERM",DATA_REAL,CELL,NONE,1);
	
	int N = 100;
	Storage::real K = 1000.0, dK = 10.0, L = 4.0;

	
	if (argc > 3) K = atof(argv[3]);
	if (argc > 4) dK = atof(argv[4]);
	if (argc > 5) L = atof(argv[5]);
	if (argc > 6) N = atoi(argv[6]);

	double* map = new double[N * N * N];
	std::fill(map, map + N * N * N, 0.0);
	init3d(map, N, 0.0, 1.0);
	rand3d(map, N, 0, N - 1, 0, N - 1, 0, N - 1, 0.5);

    
    {
        Storage::bulk_array name = m->self()->BulkArray(m->CreateTag("PROBLEMNAME",DATA_BULK,MESH,NONE));
        name.replace(name.begin(),name.end(),problem_name.begin(),problem_name.end());
    }


	for (Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); ++it)
	{
		Storage::real_array c = it->Coords();
		for (int k = 0; k < c.size(); ++k)
			c[k] = (c[k] - min[k]) / (max[k] - min[k]);
		c[2] *= L;
	}

	TagReal tag_C = m->CreateTag("C", DATA_REAL, CELL, NONE, 1);

	std::cout << "Setting randmap permeability" << std::endl;
	const Storage::real h = 1.0 / 25.0, hz = 1.0 / (L * 25.0);
	for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
	{
		Storage::real c[3];
		it->Centroid(c);
		Storage::real x = c[0], y = c[1], z = c[2] / L;

		if (z * L < 1.0)
			tag_C[*it] = 1.0;
		else tag_C[*it] = 0.0;

		x = floor(x / h) * h + 0.5 * h;
		y = floor(y / h) * h + 0.5 * h;
		z = floor(z / hz) * hz + 0.5 * hz;

		double pert = intrp3d(map, N, x, y, z) * 90.0;
		perm[*it] = K + pert * dK;
	}

	TagRealArray tag_BCP = m->CreateTag("BOUNDARY_CONDITION_PRESSURE", DATA_REAL, FACE, FACE, 3);
	TagRealArray tag_BCC = m->CreateTag("BOUNDARY_CONDITION_FRACTION", DATA_REAL, FACE, FACE, 3);
	for (Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it)
	{
		Storage::real nrm[3];
		it->UnitNormal(nrm);
		if (fabs(nrm[2] + 1.0) < 1.0e-3) //back side
		{
			//pressure equals 40
			tag_BCP[*it][0] = 1.0;
			tag_BCP[*it][1] = 0.0;
			tag_BCP[*it][0] = 40.0;
			//fraction equals 1
			tag_BCC[*it][0] = 1.0;
			tag_BCC[*it][1] = 0.0;
			tag_BCC[*it][0] = 1.0;
		}
		if (fabs(nrm[2] - 1.0) < 1.0e-3) //front side
		{
			//pressure equals 0
			tag_BCP[*it][0] = 1.0;
			tag_BCP[*it][1] = 0.0;
			tag_BCP[*it][0] = 0.0;
		}
	}

	std::cout << "Saving output to " << argv[2] << std::endl;
	m->Save(argv[2]);
	delete m;
}

