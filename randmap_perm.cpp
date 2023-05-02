

#include "inmost.h"

using namespace INMOST;
std::string problem_name = "randmap_permeability";
const double pi = 3.1415926535897932384626433832795;

#define ind(r,c) ((r)*N + (c))
const double noind = std::numeric_limits<double>::max();

void rand2d(double* arr, int N, int Nl, int Nr, int Nb, int Nt, double t)
{
	//std::cout << "Nl:Nr " << Nl <<":" << Nr << " Nb:Nt " << Nb << ":" << Nt << std::endl;
	if (Nr - Nl < 2 && Nt - Nb < 2)
	{
		//std::cout << "exit" << std::endl;
		return;
	}
	//const double t = 0.15;
	int Nk = (Nb + Nt) / 2;
	int Nm = (Nl + Nr) / 2;
	//std::cout << "Nk " << Nk << " Nm " << Nm << std::endl;
	double lb = arr[ind(Nb, Nl)];
	double rb = arr[ind(Nb, Nr)];
	double lt = arr[ind(Nt, Nl)];
	double rt = arr[ind(Nt, Nr)];
	if (lb != lb || rb != rb || lt != lt || rt != rt) throw - 1;
	if (arr[ind(Nk, Nl)] == noind) arr[ind(Nk, Nl)] = 0.5 * (lb + lt) + (2 * (rand() * 1.0 / RAND_MAX) - 1) * t;
	if (arr[ind(Nk, Nr)] == noind) arr[ind(Nk, Nr)] = 0.5 * (rb + rt) + (2 * (rand() * 1.0 / RAND_MAX) - 1) * t;
	if (arr[ind(Nb, Nm)] == noind) arr[ind(Nb, Nm)] = 0.5 * (lb + rb) + (2 * (rand() * 1.0 / RAND_MAX) - 1) * t;
	if (arr[ind(Nt, Nm)] == noind) arr[ind(Nt, Nm)] = 0.5 * (lt + rt) + (2 * (rand() * 1.0 / RAND_MAX) - 1) * t;
	arr[ind(Nk, Nm)] = 0.25 * (lb + rb + lt + rt) + (2 * (rand() * 1.0 / RAND_MAX) - 1) * t;
	rand2d(arr, N, Nl, Nm, Nb, Nk, t * 0.5); rand2d(arr, N, Nm, Nr, Nb, Nk, t * 0.5);
	rand2d(arr, N, Nl, Nm, Nk, Nt, t * 0.5); rand2d(arr, N, Nm, Nr, Nk, Nt, t * 0.5);
}

void init2d(double* arr, int N, double mint, double maxt)
{
	for (int k = 0; k < N * N; ++k) arr[k] = noind;
	arr[ind(0, 0)] = mint + (rand() * 1.0 / RAND_MAX) * (maxt - mint);
	arr[ind(0, N - 1)] = mint + (rand() * 1.0 / RAND_MAX) * (maxt - mint);
	arr[ind(N - 1, 0)] = mint + (rand() * 1.0 / RAND_MAX) * (maxt - mint);
	arr[ind(N - 1, N - 1)] = mint + (rand() * 1.0 / RAND_MAX) * (maxt - mint);
}

double intrp2d(const double* arr, int N, double x, double y)
{
	int n = (int)ceil(x * (N - 1));
	int m = (int)ceil(y * (N - 1));
	if (n == 0) n = 1;
	if (m == 0) m = 1;
	double dh = 1.0 / (double)(N - 1);
	double kx = (x - (n - 1) * dh) / dh;
	double ky = (y - (m - 1) * dh) / dh;
	//if( kx < 0 || kx > 1 ) std::cout << "bad kx: " << kx << " x is " << x << " n is " << n << " dh is " << dh << " N is " << N << std::endl;
	//if( ky < 0 || ky > 1 ) std::cout << "bad ky: " << ky << " y is " << y << " m is " << m << " dh is " << dh << " N is " << N << std::endl;
	double lb = arr[ind(m - 1, n - 1)];
	double rb = arr[ind(m - 1, n + 0)];
	double lt = arr[ind(m + 0, n - 1)];
	double rt = arr[ind(m + 0, n + 0)];
	if (lb != lb || rb != rb || lt != lt || rt != rt) throw - 1;
	return (1 - ky) * (lb * (1 - kx) + rb * kx) + ky * (lt * (1 - kx) + rt * kx);
}


int main(int argc, char ** argv)
{
	if( argc < 3 )
	{
		std::cout << "Usage: " << argv[0] << " mesh mesh_out [Kx=1000] [Ky=10] [N=100]" << std::endl;
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

	int N = 100;
	Storage::real Kx = 1000.0, Ky = 10.0, Kz = 1.0;

	
	if (argc > 3) Kx = atof(argv[3]);
	if (argc > 4) Ky = atof(argv[4]);
	if (argc > 5) N = atoi(argv[5]);

	double* map = new double[N * N];
	std::fill(map, map + N * N, 0.0);
	init2d(map, N, 0.0, 1.0);
	rand2d(map, N, 0, N - 1, 0, N - 1, 0.5);

    
    {
        Storage::bulk_array name = m->self()->BulkArray(m->CreateTag("PROBLEMNAME",DATA_BULK,MESH,NONE));
        name.replace(name.begin(),name.end(),problem_name.begin(),problem_name.end());
    }


	std::cout << "Setting randmap permeability" << std::endl;
	const Storage::real h = 1.0 / 25.0;
	for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
	{
		Storage::real_array perm = tensor[*it];
		it->Centroid(c);
		Storage::real a = (c[0]-min[0])/(max[0]-min[0]);
		Storage::real b  = (c[1]-min[1])/(max[1]-min[1]);

		a = floor(a / h) * h + 0.5 * h;
		b = floor(b / h) * h + 0.5 * h;

		double angle = intrp2d(map, N, a, b) * 90.0;
		perm[0] = Kx + angle * Ky;
		perm[3] = Kx + angle * Ky;
		perm[5] = Kz;
		
		/*
		Storage::real acos = cos(angle / 180.0 * pi);
		Storage::real asin = sin(angle / 180.0 * pi);
		perm[0] = Kx * acos * acos + Ky * asin * asin;
		perm[1] = (Kx - Ky) * acos * asin;
		perm[2] = 0.0;
		perm[3] = Ky * acos * acos + Kx * asin * asin;
		perm[4] = 0.0;
		perm[5] = Kz;
		*/
	}
	std::cout << "Saving output to " << argv[2] << std::endl;
	m->Save(argv[2]);
	delete m;
}

