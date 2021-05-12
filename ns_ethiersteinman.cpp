#include "inmost.h"

using namespace INMOST;
//setup on mesh:
//
// FORCE - 3 entries
//
// BOUNDARY_CONDITION_VELOCITY - 7 entries
// r = (a5,a6,a7)
// n.(a1*u + a2*t) = n.r
// (I - nn.)(a3*u + a4*t) = (I-nn.)r
//
// BOUNDARY_CONDITION_PRESSURE - 1 entry
//
// REFERENCE_VELOCITY - 3 entries
//
// REFERENCE_PRESSURE - 1 entry
//

typedef Storage::real real;


const real pi = 3.1415926535897932384626433832795;


real refUmean(real cnt[3], real t)
{
	real x = cnt[0], y = cnt[1], z = cnt[2];
	const real a = pi / 4.0;
	const real d = pi / 2.0;
	return -a * (exp(a * x) * sin(a * y + d * z) + exp(a * z) * cos(a * x + d * y));// *exp(-nu_global * d * d * t);
}

real refVmean(real cnt[3], real t)
{
	real x = cnt[0], y = cnt[1], z = cnt[2];
	const real a = pi / 4.0;
	const real d = pi / 2.0;
	return -a * (exp(a * y) * sin(a * z + d * x) + exp(a * x) * cos(a * y + d * z));// *exp(-nu_global * d * d * t);
}

real refWmean(real cnt[3], real t)
{
	real x = cnt[0], y = cnt[1], z = cnt[2];
	const real a = pi / 4.0;
	const real d = pi / 2.0;
	return -a * (exp(a * z) * sin(a * x + d * y) + exp(a * y) * cos(a * z + d * x));// *exp(-nu_global * d * d * t);
}


real refPmean(real cnt[3], real t)
{
	real x = cnt[0], y = cnt[1], z = cnt[2];
	const real a = pi / 4.0;
	const real d = pi / 2.0;
	return -a * a / 2.0 * (exp(2 * a * x) + exp(2 * a * y) + exp(2 * a * z) +
		2.0 * sin(a * x + d * y) * cos(a * z + d * x) * exp(a * (y + z)) +
		2.0 * sin(a * y + d * z) * cos(a * x + d * y) * exp(a * (z + x)) +
		2.0 * sin(a * z + d * x) * cos(a * y + d * z) * exp(a * (x + y)));// *exp(-2 * nu * d * d * t);
}



variable refP(real vx, real vy, real vz)//, real vt, real nu)
{
	unknown x(vx, 0);
	unknown y(vy, 1);
	unknown z(vz, 2);
	//unknown t(vt, 3);
	const real a = pi / 4.0;
	const real d = pi / 2.0;
	return -a * a / 2.0 * (exp(2 * a * x) + exp(2 * a * y) + exp(2 * a * z) +
		2.0 * sin(a * x + d * y) * cos(a * z + d * x) * exp(a * (y + z)) +
		2.0 * sin(a * y + d * z) * cos(a * x + d * y) * exp(a * (z + x)) +
		2.0 * sin(a * z + d * x) * cos(a * y + d * z) * exp(a * (x + y)));// *exp(-2 * nu * d * d * t);
}

int main(int argc, char ** argv)
{
	
	if( argc < 2 )
	{
		std::cout << "Usage: " << argv[0] << " mesh mesh_out=grid_out.pmf skew=0 setU=0 setP=0 mean=1 setGP0=1" << std::endl;
		return 0;
	}
	
	Mesh * m = new Mesh;
	m->SetFileOption("VERBOSITY","2");
	try
	{
		m->Load(argv[1]);
		Mesh::GeomParam table;
		table[ORIENTATION] = FACE;
		m->PrepareGeometricData(table);
	}
	catch(...)
	{
		std::cout << "Cannot load the mesh " << argv[1] << std::endl;
		return -1;
	}
	
	std::string fout = "grid_out.pmf";
	int setP = 0, setU = 0, mean = 1, setGP0 = 1;
	double skew = 0;
	if (argc > 2) fout = std::string(argv[2]);
	if (argc > 3) skew = atof(argv[3]);
	if (argc > 4) setU = atoi(argv[4]);
	if (argc > 5) setP = atoi(argv[5]);
	if (argc > 6) mean = atoi(argv[6]);
	if (argc > 7) setGP0 = atoi(argv[7]);

	
	if( skew > 0)
	{
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
		for(Mesh::iteratorNode n = m->BeginNode(); n != m->EndNode(); ++n)
		{
			Storage::real_array c = n->Coords();
			for(int k = 0; k < 3; ++k)
			{
				double t = (c[k]-cmin[k])/(cmax[k]-cmin[k]);
				double a = 0.5 + 0.5*(t-0.5)/sqrt(pow(t-0.5,2)+pow(skew,2))*sqrt(0.25+pow(skew,2))*2;
				c[k] = cmin[k] + (cmax[k]-cmin[k])*a;
			}
		}
	}
	
	//TagRealArray force = m->CreateTag("FORCE",DATA_REAL,CELL,NONE,3);
	TagRealArray bc = m->CreateTag("BOUNDARY_CONDITION_VELOCITY",DATA_REAL,FACE,FACE,7);
	//TagReal bcp = m->CreateTag("BOUNDARY_CONDITION_PRESSURE",DATA_REAL,FACE,FACE,1);
	TagRealArray uvw = m->CreateTag("UVW",DATA_REAL,CELL,NONE,3);
	TagReal p = m->CreateTag("P",DATA_REAL,CELL,NONE,1);
	TagReal gp0;
	if( setP )
		gp0 = m->CreateTag("GP0", DATA_REAL, FACE, FACE, 1);
	
	rMatrix n(3,1),x(3,1), gp(3,1);

	if (mean) std::cout << "Calculating mean values" << std::endl;
	if (setU) std::cout << "Set initial U" << std::endl;
	if (setP) std::cout << "Set initial P" << std::endl;
	
	real intP = 0, intV = 0, V;
	for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
	{
		if (setU)
		{
			if (mean)
			{
				uvw[*it][0] = it->Mean(refUmean, 0);
				uvw[*it][1] = it->Mean(refVmean, 0);
				uvw[*it][2] = it->Mean(refWmean, 0);
			}
			else
			{
				it->Barycenter(x.data());
				uvw[*it][0] = refUmean(x.data(), 0);
				uvw[*it][1] = refVmean(x.data(), 0);
				uvw[*it][2] = refWmean(x.data(), 0);
			}
		} 
		else uvw(*it, 3, 1).Zero();
		if (setP)
		{
			if (mean)
				p[*it] = it->Mean(refPmean, 0);
			else
			{
				it->Barycenter(x.data());
				p[*it] = refPmean(x.data(),0);
			}
			V = it->Volume();
			intP += p[*it] * V;
			intV += V;
		}
		else p[*it] = 0;
	}
	if (setP)
	{
		if (intV) intP /= intV;
		std::cout << "pressure integral: " << intP << std::endl;
		if (intP)
		{
			std::cout << "shift pressure" << std::endl;
			for (Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
				p[*it] -= intP;
		}
	}
	
	if (setGP0) std::cout << "set normal pressure gradient on boundary" << std::endl;
	
	for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it) if( it->Boundary() )
	{
		it->UnitNormal(n.data());
		bc[*it][0] = 1;
		bc[*it][1] = 0;
		bc[*it][2] = 1;
		bc[*it][3] = 0;
		if (mean)
		{
			bc[*it][4] = it->Mean(refUmean, 0);
			bc[*it][5] = it->Mean(refVmean, 0);
			bc[*it][6] = it->Mean(refWmean, 0);
		}
		else
		{
			it->Barycenter(x.data());
			bc[*it][4] = refUmean(x.data(),0);
			bc[*it][5] = refVmean(x.data(), 0);
			bc[*it][6] = refWmean(x.data(), 0);
		}
		if (setGP0)
		{
			//it->BackCell().Barycenter(x.data());
			it->Barycenter(x.data());
			variable p = refP(x(0, 0), x(1, 0), x(2, 0));
			gp(0, 0) = p.GetDerivative(0);
			gp(1, 0) = p.GetDerivative(1);
			gp(2, 0) = p.GetDerivative(2);
			gp0[*it] = n.DotProduct(gp);
		}
	}
	
	std::cout << "Saving output to " << fout << std::endl;
	
	m->Save(fout);
	
	return 0;
}
