#include <cstdlib>
#include <fstream>

#include "misc.h"
#include "parameters.h"
#include "initialization.h"
#include "propagation.h"
#include "finalization.h"


void GetParam(int argc, char **argv)
{
	if (argc < 2)
		std::exit(1);  // not enough arguments
	dir_name = std::string(argv[1]);

	local_nx = NN / size;
	min_dxdy = grid_size / NN;
	min_dt = grid_size_t / NN_t;

	ldisp = sqr(grid_size_t*1e-15) / k2 * sqr(2.*M_PI);
	ldiff = 2. * M_PI / lambda * sqr(grid_size);
	dz = 0.05*sqr(rad/grid_size);// * ldisp / ldiff * sqr(rad_t / grid_size_t);//
}

void OutParam()
{
	std::string str(dir_name + "/param.txt");
	std::ofstream f(str);
	if (!f)
	{
		std::cout << "Can't open file " << str << "!\n";
		exit(1);
	}
	f << "NN = " << NN << "\n";
	f << "rad = " << rad << " cm\n";
	f << "grid_size = " << grid_size << " cm\n";
	f << "min_dxdy = " << min_dxdy * 10000 <<" um\n\n";
	
	f << "NN_t = " << NN_t << "\n";
	f << "rad_t = " << rad_t << " fs\n";
	f << "grid_size_t = " << grid_size_t << " fs\n";
	f << "min_dt = " << min_dt << " fs\n\n";

	f << "lambda = " << lambda*1e4 << " um\n";
	f << "diffraction length (grid) = " << ldiff << " cm\n";
	f << "diffraction length (beam) = " << ldiff * sqr(rad / grid_size) << " cm\n";
	f << "dispersion length (grid) = " << ldisp << " cm\n";
	f << "dispersion length (pulse) = " << ldisp * sqr(rad_t / grid_size_t) << " cm\n\n";
	
	f << "number of processes  = " << size << "\n";
}

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	GetParam(argc, argv);
	
	if (rank == 0) 
		OutParam();

	MPI_Barrier(MPI_COMM_WORLD);
	
	Allocate();
	Propagate();
	Finalize();

	MPI_Finalize();
}