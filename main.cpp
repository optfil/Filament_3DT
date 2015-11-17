#include <mpi.h>
#include <fftw_mpi.h>
#include <complex>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "misc.h"
#include "parameters.h"
#include "initialization.h"
#include "propagation.h"
#include "finalization.h"

//==========//============//==========//==========//==========//

void GetParam(int argc, char **argv)
{
	if (argc < 2)
		std::exit(1);  // not enough arguments
	strcpy(dir_name, argv[1]);

	local_nx = NN / size;
	min_dxdy = grid_size / NN;
	min_dt = grid_size_t / NN_t;

	ldisp = sqr(grid_size_t*1e-15) / k2 * sqr(2.*M_PI);
	ldiff = 2. * M_PI / lambda * sqr(grid_size);
	dz = 0.05*sqr(rad/grid_size);// * ldisp / ldiff * sqr(rad_t / grid_size_t);//
}

void OutParam()
{
	char buf[1024];
	sprintf(buf, "%s/param.txt", dir_name);
	FILE *p_f = fopen(buf, "wt+");
	
	fprintf(p_f, "NN = %d\n", NN);
	fprintf(p_f, "rad = %f cm\n", rad);
	fprintf(p_f, "grid_size = %f cm\n", grid_size);
	fprintf(p_f, "min_dxdy = %f um\n\n", min_dxdy * 10000);
	
	fprintf(p_f, "NN_t = %d\n", NN_t);
	fprintf(p_f, "rad_t = %f fs\n", rad_t);
	fprintf(p_f, "grid_size_t = %f fs\n", grid_size_t);
	fprintf(p_f, "min_dt = %f fs\n\n", min_dt);

	fprintf(p_f, "lambda = %f um\n", lambda*1e4);
	fprintf(p_f, "diffraction length (grid) = %f cm\n", ldiff);
	fprintf(p_f, "diffraction length (beam) = %f cm\n", ldiff * sqr(rad / grid_size));
	fprintf(p_f, "dispersion length (grid) = %f cm\n", ldisp);
	fprintf(p_f, "dispersion length (pulse) = %f cm\n\n", ldisp * sqr(rad_t / grid_size_t));
	
	fprintf(p_f, "number of processes = %d\n", size);
	
	fclose(p_f);
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