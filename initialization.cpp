#include "dispersion.h"
#include "parameters.h"
#include "initialization.h"


void CheckParametres()
{/*
	int err_code = -1;
	
	if (size % 2 != 0)
	{
		p_log = fopen(log_name, "a+t");
		fprintf(p_log, "001: Odd number of processes!\n");
		fclose(p_log);
		MPI_Abort(MPI_COMM_WORLD, err_code);
	}
	
	if (NN % size != 0)
	{
		p_log = fopen(log_name, "a+t");
		fprintf(p_log, "002: Sizes of grid and communicator don't match!\n");
		fclose(p_log);
		MPI_Abort(MPI_COMM_WORLD, err_code);
	}*/
}

void InitField(std::complex<double> *Arr, size_t n_x, size_t Ntotal, size_t N_t, double gs, double rad, double gs_t, double rad_t)
{
	const size_t len = N_t*Ntotal;
	size_t idx, idx1;
	double tmp;
	if (N_t == 1)
		for (size_t i = 0; i < n_x; ++i) 
		{
			idx1 = i*len;
			for (size_t j = 0; j < Ntotal; ++j)
			{
				idx = idx1 + j*N_t;
				tmp = exp(-0.5*sqr(gs/rad/Ntotal)*(sqr(i + rank*n_x - (Ntotal - 1.)/2.) + sqr(j - (Ntotal - 1.)/2.)));
				for (size_t k = 0; k < N_t; ++k)
					Arr[idx + k] = tmp;
			}
		}
	else
		for (size_t i = 0; i < n_x; ++i) 
		{
			idx1 = i*len;
			for (size_t j = 0; j < Ntotal; ++j)
			{
				idx = idx1 + j*N_t;
				tmp = exp(-0.5*sqr(gs/rad/Ntotal)*(sqr(i + rank*n_x - (Ntotal - 1.)/2.) + sqr(j - (Ntotal - 1.)/2.)));
				for (size_t k = 0; k < N_t; ++k)
					Arr[idx + k] = tmp*exp(-0.5*sqr(gs_t/rad_t)*sqr(-0.5 + (double)k/(double)(N_t - 1)));
			}
		}
}

void Allocate()
{
	CheckParametres();
	
	EEE = new std::complex<double>[local_nx*NN*NN_t];
	sp_coeff = new std::complex<double>[NN_t];
	fluence = new double[local_nx*NN];
	
	CreateDispersionCoeff(sp_coeff, NN_t, ldiff, ldisp);
	
	plan_f = fftw3d_mpi_create_plan(MPI_COMM_WORLD, NN, NN, NN_t, FFTW_FORWARD, FFTW_ESTIMATE);
	plan_b = fftw3d_mpi_create_plan(MPI_COMM_WORLD, NN, NN, NN_t, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	MPI_Barrier(MPI_COMM_WORLD);
}