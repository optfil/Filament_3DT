#include "calculation.h"
#include "parameters.h"
#include <fftw_mpi.h>
/*
void DampingField(complex<double> *A, size_t n_x, size_t Ntotal, double rnb = 8., double degree = 4.)
{
	size_t i_glob;
	double betta;
	for (size_t i_loc = 0; i_loc < n_x; ++i_loc)
	{
		i_glob = i_loc + rank*n_x;
		for (size_t j = 0; j < Ntotal; ++j)
		{
			betta = exp(-pow((sqrt(sqr(i_glob - Ntotal/2) + sqr(j - Ntotal/2)) - Ntotal/2)/rnb, degree));
			A[i_loc*Ntotal + j] *= exp(-betta);
		}
	}
}

void DampingSpectrum(complex<double> *A, size_t n_x, size_t Ntotal, double rnb = 16., double degree = 2.)
{
	size_t i_glob;
	double betta;
	for (size_t i_loc = 0; i_loc < n_x; ++i_loc)
	{
		i_glob = i_loc + rank*n_x;
		for (size_t j = 0; j < Ntotal; ++j)
		{
			betta = exp(-pow((sqr(i_glob - Ntotal/2) + sqr(j - Ntotal/2))/sqr(rnb), degree));
			A[i_loc*Ntotal + j] *= exp(-betta);
		}
	}
}
*/
void Linear(std::complex<double> *A, size_t n_x, size_t Ntotal, size_t N_t, double z_step)
{
	size_t i_glob, idx, idx1;	
	const std::complex<double> tmp1 = 2.*sqr(M_PI)*z_step*_COMPLEX_I;
	std::complex<double> tmp;

	fftw_complex *A_fftw = reinterpret_cast<fftw_complex*>(A);
	fftwnd_mpi(plan_f, 1, A_fftw, NULL, FFTW_TRANSPOSED_ORDER);
	A = reinterpret_cast<std::complex<double>*>(A_fftw);
	
	if (rank < size/2)										// even number of processes only!
		for (size_t i_loc = 0; i_loc < n_x; ++i_loc)
		{
			i_glob = i_loc + rank*n_x;
			idx1 = i_loc*Ntotal*N_t;
			for (size_t j = 0; j < Ntotal/2; ++j) 
			{
				idx = idx1 + j*N_t;
				tmp = exp((double)(sqr(i_glob) + sqr(j))*tmp1);
				for (size_t k = 0; k < N_t; ++k)
					A[idx + k] *= tmp*exp(sp_coeff[k]*z_step);
			}
			for (size_t j = Ntotal/2; j < Ntotal; ++j) 
			{
				idx = idx1 + j*N_t;
				tmp = exp((double)(sqr(i_glob) + sqr(j - Ntotal))*tmp1);
				for (size_t k = 0; k < N_t; ++k)
					A[idx + k] *= tmp*exp(sp_coeff[k]*z_step);
			}
		}
	else
		for (size_t i_loc = 0; i_loc < n_x; ++i_loc)
		{
			i_glob = i_loc + rank*n_x - Ntotal;
			idx1 = i_loc*Ntotal*N_t;
			for (size_t j = 0; j < Ntotal/2; ++j) 
			{
				idx = idx1 + j*N_t;
				tmp = exp((double)(sqr(i_glob) + sqr(j))*tmp1);
				for (size_t k = 0; k < N_t; ++k)
					A[idx + k] *= tmp*exp(sp_coeff[k]*z_step);
			}
			for (size_t j = Ntotal/2; j < Ntotal; ++j) 
			{
				idx = idx1 + j*N_t;
				tmp = exp((double)(sqr(i_glob) + sqr(j - Ntotal))*tmp1);
				for (size_t k = 0; k < N_t; ++k)
					A[idx + k] *= tmp*exp(sp_coeff[k]*z_step);
			}
		}
	
//	DampingSpectrum(A, n_x, Ntotal);
		
	A_fftw = reinterpret_cast<fftw_complex*>(A);
	fftwnd_mpi(plan_b, 1, A_fftw, NULL, FFTW_TRANSPOSED_ORDER);
	A = reinterpret_cast<std::complex<double>*>(A_fftw);

	double buf = 1./sqr(Ntotal)/N_t;
	for (size_t idx = 0; idx < n_x*Ntotal*N_t; ++idx) 
		A[idx] *= buf;
}

