#ifndef PARAMETERS_H_INCLUDED
#define PARAMETERS_H_INCLUDED

#include <cstddef>
#include <string>
#include <mpi.h>
#include <fftw_mpi.h>

#include "misc.h"


extern std::string dir_name;									// series directory

const size_t NN = 1024;      									// whole number of points, square grid
extern size_t local_nx;											// number of rows for a process
const double rad = 0.1;										// beam radius, cm
const double grid_size = 1.;									// size of grid, cm
extern double min_dxdy;										// cross-section step of the grid, ~5 um

const size_t NN_t = 1024;										// number of time lays
const double rad_t = 500;										// pulse half-duration in fs; "1/e-time radius" for intensity
const double grid_size_t = 5000;								// size of grid in time dimension
extern double min_dt;											// time step of the grid, ~1 fs 

const long int schmax = 100;									// anyway, somewhere we must stop...
const int save_frequency = 1;									// number of steps between savings
extern double dz;												// initial step
		
const double lambda = 0.8e-4;									// wavelength, cm
extern double ldiff, ldisp;											// diffraction length, cm
const double k2 = 2.72659401330929268349444781979e-31*sqr(2*M_PI);	// d2k/df2 for f = f0, cannot be calculated from n(f) due to round-off error

extern std::complex<double> *EEE;									// EEE - E-field of pulse part
extern std::complex<double> *sp_coeff;								// coefficients for dispersion, not to recalculate them
extern double *fluence;										// fluence - integral of intensity over time in cross-section

extern int rank, size; 										// rank of process and size of communicator
extern fftwnd_mpi_plan plan_f, plan_b;							// Fourier plans

#endif  // PARAMETERS_H_INCLUDED