#include "parameters.h"

char dir_name[1024];

size_t local_nx;
double min_dxdy;
double min_dt;
double dz;
double ldiff, ldisp;

std::complex<double> *EEE;
std::complex<double> *sp_coeff;
double *fluence;

int rank, size;
fftwnd_mpi_plan plan_f, plan_b;
char param_name[1024];