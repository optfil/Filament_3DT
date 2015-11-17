#ifndef CALCULATION_H_INCLUDED
#define CALCULATION_H_INCLUDED

#include <complex>
#include <cstddef>

void Linear(std::complex<double> *A, size_t n_x, size_t Ntotal, size_t N_t, double z_step);

#endif  // CALCULATION_H_INCLUDED