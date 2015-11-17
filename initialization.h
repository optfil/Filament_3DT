#ifndef INITIALIZATION_H_INCLUDED
#define INITIALIZATION_H_INCLUDED

#include <complex>
#include <cstddef>

void InitField(std::complex<double> *Arr, size_t n_x, size_t Ntotal, size_t N_t, double gs, double rad, double gs_t, double rad_t);
void Allocate();

#endif  // INITIALIZATION_H_INCLUDED