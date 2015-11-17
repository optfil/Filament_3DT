#ifndef DISPERSION_H_INCLUDED
#define DISPERSION_H_INCLUDED

#include <complex>


void CreateDispersionCoeff(std::complex<double> *coeff, int N_t, double ldiff, double ldisp);

#endif  // DISPERSION_H_INCLUDED