#ifndef OUTPUT_H_INCLUDED
#define OUTPUT_H_INCLUDED

#include <complex>
#include <cstddef>


void OutRealBinary(const double *source, size_t len, const char *file_name);
void OutRealText(const double *source, size_t len, const char *file_name);
void OutRealParallel(const double *source, size_t len, const char *file_name);
void OutComplexBinary(const std::complex<double> *source, size_t len, const char *file_name);
void OutComplexText(const std::complex<double> *source, size_t len, const char *file_name);
void OutComplexParallel(const std::complex<double> *source, size_t len, const char *file_name);
void InComplexParallel(std::complex<double> *source, size_t len, const char *file_name);

#endif  // OUTPUT_H_INCLUDED