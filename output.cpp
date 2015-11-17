#include "output.h"
#include <cstdio>
#include <mpi.h>


void OutRealBinary(const double *source, size_t len, const char *file_name)
{
	FILE *p_f = fopen(file_name, "wb");
	if (!p_f)
		return;
	fwrite(source, sizeof(*source), len, p_f);
	fclose(p_f);
}

void OutRealText(const double *source, size_t len, const char *file_name)
{
	FILE *p_f = fopen(file_name, "wt");
	if (!p_f)
		return;
	for (size_t k = 0; k < len; ++k) 
		fprintf(p_f, "%f\n", source[k]);
	fclose(p_f);
}

void OutRealParallel(const double *source, size_t len, const char *file_name)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_File p_f;
	MPI_Status st;
	MPI_File_open(MPI_COMM_WORLD, file_name, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &p_f);
	MPI_File_set_view(p_f, rank*len*sizeof(*source), MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
	MPI_File_write_all(p_f, source, len, MPI_DOUBLE, &st);
	MPI_File_close(&p_f);
}

void OutComplexBinary(const std::complex<double> *source, size_t len, const char *file_name)
{
	FILE *p_f = fopen(file_name, "wb");
	if (!p_f)
		return;
	fwrite(source, sizeof(*source), 2*len, p_f);
	fclose(p_f);	
}

void OutComplexText(const std::complex<double> *source, size_t len, const char *file_name)
{
	FILE *p_f = fopen(file_name, "wt");
	if (!p_f)
		return;
	for (size_t k = 0; k < len; ++k) 
		fprintf(p_f, "%f\n%f\n", source[k].real(), source[k].imag());
	fclose(p_f);	
}

void OutComplexParallel(const std::complex<double> *source, size_t len, const char *file_name)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_File p_f;
	MPI_Status st;
	MPI_File_open(MPI_COMM_WORLD, file_name, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &p_f);
	MPI_File_set_view(p_f, rank*len*sizeof(*source), MPI_DOUBLE_COMPLEX, MPI_DOUBLE_COMPLEX, "native", MPI_INFO_NULL);
	MPI_File_write_all(p_f, source, len, MPI_DOUBLE_COMPLEX, &st);
	MPI_File_close(&p_f);
}

void InComplexParallel(std::complex<double> *source, size_t len, const char *file_name)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_File p_f;
	MPI_Status st;
	MPI_File_open(MPI_COMM_WORLD, file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &p_f);
	MPI_File_set_view(p_f, rank*len*sizeof(*source), MPI_DOUBLE_COMPLEX, MPI_DOUBLE_COMPLEX, "native", MPI_INFO_NULL);
	MPI_File_read_all(p_f, source, len, MPI_DOUBLE_COMPLEX, &st);
	MPI_File_close(&p_f);
}