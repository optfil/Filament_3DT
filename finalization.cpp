#include "parameters.h"


void Finalize()
{
	delete[] EEE;
	delete[] sp_coeff;
	delete[] fluence;

	fftwnd_mpi_destroy_plan(plan_f);
	fftwnd_mpi_destroy_plan(plan_b);
}