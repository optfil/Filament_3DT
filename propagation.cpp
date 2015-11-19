#include "propagation.h"
#include "initialization.h"
#include "parameters.h"
#include "calculation.h"
#include <sstream>
#include <fstream>
#include <iomanip>


static double cur_z;
static double max_phase, global_max_phase, max_I, global_max_I, max_F, global_max_F;
static long int sch;


void Tact()
{
	++sch;
	cur_z += dz;
	
	for (size_t i = 0; i < local_nx*NN; ++i) 
		fluence[i] = 0.;
	max_F = 0.;
	max_I = 0.;
	max_phase = 0.;

	Linear(EEE, local_nx, NN, NN_t, dz);

	for (size_t i = 0; i < local_nx; ++i)
	{
		size_t idx11 = i*NN, idx21 = i*NN*NN_t;
		for (size_t j = 0; j < NN; ++j)
		{
			size_t idx12 = idx11 + j, idx22 = idx21 + j*NN_t;
			for (size_t k = 0; k < NN_t; ++k)
				fluence[idx12] += norm(EEE[idx22 + k]);
			if (fluence[idx12] > max_F) 
				max_F = fluence[idx12];
		}
	}
}

std::complex<double> direct(double z, double t)
{
	double z_disp = sqr(rad_t*1e-15)/k2*sqr(2.*M_PI);
	double z_diff = sqr(rad)*2.*M_PI/lambda;
	return exp(-sqr(grid_size/2./NN/rad))*1./sqrt(1. + _COMPLEX_I*z*ldiff/z_disp)*1./(1. - _COMPLEX_I*z*ldiff/z_diff)*exp(-0.5*sqr(t*grid_size_t/rad_t)/(1. + _COMPLEX_I*z*ldiff/z_disp));
}

void Save()
{
	if (rank == size/2)
	{
		std::ostringstream ss;
		ss << dir_name << "/Field/intensity" << std::setw(6) << std::setfill('0') << sch << ".bin";
		std::ofstream f(ss.str());
		
		for (size_t k = 0; k < NN_t; ++k)
		{
			std::complex<double> field_e = EEE[NN*NN_t/2 + k], dir = direct(cur_z, (double)k/(NN_t - 1.) - 1./2.);
			f << k << ' ' << cur_z << ' ' 
		      << std::real(field_e) << ' ' << std::imag(field_e) << ' ' << std::norm(field_e) << ' ' << std::arg(field_e) << ' ' 
			  << std::real(dir)     << ' ' << std::imag(dir)     << ' ' << std::norm(dir)     << ' ' << std::arg(dir)     << '\n';
		}
	}

/*	complex<double> *buf;
	buf = new complex<double>[NN*local_nx];
	for (int k = 0; k < NN_t; k++) 
	{
		for (int i = 0; i < local_nx; i++) for (int j = 0; j < NN; j++)
			buf[i*NN + j] = EEE[(i*NN + j)*NN_t + k];
		sprintf(file_name, "%s/log/pole%.3ld%.3d.bin", dir_name, sch, k);
		OutComplexParallel(buf, local_nx*NN, file_name);
	}
	delete[] buf;*/
//	sprintf(file_name, "%s/Fluence/fluence%.6ld.bin", dir_name, sch);
//	OutRealParallel(fluence, local_nx*NN, file_name);
}

void Propagate()
{
	InitField(EEE, local_nx, NN, NN_t, grid_size, rad, grid_size_t, rad_t);
	MPI_Barrier(MPI_COMM_WORLD);
	
	std::string log_name(dir_name + "log.txt");
	
	if (rank == 0)
	{
		std::ofstream log(log_name, std::ios_base::app);
		log << "iteration cur_z z_step phase intensity fluence\n";
	}
	
	sch = 0;
	cur_z = 0.;
	Save();
	while (sch < schmax)
	{
		Tact();
//		MPI_Allreduce(&max_I, &global_max_I, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		MPI_Allreduce(&max_F, &global_max_F, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
//		MPI_Allreduce(&max_phase, &global_max_phase, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		if (sch % save_frequency == 0) 
			Save();
		
		if (rank == 0)
		{
			std::ofstream log(log_name, std::ios_base::app);
			log << sch << " " << cur_z*ldiff << " " << dz*ldiff << " " << global_max_phase << " " << global_max_I << " " << global_max_F << "\n";
		}

//		if ((global_max_phase > PHASE_THRESHOLD_HIGH) || (global_max_phase < PHASE_THRESHOLD_LOW)) dz *= PHASE_THRESHOLD_LOW/global_max_phase;
	}
}