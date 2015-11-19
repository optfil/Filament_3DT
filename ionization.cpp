#include "ionization.h"
#include <cmath>
#include <fstream>


IonizationData::IonizationData(const char * filename)
{
	std::fstream f(filename);
	if (!f)
	{
		printf("Can't open ionization data file %s!\n", filename);
		exit(1);
	}

	f >> imin_log >> imax_log >> i_step;
	imin_log -= 4.0;  // convert W/m^2 -> W/cm^2 (4 orders)
	imax_log -= 4.0;  // convert W/m^2 -> W/cm^2 (4 orders)
	
	size_t data_size = (int)((imax_log - imin_log) / i_step) + 1; 
	data.reserve(data_size);
	double tmp;
	for (int i = 0; i < data_size; ++i)
	{
		f >> tmp;
		data.push_back(std::pow(10.0, tmp - 15.0));  // convert s^-1 to fs^-1 (15 orders)
	}
}

double IonizationData::getRate(double intensity) const
{
	int idx = (int)((log10(intensity) - imin_log) / i_step + 0.5);
	if (idx < 0) 
		idx = 0;
	if (idx >= data.size())
		idx = data.size() - 1;

	return data[idx];
}

  