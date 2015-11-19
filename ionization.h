#ifndef IONIZATION_H_INCLUDED
#define IONIZATION_H_INCLUDED

#include <vector>


class IonizationData {
public:
	IonizationData(const char * filename);
	double getRate(double intensity) const;  // intensity must be in W/cm^2
private:
	std::vector<double> data;
	double imin_log;        // log10(I_min)
	double imax_log;        // log10(I_max)
	double i_step;
};

#endif  // IONIZATION_H_INCLUDED