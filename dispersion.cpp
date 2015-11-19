#include "dispersion.h"
#include "parameters.h"


const double c_light = 299792458*1e2;		// speed of light in vacuum, cm/s

double n_refract_air(double f)
// Cauchy formula, coefficients taken from Akhmanov, Nikitin. Fizicheskaya optika. M. "Nauka", 2004.
{
	double AA = 2.73*1e-4;
	double BB = 7.52*1e-11;
	return 1. + AA*(1. + BB*sqr(f/c_light));
}

/*
double n_refract_air(double f)
{
	double f0 = c_light/lambda;
	double k0 = 2.*M_PI/lambda*1.000276208;
	double k1 = 1./c_light*2.*M_PI*1.000282623;
	double k2 = 2.72659e-31*sqr(2.*M_PI);
	double k3 = 2.8724303e-44;
	
	return c_light/2./M_PI/f*(k0 + k1*(f - f0) + k2/2.*sqr(f - f0) + k3/6.*pow(f - f0, 3));
}
*/
/*
double n_refract_air(double f)
// based on formula in Peck, Reeder. Dispersion of Air. JOSA, v.62, 8, 1972.
// http://www.opticsinfobase.org/view_article.cfm?gotourl=http%3A%2F%2Fwww%2Eopticsinfobase%2Eorg%2FDirectPDFAccess%2F0883FD76%2DCBB5%2D5776%2D85B864A52FA014F8%5F54743%2Epdf%3Fda%3D1%26id%3D54743%26seq%3D0%26mobile%3Dno&org=Moscow%20State%20University%20Scientific%20Lib
{
	double c_l = c_light*1e4;	// c_light must be expressed in um/s
	return 1. + 1e-8*(8060.51 + 2480990/(132.274 - sqr(f/c_l)) + 17455.7/(39.32957 - sqr(f/c_l)));
}
*/

double k_number(double f)
{
	return 2.*M_PI*f/c_light*n_refract_air(f);
}

void CreateDispersionCoeff(std::complex<double> *coeff, int N_t, double ldiff, double ldisp)
// ldiff, ldisp in cm
{
	double freq[N_t];
	double f0 = c_light/lambda;
	double k0 = 2.*M_PI/lambda*1.00027620775000003483512500679;
	double k1 = 1./c_light*2.*M_PI*1.00028262324999994703489392123;
	
//	complex<double> tmp = -_COMPLEX_I*ldiff/ldisp*sqr(M_PI)*2.;
//	for (k = 0; k < N_t/2; k++) freq[k] = k;
//	for (k = N_t/2; k < N_t; k++) freq[k] = k - N_t;
//	for (k = 0; k < N_t; k++) coeff[k] = tmp*sqr(freq[k]);
	
	std::complex<double> tmp = -_COMPLEX_I/2./k0*ldiff;
	for (size_t k = 0; k < N_t/2; ++k) freq[k] = k/grid_size_t*1e15;
	for (size_t k = N_t/2; k < N_t; ++k) freq[k] = (k - N_t)/grid_size_t*1e15;
	for (size_t k = 0; k < N_t; ++k) coeff[k] = tmp*(sqr(k_number(freq[k] + f0)) - sqr(k0 + k1*freq[k]));
}