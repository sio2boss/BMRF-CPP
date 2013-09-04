// tcdf and dependent functions taken from:
// 	http://www-stat.stanford.edu/~serban/gxna/src/cdf.cpp

#ifndef CDF_FUNCTIONS_H_
#define CDF_FUNCTIONS_H_

#include "matvec.h"

/**
 * T-test CDF
 */
extern "C" void tcdf(Float1dArray& in, double dfe, Double1dArray& out);


/**
 * Generate CDF
 */
extern "C" void make_cdf(Uint1dArray &freq, Float1dArray &cdf_array, const unsigned int &bootstraps);

/**
 * icdf using numerical recipies
 */
extern "C" void icdf_normal(Double1dArray& in, Double1dArray& out, float mu, float sigma);


/**
 * icdf using numerical recipies
 */
extern "C" void icdf_normal(Double1dArray& in, Double1dArray& out, float mu, float sigma);

extern "C" {
	extern const double cof[];
}

#endif /* CDF_FUNCTIONS_H_ */
