/*
 * rand_functions.h
 *
 *  Created on: Feb 10, 2013
 *      Author: sio2
 */

#ifndef RAND_FUNCTIONS_H_
#define RAND_FUNCTIONS_H_

#include <cmath>

inline float frand()
{
	return (float) rand() / (float) RAND_MAX;
}

inline bool brand(float prob = 0.5)
{
	if (std::isnan(prob) || std::isinf(prob)) return false;
	return (rand() <=  (prob * ((float)RAND_MAX)) );
}

inline unsigned int rrand(unsigned int range)
{
	return rand() % range;
}

#endif /* RAND_FUNCTIONS_H_ */
