/*
 * matvec_functions.h
 *
 *  Created on: Feb 8, 2013
 *      Author: sio2
 */

#ifndef MATVEC_FUNCTIONS_H_
#define MATVEC_FUNCTIONS_H_

#include "matvec.h"

/**
 * pretty printing
 */
void prettyPrintArray(const std::string &var, Uint1dArray& array);

/**
 * pretty printing
 */
void prettyPrintArray(const std::string &var, Float1dArray& array);

/**
 * Test set if empty
 */
bool isempty(Uint1dArray &a);

/**
 * copy (replace
 */
int copyArray(Uint1dArray &from, Uint1dArray &to);

/**
 * copy (replace
 */
int copyArray(Float1dArray &from, Float1dArray &to);

/**
 * copy (replace
 */
int copyArray(Float2dArray &from, Float2dArray &to);

/**
 * copy (replace
 */
int copyArray(Double1dArray &from, Double1dArray &to);

/**
 * Extract row
 */
void extractRow(Uint2dArray &from, const unsigned int &row, const unsigned int &length, Uint1dArray &to);

/**
 * 2d to 1d
 */
void convert2dTo1d(Uint2dArray &a, Uint1dArray &b);

/**
 * 2d to 1d
 */
void convert2dTo1d(Float2dArray &a, Float1dArray &b);

void convert2dTo1dFloatToUint(Float2dArray &a, Uint1dArray &b);

/**
 * Transpose matrix: m -> m'
 */
void transpose(Float2dArray &m, Float2dArray& m_trans);


/**
 * Transpose matrix: m -> m'
 */
void transpose(Double2dArray &m, Double2dArray& m_trans);

/**
 * Transpose matrix: m -> m'
 */
void transposeFloatToUint(Float2dArray &m, Uint2dArray& m_trans);

/**
 * Standard matrix multiply
 */
void multiplyMatrix(Float2dArray &a, Float2dArray &b, Float2dArray &c);

/**
 * Standard matrix multiply
 */
void multiplyMatrix(Double2dArray &a, Double2dArray &b, Double2dArray &c);


/**
 * Standard scalar divide (c = a./b)
 */
int scalarDivideArray(Float1dArray &a, Float1dArray &b, Float1dArray &c);


/**
 * 1D mean
 */
unsigned int max(Uint1dArray& p);

/**
 * 1D mean
 */
float mean(Double1dArray& p);

/**
 * 1D mean
 */
float mean(Float1dArray& p);

/**
 * 2d mean (row-wise)
 */
void mean(Float2dArray& in, Float1dArray &out);

/**
 * Standard Deviation 1d
 */
float stddev(Float1dArray &p, const float arraymean);


/**
 * Standard Deviation 1d
 */
float stddev(Double1dArray &p, const double arraymean);

/**
 * Standard Deviation 2d  (row-wise)
 */
void stddev(Float2dArray &in, Float1dArray &arraymean, Float1dArray &out);

/**
 * Variance 1d
 */
float var(Float1dArray &p, const float arraymean);

/**
 * Variance 2d  (row-wise)
 */
void var(Float2dArray& in, Float1dArray &arraymean, Float1dArray &out);

/**
 * Array subtraction
 */
void subtract(Float1dArray &x, Float1dArray &y, Float1dArray &difference);


/**
 * Array subtraction
 */
void subtract(Double1dArray &x, Double1dArray &y, Double1dArray &difference);


/**
 * Use SVD for pseudo inverse of matrix a
 */
void pinv(Double2dArray &a, Double2dArray &a_pinv);

#endif /* MATVEC_FUNCTIONS_H_ */
