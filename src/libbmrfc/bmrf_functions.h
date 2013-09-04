#ifndef BMRF_FUNCTIONS_H_
#define BMRF_FUNCTIONS_H_

#include <cassert>
#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <climits>

#include "matvec.h"
#include "matvec_functions.h"
#include "set_functions.h"
#include "cdf_functions.h"
#include "rand_functions.h"

/*
inline int nextPow2(int x) {
    if (x < 0)
        return 0;
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return x + 1;
}*/

inline int nextpow2(int x)
{
	return 1 << int(ceil(log2(double(x))));
}

/**
 * Take ppiArray and produce a subset based on the list of gids
 */
void genescore(Float2dArray &x, Float2dArray &y, Double1dArray &zscore);

/**
 * ttest, some fixes, and icdf, normalize (will reduce a NxM matrix to 1xM)
 */
void getppisubnet(Uint2dArray &ppi, Uint1dArray &gid, Uint2dArray &sppi);

/**
 * Get scores from genes based on unique sppi, gid should be unique
 */
void extractScores(Uint2dArray &sppi, Uint1dArray &gid, Double1dArray &zscore,
		Uint1dArray &ugeneid, Double1dArray &subscore);

/**
 * Convert ppi array to unique gene ids
 */
void getUniqueGeneIds(Uint2dArray &sppi, Uint1dArray &ugeneid);

/**
 * find connections
 */
void g_conn(Uint2dArray &ppi, Uint1dArray &sgeneid, unsigned int index,
		Uint1dArray &igconn);

/**
 * x = pinv( (coef1*a1) + (coef2*a2) ). a1, a2 and x must be of same size.
 */
void pinvArrayAddWithCoefMul(const double &coef1, Double2dArray &a1,
		const double &coef2, Double2dArray &a2, Double2dArray &x);

void multiplyVectorMatrix(Double1dArray &a, Double2dArray &b, Double1dArray &c);

void multiplyMatrixVector(Double2dArray &a, Double1dArray &b, Double1dArray &c);

void scoreUpdate(Double2dArray &x2d, float gamma, Double1dArray &score,
		Double1dArray& x);

/**
 * (a-b)'*(a-b)
 */
float squareTransposeVectorSubtract(Double1dArray &f, Double1dArray &score);

/**
 * a'*b*a
 */
float squareTransposeMatrixVectorMultiply(Double1dArray &a, Double2dArray &b);

/**
 * Compute the subnetwork score
 * Inputs:
 * gid: gene id vector
 * cand_network_id: gene ids in the interested network
 * zscore: z-score of gene ids in geneid
 * ppi: protein-portein interaction network
 * @return network potentials, negative of network score
 */
double mrfnetscore(Uint1dArray &geneid, Uint1dArray &cand_network_id,
		Double1dArray &zscore, Uint2dArray &ppi);

void populateClasses(Float2dArray &labels, Uint1dArray &index_class1,
		Uint1dArray &index_class2, Float1dArray &classes);

/**
 * Random resampling
 */
void randsample(Float2dArray &data, Uint1dArray &class_labels,
		Float2dArray &randsampled_data);

/**
 * get portion of data based on classes (no random resampling)
 * data [ 367 x 40 ], class_labels [ 20 ] or [ 40 ], sampled_data [ 367 x 20 ];
 * data.x and sampled_data.x must match
 * class_labels.x <= data.y and class_labels.x == sampled_data.y
 */
void classsample(Float2dArray &data, Uint1dArray &class_labels,
		Float2dArray &sampled_data);

/**
 * Generate a candidate network
 */
void netcand(Uint2dArray &ppi, Uint1dArray &geneid,
		Uint1dArray &prev_network_ids, Uint1dArray &prev_network_distances,
		const unsigned int &upper_distance,
		Uint1dArray &candidate_network_id_Array,
		Uint1dArray &candidate_network_distance_Array);

/**
 *
 */
void mrfsearchnet(Uint2dArray &ppi, Uint1dArray &geneids, Double1dArray &zscore,
		const unsigned int &seedgene, const unsigned int &distance,
		float temperature, Uint1dArray &sub_network, double &sub_netscore);

void printProgress(unsigned int count);

void printProgress2(unsigned int count, unsigned int seed);

#endif
