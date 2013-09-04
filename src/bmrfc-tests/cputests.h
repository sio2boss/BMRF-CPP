/*
 * tests.h
 *
 *  Created on: Feb 8, 2013
 *      Author: sio2
 */

#ifndef CPU_TESTS_H_
#define CPU_TESTS_H_

#include "bmrf_functions.h"
#include "rand_functions.h"
#include "MatlabFile.hpp"

/**
 *
 */
void test_matvec();

/**
 * Test union, must handle no duplicates and duplicates
 */
void test_union();

/**
 * Test intersect, three cases: full, partial and no overlap
 */
void test_intersect();


void test_unique();

/**
 * Test union, must handle no duplicates and duplicates
 */
void test_setdiff();

void test_getppisubnet();

void test_matops();

void test_basicstats();

void test_cdfs();

void test_genescore(Float2dArray &normgenes, Float2dArray &geneLabelArray,
		Double1dArray &zscore_genes);
void test_g_conn();

void test_rands();

void test_scoreNetwork(Uint2dArray &ppi, Uint1dArray &geneid,
		Double1dArray & zscore);

void test_randsample(Float2dArray &genes, Float2dArray &labels);

void test_netcand(Uint2dArray &ppi, Uint1dArray &geneids,
		const unsigned int &seedgene);

void test_mrfsearchnet(Uint2dArray &ppi, Uint1dArray &geneids,
		const unsigned int &seedgene, Double1dArray &zscore_genes);

#endif /* TESTS_H_ */
