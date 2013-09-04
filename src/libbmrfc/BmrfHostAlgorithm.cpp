#include "BmrfHostAlgorithm.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cstdio>
#include <cstring>
#include <algorithm>

#include <hdf5.h>

#include <boost/threadpool.hpp>
#include <boost/shared_ptr.hpp>

#include "matvec.h"
#include "matvec_functions.h"
#include "set_functions.h"
#include "bmrf_functions.h"
#include "bmrf_constants.h"
#include "IoFactory.hpp"

using namespace std;

class Job {
public:
	Job(const int bootstrap, const int seed,
			const unsigned int zscore_block_size, Uint2dArray* ppi,
			Uint1dArray* geneIdArray1d, double *zscore_block,
			Float1dArray* seedGeneIdArray1d, const int distance,
			const int temperature, Uint1dArray* bootstrap_network_array,
			Uint2dArray* bootstrap_network_length_array,
			Float2dArray* bootstrap_netscore_array) :
			bootstrap(bootstrap), seed(seed), zscore_block_size(
					zscore_block_size), ppi(ppi), geneIdArray1d(geneIdArray1d), zscore_block(
					zscore_block), seedGeneIdArray1d(seedGeneIdArray1d), distance(
					distance), temperature(temperature), bootstrap_network_array(
					bootstrap_network_array), bootstrap_network_length_array(
					bootstrap_network_length_array), bootstrap_netscore_array(
					bootstrap_netscore_array) {
	}
	~Job() {
	}
	void run() {

		printProgress2(bootstrap, seed);

		// Use zscore block
		double* loc = &zscore_block[bootstrap * zscore_block_size];
		Double1dArray zscore_this(zscore_block_size);
		double* prev = zscore_this.data;
		zscore_this.data = loc;

		Uint1dArray sub_network(MAX_SUBNETWORK_SIZE);
		double sub_netscore = 0.0f;
		mrfsearchnet(*ppi, *geneIdArray1d, zscore_this,
				(*seedGeneIdArray1d)[seed], distance, temperature, sub_network,
				sub_netscore);

		// Copy to big array
		(*bootstrap_netscore_array)(bootstrap, seed) = sub_netscore;
		(*bootstrap_network_length_array)(bootstrap, seed) =
				sub_network.length();
		memcpy(
				&((*bootstrap_network_array).data[(bootstrap
						* (*seedGeneIdArray1d).x * MAX_SUBNETWORK_SIZE)
						+ (seed * MAX_SUBNETWORK_SIZE)]), sub_network.data,
				sub_network.size());
//			 prettyPrintArray("sub_network", sub_network);

// Revert to keep data around
		zscore_this.data = prev;
	}

private:
	const int bootstrap;
	const int seed;
	const unsigned int zscore_block_size;
	Uint2dArray* ppi;
	Uint1dArray* geneIdArray1d;
	double *zscore_block;
	Float1dArray* seedGeneIdArray1d;
	const int distance;
	const int temperature;
	Uint1dArray* bootstrap_network_array;
	Uint2dArray* bootstrap_network_length_array;
	Float2dArray* bootstrap_netscore_array;

};

/**
 *
 */
BmrfHostAlgorithm::BmrfHostAlgorithm(BmrfParams* params) :
		IBmrfAlgorithm(params) {
	timer = new Timer();
}

/**
 *
 */
BmrfHostAlgorithm::~BmrfHostAlgorithm() {
	if (ppiArray != NULL)
		delete ppiArray;
	if (geneArray != NULL)
		delete geneArray;
	if (geneIdArray != NULL)
		delete geneIdArray;
	if (geneLabelArray != NULL)
		delete geneLabelArray;
	if (seedgeneIdArray != NULL)
		delete seedgeneIdArray;
	if (bmrfNetworkScoreArray != NULL)
		delete bmrfNetworkScoreArray;
	if (bmrfNetworkIdArray != NULL)
		delete bmrfNetworkIdArray;
	if (bmrfNetworkLengthArray != NULL)
		delete bmrfNetworkLengthArray;

	delete timer;
}

/*
 *
 */
bool BmrfHostAlgorithm::loadData() {

	// Read HDF5 formatted files
	bool read_success = true;

	printf("\n");
	read_success = read_success
			&& IoFactory::readArrayFromFile(params->ppiFilename, "ppiArray",
					ppiArray);
	read_success = read_success
			&& IoFactory::readArrayFromFile(params->genesFilename, "geneArray",
					geneArray);
	read_success = read_success
			&& IoFactory::readArrayFromFile(params->genesFilename,
					"geneIdArray", geneIdArray);
	read_success = read_success
			&& IoFactory::readArrayFromFile(params->genesFilename,
					"geneLabelArray", geneLabelArray);
	read_success = read_success
			&& IoFactory::readArrayFromFile(params->seedsFilename,
					"seedGeneIdArray", seedgeneIdArray);

	return read_success;
}

/**
 *
 */
bool BmrfHostAlgorithm::run() {

	std::srand(std::time(0));

	if (params->multi == true)
		return runHostMulti();
	else
		return runHost();
}

/**
 *
 */
bool BmrfHostAlgorithm::runHost() {

	printf("\nBMRF on Host (CPU):\n");

	// Procedure: 1) Normalize / compute initial zscore
	printf("  Normalize / Zscore Calculation...\n");
	timer->start();

	// Transpose
	Uint2dArray ppi(ppiArray->y, ppiArray->x);
	transposeFloatToUint((*ppiArray), ppi);

	// normalize and log2
	Float2dArray normgeneArray(geneArray->y, geneArray->x);
	transpose((*geneArray), normgeneArray);
	for (unsigned int i = 0; i < normgeneArray.x; ++i) {
		for (unsigned int j = 0; j < normgeneArray.y; ++j) {
			normgeneArray(i, j) = log2(normgeneArray(i, j) + 4);
		}
	}

	Float1dArray geneMean(normgeneArray.x);
	Float1dArray geneStddev(normgeneArray.x);
	mean(normgeneArray, geneMean);
	stddev(normgeneArray, geneMean, geneStddev);
	for (int i = 0; i < normgeneArray.x; ++i) {
		for (int j = 0; j < normgeneArray.y; ++j) {
			normgeneArray(i, j) = (normgeneArray(i, j) - geneMean[i])
					/ geneStddev[i];
		}
	}

	// Divide by class
	Float1dArray classes(2);
	classes[0] = 1.0f;
	classes[1] = 2.0f;
	Uint1dArray index_class1((*geneLabelArray).x);
	Uint1dArray index_class2((*geneLabelArray).x);
	populateClasses((*geneLabelArray), index_class1, index_class2, classes);

	Uint1dArray index_classall(index_class1.length() + index_class2.length());
	unionSets(index_class1, index_class2, index_classall);

	Float2dArray normGeneArrayClass1_cpu(normgeneArray.x,
			index_class1.length());
	classsample(normgeneArray, index_class1, normGeneArrayClass1_cpu);

	Float2dArray normGeneArrayClass2_cpu(normgeneArray.x,
			index_class2.length());
	classsample(normgeneArray, index_class2, normGeneArrayClass2_cpu);

	// Zscore
	Double1dArray zscore(normgeneArray.x);
	genescore(normGeneArrayClass1_cpu, normGeneArrayClass2_cpu, zscore);
	Double1dArray zscore0(zscore.max);
	copyArray(zscore, zscore0);

	// Convert gene ids
	Uint1dArray geneIdArray1d(geneIdArray->y);
	convert2dTo1dFloatToUint((*geneIdArray), geneIdArray1d);

	// netscores
	Float1dArray seedGeneIdArray1d(seedgeneIdArray->y);
	convert2dTo1d((*seedgeneIdArray), seedGeneIdArray1d);
	timer->mark("normalize");

	// Procedure: 2) Bootstrap
	Uint1dArray bootstrap_network_array(
			params->bootstraps * seedGeneIdArray1d.x * MAX_SUBNETWORK_SIZE);
	Uint2dArray bootstrap_network_length_array(params->bootstraps,
			seedGeneIdArray1d.x);
	Float2dArray bootstrap_netscore_array(params->bootstraps,
			seedGeneIdArray1d.x);

	printf("  Bootstrap: (%i, %i): ", params->bootstraps, seedGeneIdArray1d.x);
	for (int bootstrap = 0; bootstrap < params->bootstraps; ++bootstrap) {
		printProgress(bootstrap);

		// Compute Zscore with random sampling
		Float2dArray rand_class1(normgeneArray.x, index_class1.length());
		Float2dArray rand_class2(normgeneArray.x, index_class2.length());
		randsample(normgeneArray, index_class1, rand_class1);
		randsample(normgeneArray, index_class2, rand_class2);
		genescore(rand_class1, rand_class2, zscore);
		for (int seed = 0; seed < seedGeneIdArray1d.x; ++seed) {

			Uint1dArray sub_network(MAX_SUBNETWORK_SIZE);
			double sub_netscore = 0.0f;
			mrfsearchnet(ppi, geneIdArray1d, zscore, seedGeneIdArray1d[seed],
					params->distance, params->temperature, sub_network,
					sub_netscore);

			// Copy to big array
			bootstrap_netscore_array(bootstrap, seed) = sub_netscore;
			bootstrap_network_length_array(bootstrap, seed) =
					sub_network.length();
			memcpy(
					&(bootstrap_network_array.data[(bootstrap
							* seedGeneIdArray1d.x * MAX_SUBNETWORK_SIZE)
							+ (seed * MAX_SUBNETWORK_SIZE)]), sub_network.data,
					sub_network.size());
//			prettyPrintArray("sub_network", sub_network);
		}
	}
	printf("\n");
	timer->mark("bootstrap");

	// Procedure: 3) Control group bootstrap
	Uint1dArray baseline_network_array(
			params->bootstraps * seedGeneIdArray1d.x * MAX_SUBNETWORK_SIZE);
	Uint2dArray baseline_network_length_array(params->bootstraps,
			seedGeneIdArray1d.x);
	Float2dArray baseline_netscore_array(params->bootstraps,
			seedGeneIdArray1d.x);
	printf("  Baseline: (%i, %i): ", params->bootstraps, seedGeneIdArray1d.x);
	fflush(stdout);
	for (int bootstrap = 0; bootstrap < params->bootstraps; ++bootstrap) {
		printProgress(bootstrap);

		// Compute Zscore with random sampling without respect to class
		Float2dArray rand_class1(normgeneArray.x, index_class1.length());
		Float2dArray rand_class2(normgeneArray.x, index_class2.length());
		randsample(normgeneArray, index_classall, rand_class1);
		randsample(normgeneArray, index_classall, rand_class2);
		genescore(rand_class1, rand_class2, zscore);

		for (int seed = 0; seed < seedGeneIdArray1d.x; ++seed) {

			Uint1dArray sub_network(MAX_SUBNETWORK_SIZE);
			double sub_netscore = 0.0f;
			mrfsearchnet(ppi, geneIdArray1d, zscore, seedGeneIdArray1d[seed],
					params->distance, params->temperature, sub_network,
					sub_netscore);

			// Copy to big array
			baseline_netscore_array(bootstrap, seed) = sub_netscore;
			baseline_network_length_array(bootstrap, seed) =
					sub_network.length();
			memcpy(
					&(baseline_network_array.data[(bootstrap
							* seedGeneIdArray1d.x * MAX_SUBNETWORK_SIZE)
							+ (seed * MAX_SUBNETWORK_SIZE)]), sub_network.data,
					sub_network.size());
//			prettyPrintArray("sub_network", sub_network);

		}
	}
	printf("\n");
	timer->mark("baseline");

	// Procedure: 4) Credibility analysis
	bmrfNetworkScoreArray = new Float2dArray(seedGeneIdArray1d.x, 1);
	bmrfNetworkLengthArray = new Float2dArray(seedGeneIdArray1d.x, 1);
	bmrfNetworkIdArray = new Float2dArray(seedGeneIdArray1d.x, ppi.x);
	printf("  Credibility Analysis...\n");
	for (int seed = 0; seed < seedGeneIdArray1d.x; ++seed) {
		Uint1dArray temp_network(ppi.x);
		temp_network.shrink(0);
		Uint1dArray temp_gid(ppi.x);
		temp_gid.shrink(0);

		// Gid
		Uint1dArray dect_gid_all(ppi.x);
		dect_gid_all.shrink(0);
		Uint1dArray dect_gid_all_baseline(ppi.x);
		dect_gid_all_baseline.shrink(0);
		unsigned int length = 0;
		for (int bootstrap = 0; bootstrap < params->bootstraps; ++bootstrap) {

			length = bootstrap_network_length_array(bootstrap, seed);
			memcpy(temp_network.data,
					&(bootstrap_network_array.data[(bootstrap
							* seedGeneIdArray1d.x * MAX_SUBNETWORK_SIZE)
							+ (seed * MAX_SUBNETWORK_SIZE)]),
					length * sizeof(unsigned int));
			temp_network.shrink(length);
			unionSets(dect_gid_all, temp_network, temp_gid);
			copyArray(temp_gid, dect_gid_all);

			length = baseline_network_length_array(bootstrap, seed);
			memcpy(temp_network.data,
					&(baseline_network_array.data[(bootstrap
							* seedGeneIdArray1d.x * MAX_SUBNETWORK_SIZE)
							+ (seed * MAX_SUBNETWORK_SIZE)]),
					length * sizeof(unsigned int));
			unionSets(dect_gid_all_baseline, temp_network, temp_gid);
			copyArray(temp_gid, dect_gid_all_baseline);
		}

		// Frequency
		Uint1dArray freq_dect(dect_gid_all.x);
		memset(freq_dect.data, 0, freq_dect.size());
		Uint1dArray freq_dect_baseline(dect_gid_all_baseline.x);
		memset(freq_dect_baseline.data, 0, freq_dect_baseline.size());
		Uint1dArray b(dect_gid_all.x);
		Uint1dArray b_baseline(dect_gid_all_baseline.x);
		for (int bootstrap = 0; bootstrap < params->bootstraps; ++bootstrap) {
			length = bootstrap_network_length_array(bootstrap, seed);
			memcpy(temp_network.data,
					&(bootstrap_network_array.data[(bootstrap
							* seedGeneIdArray1d.x * MAX_SUBNETWORK_SIZE)
							+ (seed * MAX_SUBNETWORK_SIZE)]),
					length * sizeof(unsigned int));
			temp_network.shrink(length);
			intersectSetsReturnIndex(dect_gid_all, temp_network, b);
			for (int i = 0; i < b.x; ++i)
				freq_dect[b[i]] = freq_dect[b[i]] + 1;

			length = baseline_network_length_array(bootstrap, seed);
			memcpy(temp_network.data,
					&(baseline_network_array.data[(bootstrap
							* seedGeneIdArray1d.x * MAX_SUBNETWORK_SIZE)
							+ (seed * MAX_SUBNETWORK_SIZE)]),
					length * sizeof(unsigned int));
			temp_network.shrink(length);
			intersectSetsReturnIndex(dect_gid_all_baseline, temp_network,
					b_baseline);
			for (int i = 0; i < b_baseline.x; ++i)
				freq_dect_baseline[b_baseline[i]] =
						freq_dect_baseline[b_baseline[i]] + 1;
		}

		// calculate the CDF of condidence, discard the seed gene
		Float1dArray cdf_obs(100);
		Float1dArray cdf_baseline(100);
		make_cdf(freq_dect, cdf_obs, params->bootstraps);
		make_cdf(freq_dect_baseline, cdf_baseline, params->bootstraps);

		// calculate FDR
		Float1dArray fdr(cdf_baseline.x);
		scalarDivideArray(cdf_baseline, cdf_obs, fdr);
		unsigned int thrd_freq = 15;
		for (int i = 0; i < fdr.x; ++i) {
			if (fdr[i] <= 0.05f) {
				thrd_freq = max(i, 15);
				break;
			}
		}

		Uint1dArray BMRF_network_ID(MAX_SUBNETWORK_SIZE);
		BMRF_network_ID.shrink(0);
		for (int i = 0; i < dect_gid_all.x; ++i) {
			if (freq_dect[i] >= thrd_freq)
				BMRF_network_ID.append(dect_gid_all[i]);
		}
		float BMRF_network_score = -1.0f
				* mrfnetscore(geneIdArray1d, BMRF_network_ID, zscore0, ppi);

// Roll-up all good subnetworks
		printf("BMRF_network_ID @ %f = ", BMRF_network_score);
		BMRF_network_ID.print(BMRF_network_ID.x);

		(*bmrfNetworkScoreArray)(seed, 0) = BMRF_network_score;
		(*bmrfNetworkLengthArray)(seed, 0) = BMRF_network_ID.x;
		for (int i = 0; i < BMRF_network_ID.x; ++i)
			(*bmrfNetworkIdArray)(seed, i) = float(BMRF_network_ID[i]);

	}

	timer->mark("credibility");
	timer->stop();

	// Display performance
	timer->print();

	return true;
}

/**
 *
 */
bool BmrfHostAlgorithm::runHostMulti() {

	printf("\nBMRF on Host (CPU) Multithreaded:\n");

	// Procedure: 1) Normalize / compute initial zscore
	printf("  Normalize / Zscore Calculation...\n");
	timer->start();

	// Transpose
	Uint2dArray ppi(ppiArray->y, ppiArray->x);
	transposeFloatToUint((*ppiArray), ppi);

	// normalize and log2
	Float2dArray normgeneArray(geneArray->y, geneArray->x);
	transpose((*geneArray), normgeneArray);
	for (unsigned int i = 0; i < normgeneArray.x; ++i) {
		for (unsigned int j = 0; j < normgeneArray.y; ++j) {
			normgeneArray(i, j) = log2(normgeneArray(i, j) + 4);
		}
	}

	Float1dArray geneMean(normgeneArray.x);
	Float1dArray geneStddev(normgeneArray.x);
	mean(normgeneArray, geneMean);
	stddev(normgeneArray, geneMean, geneStddev);
	for (int i = 0; i < normgeneArray.x; ++i) {
		for (int j = 0; j < normgeneArray.y; ++j) {
			normgeneArray(i, j) = (normgeneArray(i, j) - geneMean[i])
					/ geneStddev[i];
		}
	}

	// Divide by class
	Float1dArray classes(2);
	classes[0] = 1.0f;
	classes[1] = 2.0f;
	Uint1dArray index_class1((*geneLabelArray).x);
	Uint1dArray index_class2((*geneLabelArray).x);
	populateClasses((*geneLabelArray), index_class1, index_class2, classes);

	Uint1dArray index_classall(index_class1.length() + index_class2.length());
	unionSets(index_class1, index_class2, index_classall);

	Float2dArray normGeneArrayClass1_cpu(normgeneArray.x,
			index_class1.length());
	classsample(normgeneArray, index_class1, normGeneArrayClass1_cpu);

	Float2dArray normGeneArrayClass2_cpu(normgeneArray.x,
			index_class2.length());
	classsample(normgeneArray, index_class2, normGeneArrayClass2_cpu);

	// Zscore
	Double1dArray zscore(normgeneArray.x);
	genescore(normGeneArrayClass1_cpu, normGeneArrayClass2_cpu, zscore);
	Double1dArray zscore0(zscore.max);
	copyArray(zscore, zscore0);

	// Convert gene ids
	Uint1dArray geneIdArray1d(geneIdArray->y);
	convert2dTo1dFloatToUint((*geneIdArray), geneIdArray1d);

	// Seeds
	Float1dArray seedGeneIdArray1d(seedgeneIdArray->y);
	convert2dTo1d((*seedgeneIdArray), seedGeneIdArray1d);
	timer->mark("normalize");

	// Build zscores for bootstrams and baselines
	const unsigned int zscore_block_size = nextpow2(normgeneArray.x);
	printf("  Generate all Zscores, zscore_block_size = %i\n",
			zscore_block_size);
	double* zscore_block;
	allocateCPUMemory(zscore_block,
			params->bootstraps * 2 * zscore_block_size * sizeof(double));

	for (int bootstrap = 0; bootstrap < params->bootstraps; ++bootstrap) {

		// Index into block
		double* loc = &zscore_block[bootstrap * zscore_block_size];
		Double1dArray zscore_this(normgeneArray.x);
		double* prev = zscore_this.data;
		zscore_this.data = loc;

		// Compute Zscore with random sampling
		Float2dArray rand_class1(normgeneArray.x, index_class1.length());
		Float2dArray rand_class2(normgeneArray.x, index_class2.length());
		randsample(normgeneArray, index_class1, rand_class1);
		randsample(normgeneArray, index_class2, rand_class2);
		genescore(rand_class1, rand_class2, zscore_this);

		// Revert to keep data around
		zscore_this.data = prev;

	}
	for (int bootstrap = 0; bootstrap < params->bootstraps; ++bootstrap) {

		// Index into block
		double* loc = &zscore_block[params->bootstraps * zscore_block_size
				+ bootstrap * zscore_block_size];
		Double1dArray zscore_this(normgeneArray.x);
		double* prev = zscore_this.data;
		zscore_this.data = loc;

		// Compute Zscore with random sampling without respect to class
		Float2dArray rand_class1(normgeneArray.x, index_class1.length());
		Float2dArray rand_class2(normgeneArray.x, index_class2.length());
		randsample(normgeneArray, index_classall, rand_class1);
		randsample(normgeneArray, index_classall, rand_class2);
		genescore(rand_class1, rand_class2, zscore_this);

		// Revert to keep data around
		zscore_this.data = prev;

	}
	timer->mark("zscores");

	// Procedure: 2) Bootstrap
	Uint1dArray bootstrap_network_array(
			params->bootstraps * seedGeneIdArray1d.x * MAX_SUBNETWORK_SIZE);
	Uint2dArray bootstrap_network_length_array(params->bootstraps,
			seedGeneIdArray1d.x);
	Float2dArray bootstrap_netscore_array(params->bootstraps,
			seedGeneIdArray1d.x);

	printf("  Bootstrap: (%i, %i): ", params->bootstraps, seedGeneIdArray1d.x);
	fflush(stdout);

	// Add some tasks to the pool.
	boost::threadpool::pool bootstrap_tp(8);
	for (int bootstrap = 0; bootstrap < params->bootstraps; ++bootstrap) {
		for (int seed = 0; seed < seedGeneIdArray1d.x; ++seed) {
			boost::shared_ptr<Job> job(
					new Job(bootstrap, seed, zscore_block_size, &ppi,
							&geneIdArray1d, zscore_block, &seedGeneIdArray1d,
							params->distance, params->temperature,
							&bootstrap_network_array,
							&bootstrap_network_length_array,
							&bootstrap_netscore_array));
			bootstrap_tp.schedule(boost::bind(&Job::run, job));
		}

	}
	bootstrap_tp.wait();
	printf("\n");
	timer->mark("bootstrap");

	// Procedure: 3) Control group bootstrap
	Uint1dArray baseline_network_array(
			params->bootstraps * seedGeneIdArray1d.x * MAX_SUBNETWORK_SIZE);
	Uint2dArray baseline_network_length_array(params->bootstraps,
			seedGeneIdArray1d.x);
	Float2dArray baseline_netscore_array(params->bootstraps,
			seedGeneIdArray1d.x);

	printf("  Baseline: (%i, %i): ", params->bootstraps, seedGeneIdArray1d.x);
	fflush(stdout);

	boost::threadpool::pool baseline_tp(8);
	for (int baseline = 0; baseline < params->bootstraps; ++baseline) {
		for (int seed = 0; seed < seedGeneIdArray1d.x; ++seed) {
			boost::shared_ptr<Job> job(
					new Job(baseline, seed, zscore_block_size, &ppi,
							&geneIdArray1d,
							zscore_block
									+ (params->bootstraps * zscore_block_size),
							&seedGeneIdArray1d, params->distance,
							params->temperature, &baseline_network_array,
							&baseline_network_length_array,
							&baseline_netscore_array));
			baseline_tp.schedule(boost::bind(&Job::run, job));
		}
	}
	baseline_tp.wait();
	printf("\n");
	timer->mark("baseline");

	// Procedure: 4) Credibility analysis
	bmrfNetworkScoreArray = new Float2dArray(seedGeneIdArray1d.x, 1);
	bmrfNetworkLengthArray = new Float2dArray(seedGeneIdArray1d.x, 1);
	bmrfNetworkIdArray = new Float2dArray(seedGeneIdArray1d.x, ppi.x);
	printf("  Credibility Analysis...\n");
	for (int seed = 0; seed < seedGeneIdArray1d.x; ++seed) {
		Uint1dArray temp_network(ppi.x);
		temp_network.shrink(0);
		Uint1dArray temp_gid(ppi.x);
		temp_gid.shrink(0);

		// Gid
		Uint1dArray dect_gid_all(ppi.x);
		dect_gid_all.shrink(0);
		Uint1dArray dect_gid_all_baseline(ppi.x);
		dect_gid_all_baseline.shrink(0);
		unsigned int length = 0;
		for (int bootstrap = 0; bootstrap < params->bootstraps; ++bootstrap) {

			length = bootstrap_network_length_array(bootstrap, seed);
			memcpy(temp_network.data,
					&(bootstrap_network_array.data[(bootstrap
							* seedGeneIdArray1d.x * MAX_SUBNETWORK_SIZE)
							+ (seed * MAX_SUBNETWORK_SIZE)]),
					length * sizeof(unsigned int));
			temp_network.shrink(length);
			unionSets(dect_gid_all, temp_network, temp_gid);
			copyArray(temp_gid, dect_gid_all);

			length = baseline_network_length_array(bootstrap, seed);
			memcpy(temp_network.data,
					&(baseline_network_array.data[(bootstrap
							* seedGeneIdArray1d.x * MAX_SUBNETWORK_SIZE)
							+ (seed * MAX_SUBNETWORK_SIZE)]),
					length * sizeof(unsigned int));
			unionSets(dect_gid_all_baseline, temp_network, temp_gid);
			copyArray(temp_gid, dect_gid_all_baseline);
		}

		// Frequency
		Uint1dArray freq_dect(dect_gid_all.x);
		memset(freq_dect.data, 0, freq_dect.size());
		Uint1dArray freq_dect_baseline(dect_gid_all_baseline.x);
		memset(freq_dect_baseline.data, 0, freq_dect_baseline.size());
		Uint1dArray b(dect_gid_all.x);
		Uint1dArray b_baseline(dect_gid_all_baseline.x);
		for (int bootstrap = 0; bootstrap < params->bootstraps; ++bootstrap) {
			length = bootstrap_network_length_array(bootstrap, seed);
			memcpy(temp_network.data,
					&(bootstrap_network_array.data[(bootstrap
							* seedGeneIdArray1d.x * MAX_SUBNETWORK_SIZE)
							+ (seed * MAX_SUBNETWORK_SIZE)]),
					length * sizeof(unsigned int));
			temp_network.shrink(length);
			intersectSetsReturnIndex(dect_gid_all, temp_network, b);
			for (int i = 0; i < b.x; ++i)
				freq_dect[b[i]] = freq_dect[b[i]] + 1;

			length = baseline_network_length_array(bootstrap, seed);
			memcpy(temp_network.data,
					&(baseline_network_array.data[(bootstrap
							* seedGeneIdArray1d.x * MAX_SUBNETWORK_SIZE)
							+ (seed * MAX_SUBNETWORK_SIZE)]),
					length * sizeof(unsigned int));
			temp_network.shrink(length);
			intersectSetsReturnIndex(dect_gid_all_baseline, temp_network,
					b_baseline);
			for (int i = 0; i < b_baseline.x; ++i)
				freq_dect_baseline[b_baseline[i]] =
						freq_dect_baseline[b_baseline[i]] + 1;
		}

		// calculate the CDF of condidence, discard the seed gene
		Float1dArray cdf_obs(100);
		Float1dArray cdf_baseline(100);
		make_cdf(freq_dect, cdf_obs, params->bootstraps);
		make_cdf(freq_dect_baseline, cdf_baseline, params->bootstraps);

		// calculate FDR
		Float1dArray fdr(cdf_baseline.x);
		scalarDivideArray(cdf_baseline, cdf_obs, fdr);
		unsigned int thrd_freq = 15;
		for (int i = 0; i < fdr.x; ++i) {
			if (fdr[i] <= 0.05f) {
				thrd_freq = max(i, 15);
				break;
			}
		}

		Uint1dArray BMRF_network_ID(ppi.x);
		BMRF_network_ID.shrink(0);
		for (int i = 0; i < dect_gid_all.x; ++i) {
			if (freq_dect[i] >= thrd_freq)
				BMRF_network_ID.append(dect_gid_all[i]);
		}
		float BMRF_network_score = -1.0f
				* mrfnetscore(geneIdArray1d, BMRF_network_ID, zscore0, ppi);

//        printf("cdf_obs = "); cdf_obs.print(cdf_obs.x);
//        printf("\n");
//        printf("cdf_baseline = "); cdf_baseline.print(cdf_baseline.x);
//        printf("\n");
//        printf("fdr = "); fdr.print(fdr.x);

// Roll-up all good subnetworks
		printf("BMRF_network_ID @ %f = ", BMRF_network_score);
		BMRF_network_ID.print(BMRF_network_ID.x);

		(*bmrfNetworkScoreArray)(seed, 0) = BMRF_network_score;
		(*bmrfNetworkLengthArray)(seed, 0) = BMRF_network_ID.x;
		for (int i = 0; i < BMRF_network_ID.x; ++i)
			(*bmrfNetworkIdArray)(seed, i) = float(BMRF_network_ID[i]);

	}

	timer->mark("credibility");
	timer->stop();

	// Display performance
	timer->print();

	return true;
}

/**
 *
 */
bool BmrfHostAlgorithm::writeResults() {
	// Read HDF5 formatted files if we have data
	if (bmrfNetworkIdArray == NULL)
		return true;
	bool write_success = true;

	printf("\n");

	write_success = write_success
			&& IoFactory::writeArraysToFile(params->outputFilename,
					"bmrfNetworkIdArray", bmrfNetworkIdArray,
					"bmrfNetworkLengthArray", bmrfNetworkLengthArray,
					"bmrfNetworkScoreArray", bmrfNetworkScoreArray);

	return write_success;
}

