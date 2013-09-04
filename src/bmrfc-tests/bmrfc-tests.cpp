//#include "MatlabFile.hpp"
#include <algorithm>
#include <iostream>
#include <string>

#include "cputests.h"
#include "IoFactory.hpp"

using namespace std;



/**
 *
 */
int main(int argc, char** argv) {

	// Need to test all c++ CPU code
	std::srand(std::time(0));

	// Pull stuff from disk
	Float2dArray* ppiArray;
	Float2dArray* geneArray;
	Float2dArray* geneIdArray;
	Float2dArray* geneLabelArray;
	Float2dArray* seedgeneIdArray;

	bool read_success = true;
	read_success = read_success
			&& IoFactory::readArrayFromFile("sim/ppi.mat", "ppiArray",
					ppiArray);
	read_success = read_success
			&& IoFactory::readArrayFromFile("sim/genes.mat", "geneArray",
					geneArray);
	read_success = read_success
			&& IoFactory::readArrayFromFile("sim/genes.mat", "geneIdArray",
					geneIdArray);
	read_success = read_success
			&& IoFactory::readArrayFromFile("sim/genes.mat",
					"geneLabelArray", geneLabelArray);
	read_success = read_success
			&& IoFactory::readArrayFromFile("sim/seed_gene_ids.mat",
					"seedGeneIdArray", seedgeneIdArray);

	if (read_success == false)
		return EXIT_FAILURE;


	// Mutations -------------------------------------------------
	Uint2dArray ppi(ppiArray->y, ppiArray->x);
	transposeFloatToUint((*ppiArray), ppi);

	Uint1dArray geneIdArray1d(geneIdArray->y);
	convert2dTo1dFloatToUint((*geneIdArray), geneIdArray1d);

	Float2dArray geneArrayT_cpu(geneArray->y, geneArray->x);
	transpose(*geneArray, geneArrayT_cpu);

	Double1dArray zscore_genes(geneArrayT_cpu.x);

	// CPU -------------------------------------------------------
	//test matrix/vector stuff
	test_matvec();

	// intersect
	test_intersect();

	// union
	test_union();

	// union
	test_unique();

	// setdiff
	test_setdiff();

	//getppisubnet
	test_getppisubnet();

	// matrix ops
	test_matops();

	//test basic stats
	test_basicstats();

	//test basic stats
	test_cdfs();

	// test gene score
	test_genescore(geneArrayT_cpu, *geneLabelArray, zscore_genes);

	//test connection extraction
	test_g_conn();

	// mrfNetworkScore
	test_scoreNetwork(ppi, geneIdArray1d, zscore_genes);

	// bootstrap
	test_randsample(geneArrayT_cpu, *geneLabelArray);

	// random wrappers
	test_rands();

	//netcand
	test_netcand(ppi, geneIdArray1d, (unsigned int)(*seedgeneIdArray)(0,0));

	// mrfsearchnet
	test_mrfsearchnet(ppi, geneIdArray1d, (unsigned int)(*seedgeneIdArray)(0,0), zscore_genes);


	// Cleanup
	delete ppiArray;
	delete geneArray;
	delete geneIdArray;
	delete geneLabelArray;
	delete seedgeneIdArray;

	return EXIT_SUCCESS;
}

