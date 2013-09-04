#include "./docopt.c"
#include "errno.h"
#include <string>
#include <iostream>

#include "BmrfParams.hpp"
#include "IBmrfAlgorithm.h"
#include "BmrfHostAlgorithm.hpp"

using namespace std;

// Default parameters
static int DEFAULT_DISTANCE = 2;
static int DEFAULT_TEMPERATURE = 1;
static int DEFAULT_BOOTSTRAPS = 100;

/**
 * Handle commandline arguments, read datafiles, run gpu code
 */
int main(int argc, char *argv[]) {

	DocoptArgs args = docopt(argc, argv, 1, "Use the force Luke\n");

	if (args.help == true || args.version == true) {
		exit(EXIT_SUCCESS);
	}

	if (args.ppi == NULL || args.genes == NULL || args.output == NULL) {
		printf("%s\n", args.help_message);
		exit(EXIT_SUCCESS);
	}

	// Populate params with defaults
	BmrfParams* params = new BmrfParams();
	params->ppiFilename = args.ppi;
	params->genesFilename = args.genes;
	params->seedsFilename = args.seeds;
	params->outputFilename = args.output;
	if (args.distance != NULL)
		params->distance = atoi(args.distance);
	else
		params->distance = DEFAULT_DISTANCE;
	if (args.temperature != NULL)
		params->temperature = atoi(args.temperature);
	else
		params->temperature = DEFAULT_TEMPERATURE;
	if (args.bootstraps != NULL)
		params->bootstraps = atoi(args.bootstraps);
	else
		params->bootstraps = DEFAULT_BOOTSTRAPS;
	params->device = 0;
	params->multi = args.multi;

	// Make sure files exist and are readable
	if (params->validate() == false) {
		cout << "Input validation failed!" << endl;
		exit(EXIT_FAILURE);
	}

	// Start up the algorithm
	IBmrfAlgorithm* bmrf = new BmrfHostAlgorithm(params);

	if (bmrf->loadData() == false)
	{
		cout << "Unable to read input files." << endl;
		exit(EXIT_FAILURE);
	}

	// Do it!
	bmrf->run();

	// Write outputs
	bmrf->writeResults();

	// Clean up
	delete bmrf;
	delete params;

	exit(EXIT_SUCCESS);
}
