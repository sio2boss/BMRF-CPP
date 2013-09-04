#ifndef BMRF_PARAMS_H_
#define BMRF_PARAMS_H_

#include <string>
#include <iostream>

typedef struct {
	std::string ppiFilename;
	std::string genesFilename;
	std::string seedsFilename;
	std::string outputFilename;
	int device;
	int distance;
	int temperature;
	int bootstraps;
	bool multi;

	bool validate() {
		std::cout << "Bmrf Parameters:" << std::endl;
		std::cout << "  PPI file             : " << this->ppiFilename << std::endl;
		std::cout << "  Gene Expression file : " << this->genesFilename
				<< std::endl;
		std::cout << "  Seed Genes file      : " << this->seedsFilename << std::endl;
		std::cout << "  Output file          : " << this->outputFilename << std::endl;
		std::cout << "  Device               : " << this->device << std::endl;
		std::cout << "  Distance             = " << this->distance << std::endl;
		std::cout << "  Temperature          = " << this->temperature << std::endl;
		std::cout << "  Bootstraps           = " << this->bootstraps << std::endl;
		std::cout << "  Multiple Threads     = " << std::boolalpha << this->multi << std::endl;

		return true;
	}
	;

} BmrfParams;

#endif

