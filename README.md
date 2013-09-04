# BMRF-CPP

Bagging Markov Random Field (BMRF) developed by Li Chen, Jianhua Xuan, Rebecca B. Riggins, Yue Wang, and Robert Clarke is a new method of identifying differentially expressed subnetworks from protein-protein interaction (PPI) networks. Published November 17 2012 in Nucleic Acids Research.

This is the CPP implementation with 1D and 2D parallelism.  Intel compiler is optional.

See doc/Barnes_RO_T_2013.pdf for thesis.

## Build

You will need a few dependencies: HDF5, log4cxx, Boost, threadpool, cmake

Simply running the following will run cmake and build the code

    make

## Running

Open MATLAB and cd to ./matlab

	demoCpu

Additionally, a slighly modified version of BMRF (matlab release: http://www.cbil.ece.vt.edu/software.htm) can be run with

	demoMatlab


