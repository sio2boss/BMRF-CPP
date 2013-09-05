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

or run manually from ./matlab directory

    ../build/bmrf-cpu --ppi sim/ppi.mat --genes sim/genes.mat --seeds sim/seed_gene_ids.mat --output asdf.mat --multi

Additionally, a slighly modified version of BMRF (matlab release: http://www.cbil.ece.vt.edu/software.htm) can be run with

	demoMatlab

## Using EC2

This repo includes portions of https://github.com/sio2boss/VagrantSeed which will allow you to run easily on EC2.  Follow the "Setting up Vagrant on your machine" step on the VagrantSeed github site.  Or you could setup your own Ubuntu machine and use the provision script in ./devops

