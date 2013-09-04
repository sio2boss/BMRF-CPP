/*
 * IoFactory.cpp
 *
 *  Created on: Jan 25, 2013
 *      Author: sio2
 */

#include "IoFactory.hpp"

#include "MatlabFile.hpp"

/**
 *
 */
bool IoFactory::readArrayFromFile(std::string filename, std::string variable,
		Float2dArray* &array) {

	MatlabFile f(READ);
	f.open(filename);
	long status = f.readArray(variable, array);
	f.close();

	return (status == 0);
}

/**
 *
 */
bool IoFactory::writeArraysToFile(std::string filename, std::string variable,
		Float2dArray* &array, std::string variable2, Float2dArray* &array2,
		std::string variable3, Float2dArray* &array3) {

	long status;
	MatlabFile f(WRITE);

	f.open(filename);
	status = f.writeArray(variable, array);
	status = f.writeArray(variable2, array2);
	status = f.writeArray(variable3, array3);

	f.close();
	return (status == 0);
}
