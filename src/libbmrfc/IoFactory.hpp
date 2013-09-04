/*
 * IoFactory.hpp
 *
 *  Created on: Jan 25, 2013
 *      Author: sio2
 */

#ifndef IOFACTORY_HPP_
#define IOFACTORY_HPP_

#include "matvec.h"

class IoFactory {
public:

	static bool readArrayFromFile(std::string filename, std::string variable,
			Float2dArray* &array);
	static bool writeArraysToFile(std::string filename, std::string variable,
			Float2dArray* &array, std::string variable2, Float2dArray* &array2,
			std::string variable3, Float2dArray* &array3);

};

#endif /* IOFACTORY_HPP_ */
