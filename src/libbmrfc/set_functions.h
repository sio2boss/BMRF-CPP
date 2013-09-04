/*
 * set_functions.h
 *
 *  Created on: Feb 8, 2013
 *      Author: sio2
 */

#ifndef SET_FUNCTIONS_H_
#define SET_FUNCTIONS_H_

#include "matvec.h"

/**
 * Return indices of ppi in col where gid == ppi value
 */
void findi(Uint2dArray &ppi, unsigned int col, Uint1dArray &gid,
		Uint1dArray &index1);

/**
 * Return indices of ppi in col where gid == ppi value
 */
long findi(Uint1dArray &ppi, unsigned int &val);

/**
 * Return indices of ppi in col where gid == ppi value
 */
void findv(Uint2dArray &ppi, unsigned int col, Uint1dArray &gid,
		Uint1dArray &values);

/**
 * Return indices of ppi from both columns
 */
void finddoublev(Uint2dArray &ppi, Uint1dArray &gid, Uint1dArray &values);

/**
 *
 */
void uniqueSet(Uint1dArray &a, Uint1dArray &u_array);


/**
 * if sorted
 */
void uniqueSetInplaceMustBeSorted(Uint1dArray &a);

/**
 * Only return values in both lists
 */
void intersectSets(Uint1dArray &index1, Uint1dArray &index2,
		Uint1dArray &index3);

/**
 * Only return values in both lists
 */
void intersectSetsAlreadyUnique(Uint1dArray &index1, Uint1dArray &index2,
		Uint1dArray &index3);

/**
 * Only return values in both lists
 */
void intersectSetsReturnIndex(Uint1dArray &set1, Uint1dArray &set2,
		Uint1dArray &set3);

/**
 * Only return values in both lists, both lists must be unique
 */
void intersectSetsReturnIndexAlreadyUnique(Uint1dArray &set1, Uint1dArray &set2,
		Uint1dArray &set3);

/**
 * Bubble sort set
 */
void sortSet(Uint1dArray &to_be_sorted);

/**
 * Add lists and filter duplicates
 */
void unionSets(Uint1dArray &index1, Uint1dArray &index2, Uint1dArray &index3);

/**
 * Return all elementes in a that are not in b
 */
void diffSets(Uint1dArray &a, Uint1dArray &b, Uint1dArray &d);

#endif /* SET_FUNCTIONS_H_ */
