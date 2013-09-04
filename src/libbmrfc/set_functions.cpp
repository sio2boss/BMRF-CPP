#include "set_functions.h"

/**
 * Return indices of ppi in col where gid == ppi value
 */
void findi(Uint2dArray &ppi, unsigned int col, Uint1dArray &gid,
		Uint1dArray &index1) {
	unsigned int size_index1 = 0;
	for (int id = 0; id < gid.x; ++id) {
		for (int i = 0; i < ppi.x; ++i) {
			if (ppi(i, col) == gid[id]) {
				index1[size_index1] = i;
				size_index1 += 1;
			}
		}
	}

	index1.shrink(size_index1);
}

/**
 * Return indices of ppi in col where gid == ppi value
 */
long findi(Uint1dArray &ppi, unsigned int &val) {
	for (int i = 0; i < ppi.x; ++i) {
		if (ppi[i] == val) {
			return i;
		}
	}
	return -1;
}

/**
 * Return corresponding values of ppi where gid == ppi value
 */
void findv(Uint2dArray &ppi, unsigned int col, Uint1dArray &gid,
		Uint1dArray &values) {
	unsigned int size_values = 0;
	for (int i = 0; i < ppi.x; ++i) {
		for (int id = 0; id < gid.x; ++id) {
			if (ppi(i, col) == gid[id]) {
				values[size_values] = ppi(i, ppi.y - col - 1);
				size_values += 1;
			}
		}
	}

	values.shrink(size_values);
}

/**
 * Return corresponding values of ppi where gid == ppi value
 */
void finddoublev(Uint2dArray &ppi, Uint1dArray &gid, Uint1dArray &values) {
	unsigned int size_values = 0;
	for (int i = 0; i < ppi.x; ++i) {
		for (int id = 0; id < gid.x; ++id) {
			if (ppi(i, 0) == gid[id]) {
				values[size_values] = ppi(i, 1);
				size_values += 1;
			}
			if (ppi(i, 1) == gid[id]) {
				values[size_values] = ppi(i, 0);
				size_values += 1;
			}
		}
	}

	values.shrink(size_values);
}

/**
 * if sorted
 */
void uniqueSetAlreadySorted(Uint1dArray &a, Uint1dArray &u_array) {
	// Add first element
	int size_u_array = 1;
	u_array[0] = a[0];

	// Add others based on sorted order
	for (int i = 1; i < a.x; ++i) {
		if (a[i] != a[i - 1]) {
			u_array[size_u_array] = a[i];
			size_u_array += 1;
		}
	}
	u_array.shrink(size_u_array);
}

/**
 *
 */
void uniqueSet(Uint1dArray &a, Uint1dArray &u_array) {

//	printf("uniqueSet_size, %i\n", a.x);

// Add first element
	unsigned int size_u_array = 1;
	u_array[0] = a[0];

	// Add others if unsorted
	for (int i = 1; i < a.x; ++i) {
		bool isUnique = true;
		for (int j = 0; j < size_u_array; ++j) {
			if (a[i] == u_array[j]) {
				isUnique = false;
				break;
			}
		}
		if (isUnique == true) {
			u_array[size_u_array] = a[i];
			size_u_array += 1;
		}
	}
	u_array.shrink(size_u_array);
}

/**
 * uniqueSetInplaceMustBeSorted
 */
void uniqueSetInplaceMustBeSorted(Uint1dArray &a) {

	// One element array is already unique
	if (a.x <= 1)
		return;

	// Assume two element array
	unsigned int size_u_array = 0;
	if (a[1] != a[0]) {
		++size_u_array;
	}

	if (a.x == 2) {
		a.shrink(size_u_array+1);
		return;
	}

	// Add others if unsorted
	for (unsigned int i = 2; i < a.x; ++i) {
		if (a[i] != a[size_u_array]) {
			++size_u_array;
			a[size_u_array] = a[i];
		}
	}
		a.shrink(size_u_array+1);
}

/**
 * Only return values in both lists
 */
void intersectSets(Uint1dArray &index1, Uint1dArray &index2,
		Uint1dArray &index3) {

//	printf("intersectSets_size, %i, %i\n", index1.x, index2.x);
	Uint1dArray new_index1(index1.x);
	Uint1dArray new_index2(index2.x);
	unsigned int size_new_index1 = 0;
	unsigned int size_new_index2 = 0;
	uniqueSet(index1, new_index1);
	uniqueSet(index2, new_index2);
//	printf("new_index1: ");
//	new_index1.print(size_new_index1);
//	printf("new_index2: ");
//	new_index2.print(size_new_index2);

	unsigned int size_index3 = 0;
	for (int i = 0; i < new_index1.x; ++i) {
		for (int j = 0; j < new_index2.x; ++j) {
			if (new_index1[i] == new_index2[j]) {
				index3[size_index3] = new_index1[i];
				size_index3 += 1;
			}
		}
	}
	index3.shrink(size_index3);
}

/**
 * Only return values in both lists
 */
void intersectSetsAlreadyUnique(Uint1dArray &index1, Uint1dArray &index2,
		Uint1dArray &index3) {
	unsigned int size_index3 = 0;
	for (int i = 0; i < index1.x; ++i) {
		for (int j = 0; j < index2.x; ++j) {
			if (index1[i] == index2[j]) {
				index3[size_index3] = index1[i];
				size_index3 += 1;
			}
		}
	}
	index3.shrink(size_index3);
}

/**
 * Only return values in both lists
 */
void intersectSetsReturnIndex(Uint1dArray &set1, Uint1dArray &set2,
		Uint1dArray &set3) {

	Uint1dArray new_set1(set1.x);
	Uint1dArray new_set2(set2.x);
	unsigned int size_new_set1 = 0;
	unsigned int size_new_set2 = 0;
	uniqueSet(set1, new_set1);
	uniqueSet(set2, new_set2);
//	printf("new_set1: ");
//	new_set1.print(size_new_set1);
//	printf("new_set2: ");
//	new_set2.print(size_new_set2);

	unsigned int size_set3 = 0;
	for (int i = 0; i < new_set1.x; ++i) {
		for (int j = 0; j < new_set2.x; ++j) {
			if (new_set1[i] == new_set2[j]) {
				set3[size_set3] = i;
				size_set3 += 1;
			}
		}
	}
	set3.shrink(size_set3);
}

/**
 * Only return values in both lists
 */
void intersectSetsReturnIndexAlreadyUnique(Uint1dArray &set1, Uint1dArray &set2,
		Uint1dArray &set3) {

	unsigned int size_set3 = 0;
	for (int i = 0; i < set1.x; ++i) {
		for (int j = 0; j < set2.x; ++j) {
			if (set1[i] == set2[j]) {
				set3[size_set3] = i;
				size_set3 += 1;
			}
		}
	}
	set3.shrink(size_set3);
}

int qsortCompare(const void* a, const void *b) {
	return (*(unsigned int*) a - *(unsigned int*) b);
}

/**
 * Bubble sort set
 */
void sortSet(Uint1dArray &to_be_sorted) {
//	int i, j, flag = 1;
//	int temp;
//	int numLength = to_be_sorted.x;
//	for (i = 1; (i <= numLength) && flag; i++)
//	{
//		flag = 0;
//		for (j = 0; j < (numLength - 1); j++)
//		{
//			// ascending order simply changes to <
//			if (to_be_sorted[j + 1] < to_be_sorted[j])
//			{
//				// swap elements
//				temp = to_be_sorted[j];
//				to_be_sorted[j] = to_be_sorted[j + 1];
//				to_be_sorted[j + 1] = temp;
//				flag = 1;
//			}
//		}
//	}

//	printf("sortSet_size, %i\n", to_be_sorted.x);
	qsort(to_be_sorted.data, to_be_sorted.x, sizeof(unsigned int),
			qsortCompare);
}

/**
 * Add lists and filter duplicates
 */
void unionSets(Uint1dArray &index1, Uint1dArray &index2, Uint1dArray &index3) {

	if (index1.x == 0 && index2.x == 0) {
		index3.shrink(0);
		return;
	}

	// Compact
	memcpy(index3.data, index1.data, index1.x * sizeof(unsigned int));
	memcpy(&index3.data[index1.x], index2.data,
			index2.x * sizeof(unsigned int));
	index3.shrink(index1.x + index2.x);

	// Sort
	sortSet(index3);

	// Remove duplicates
	uniqueSetInplaceMustBeSorted(index3);
}

/**
 * Return all elements in a that are not in b
 */
void diffSets(Uint1dArray &a, Uint1dArray &b, Uint1dArray &d) {

	unsigned int size_index3 = 0;
	for (int i = 0; i < a.x; ++i) {

		// Search b for element in a
		bool b_has_element = false;
		for (int j = 0; j < b.x; ++j) {
			if (a[i] == b[j]) {
				b_has_element = true;
				break;
			}
		}

		// Copy if does not exist
		if (b_has_element == false) {
			d[size_index3] = a[i];
			size_index3 += 1;
		}
	}
	d.shrink(size_index3);
}
