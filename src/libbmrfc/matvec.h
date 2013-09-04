/*
 * matvec.h
 *
 *  Created on: Feb 3, 2013
 *      Author: sio2
 */

#ifndef MATVEC_H_
#define MATVEC_H_

#include "Memory.h"
#include "bmrf_constants.h"

#include <math.h>
#include <cmath>
#include <algorithm>

/**
 *
 */
template <typename T>
struct Array2d {
	T* data;
	unsigned int x;
	unsigned int y;
	const unsigned int x_max;
	const unsigned int y_max;

	Array2d(unsigned int nx, unsigned int ny) :
		x_max(nx), y_max(ny)
	{
		x = nx;
		y = ny;
		data = (T*)malloc(x * y * sizeof(T));
	}

	~Array2d()
	{
		free(data);
	}

	T& operator()(unsigned int i, unsigned int j)
	{
		return data[i * y + j];
	}

	void print(long max_i = -1)
	{
		if (max_i == -1)
			max_i = std::min(MAX_PRINT_LEN, x);

		printf("[");
		for (int i = 0; i < max_i; ++i)
		{
			for (int j = 0; j < y; ++j)
			{
				printf(" %g", float(data[i * y + j]));
			}
			if (i != max_i - 1)
				printf(";\n ");
		}
		if (max_i != x)
			printf("...");
		printf(" ];");
		printf("\n");
		fflush (stdout);
	}

	void shrink(unsigned int new_x, long new_y = -1)
	{
		if (new_y == -1)
		{
			x = new_x;
			return;
		}

		T* new_data;
		new_data = (T*)malloc(new_x * new_y * sizeof(T));
		unsigned int copy_size = new_y * sizeof(T);
		for (int c = 0; c < new_x; c++)
		{
			memcpy(&new_data[c * new_y], &data[c * y], copy_size);
		}
		free(data);
		data = new_data;
		x = new_x;
		y = new_y;
	}

	unsigned int length() const
	{
		return x * y;
	}

	unsigned int size() const
	{
		return length() * sizeof(T);
	}

	Array2d<T>& operator=(const Array2d<T> &rhs) {

		// make sure rhs has dimensions defined
		assert(rhs.x_max == this->x_max);
		assert(rhs.y_max == this->y_max);
		assert(rhs.data != NULL);
		assert(data != NULL);

		// Deal with copy memory there are 4 cases
		memcpy(data, rhs.data, rhs.size());

		return *this;
	}
};

typedef Array2d<unsigned int> Uint2dArray;
typedef Array2d<float> Float2dArray;
typedef Array2d<double> Double2dArray;

/**
 *
 */

template <typename T>
struct Array1d
{
	T* data;
	unsigned int x;
	const unsigned int max;

	Array1d(unsigned int n) :
			max(n)
	{
		x = n;
		data = (T*) malloc(x * sizeof(T));
	}

	~Array1d()
	{
		free(data);
	}

	void print(long max_i = -1)
	{
		if (max_i == -1)
			max_i = std::min(MAX_PRINT_LEN, x);
		printf("[");
		for (int i = 0; i < max_i; ++i)
			printf(" %g", float(data[i]));
		if (max_i < x)
			printf("...");
		printf(" ];");
//		printf("; {x: %i, hid: %i, alloc:%i}", x, x - max_i, this->max);
		printf("\n");
		fflush (stdout);
	}

	void shrink(unsigned int new_size)
	{
		x = new_size;
	}

	T& operator[](unsigned int i)
	{
		return data[i];
	}

	unsigned int length() const
	{
		return x;
	}

	unsigned int size() const
	{
		return length() * sizeof(T);
	}

	bool append(const T &val)
	{
		if (x + 1 > max)
		{
			printf("append error\n");
			return false;
		}

		data[x] = val;
		x++;
		return true;
	}

	Array1d<T>& operator=(const Array1d<T> &rhs) {

		// make sure rhs has dimensions defined
		assert(rhs.max == this->max);
		assert(rhs.data != NULL);
		assert(data != NULL);

		// Deal with copy memory there are 4 cases
		memcpy(data, rhs.data, rhs.size());

		return *this;
	}
};

typedef Array1d<unsigned int> Uint1dArray;
typedef Array1d<float> Float1dArray;
typedef Array1d<double> Double1dArray;

#endif /* MATVEC_H_ */
