#ifndef MEMORY_H_
#define MEMORY_H_

#include <cstdio>
#include <stdlib.h>
#include <cstring>

template<typename T> void allocateCPUMemory(T*& cpu, std::size_t size_bytes)
{
	cpu = (T*)malloc(size_bytes);
}

template<typename T> void freeCPUMemory(T*& cpu)
{
	free(cpu);
	cpu = NULL;
}

// Define a CPU vs GPU type of memory
enum MemoryType
{
	CPU = 1
};

// Now we can mask the four different memory transfer types
template<typename T> void copy(T* &a, T* const b, const MemoryType &a_type,
		const MemoryType &b_type, const unsigned int &size)
{
	
		memcpy(a, b, size);
}

// Wrap allocate and free using type
template<typename T> void allocateMemory(T*& data, std::size_t size_bytes,
		const MemoryType &type)
{
		allocateCPUMemory(data, size_bytes);
}

template<typename T> void freeMemory(T*& data, const MemoryType &type)
{
		freeCPUMemory(data);
	data = NULL;
}

template<typename T> void wipeMemory(T*& data, std::size_t size_bytes, const MemoryType &type, int initval =
		0)
{
		std::memset((void*) data, initval, size_bytes);
}

#endif
