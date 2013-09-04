#ifndef BMRF_HOST_ALGORITHM_H_
#define BMRF_HOST_ALGORITHM_H_

#include "BmrfParams.hpp"
#include "IBmrfAlgorithm.h"

#include "hdf5.h"
#include "matvec.h"
#include "Timer.h"

/**
 * @class BmrfHostAlgorithm
 * @brief BMRF implementation using CPU
 */
class BmrfHostAlgorithm : public IBmrfAlgorithm {
public:

	BmrfHostAlgorithm(BmrfParams* params);
	virtual ~BmrfHostAlgorithm();

	virtual bool loadData();
	virtual bool run();
	bool runHost();
	bool runHostMulti();
	virtual bool writeResults();

public:

    Timer* timer;

};

#endif
