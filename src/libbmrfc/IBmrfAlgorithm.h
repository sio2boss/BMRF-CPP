#ifndef I_BMRF_ALGORITHM_H_
#define I_BMRF_ALGORITHM_H_

#include "BmrfParams.hpp"

#include "matvec.h"


class IBmrfAlgorithm {
public:

	IBmrfAlgorithm(BmrfParams* params): params(params) {}
	virtual ~IBmrfAlgorithm() {}

	virtual bool loadData() = 0;
	virtual bool run() = 0;
	virtual bool writeResults() = 0;

public:

	BmrfParams* params;

	Float2dArray* ppiArray;
	Float2dArray* geneArray;
	Float2dArray* geneIdArray;
	Float2dArray* geneLabelArray;
	Float2dArray* seedgeneIdArray;

	Float2dArray* bmrfNetworkScoreArray;
	Float2dArray* bmrfNetworkLengthArray;
	Float2dArray* bmrfNetworkIdArray;
};

#endif
