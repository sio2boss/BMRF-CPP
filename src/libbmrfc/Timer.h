#ifndef TIMER_H_
#define TIMER_H_

#include <iostream>
#include <string.h>
#include <queue>
#include <vector>
#include <string>
#include <stdio.h>
#include <stdlib.h>

#include "Logging.h"

/**
 * @class Timer
 * This class is to facilitate timing benchmarks where you have successive calls
 * the each should be evaluated for execution time (such as in any algorithm
 * implementation)
 *
 * Usage should be list this:
 * Timer ctimer;
 * ctimer.start();
 *   Do stuff
 * ctimer.mark("I did stuff");
 *   More stuff that is a lot cooler
 * ctimer.mark("I did more cooler stuff");
 * ctimer.print();
 *
 */
class Timer {

public:

	Timer();
    virtual ~Timer();

    void start();
    void stop();
    void mark(std::string benchmark_point_name);
    void print();

private:

    clock_t cpu_currenttimer;
    std::vector<std::pair<std::string, clock_t> > cpu_timervector;

    LOG4CXX_DECLARE(logger)
};

#endif
