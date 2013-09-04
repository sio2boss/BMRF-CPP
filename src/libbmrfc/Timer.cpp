#include "Timer.h"
#include <time.h>

LOG4CXX_INIT(Timer::logger, "Timer");

/**
 * Constructor.  cpu_only argument disables cuda calls
 */
Timer::Timer() {
}

/**
 * Destructor
 */
Timer::~Timer() {
	stop();
}

/**
 * Start recording time
 */
void Timer::start() {

	if (cpu_timervector.empty() == false)
		return;

	std::string name = "Start";
	cpu_currenttimer = clock();
	std::pair<std::string, clock_t> pair(name, cpu_currenttimer);
	cpu_timervector.push_back(pair);

}

/**
 * Record a new time, allows diffing against previous time to determine
 * how long benchmark_point_name took to execute.
 */
void Timer::mark(std::string benchmark_point_name) {
	// Stop current benchmark point
	cpu_currenttimer = clock();
	std::pair<std::string, clock_t> pair(benchmark_point_name,
			cpu_currenttimer);
	cpu_timervector.push_back(pair);
}

void Timer::stop() {
}

void Timer::print() {
	double time;

	// Iterate printing elapsed time as we go
	printf("\nBenchmarkPoint, Time(ms)\n");
	fflush(stdout);

	clock_t a;
	clock_t b;
	std::vector<std::pair<std::string, clock_t> >::iterator i;

	// Get start time
	i = cpu_timervector.begin();
	a = (*i).second;
	++i;

	for (; i != cpu_timervector.end(); ++i) {
		b = (*i).second;
		time = double(b - a) * 1000.0f / double(CLOCKS_PER_SEC);
		printf("%s, %.3f\n", (*i).first.c_str(), time);
		a = b;
	}

	fflush(stdout);

}
