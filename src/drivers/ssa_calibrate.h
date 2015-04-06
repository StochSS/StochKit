#ifndef _SSA_CALIBRATE_H_
#define _SSA_CALIBRATE_H_

#include <vector>
#include "StandardDriverTypes.h"

#include <utility>
#ifdef WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

namespace STOCHKIT
{

/**
    @brief function which interacts with solver to calibrate and interacts with output
    to output the calibrator results
*/
template<typename _SSASolverType, typename _OutputType>
bool calibrate(_SSASolverType &solver, std::size_t realizations, double startTime, double endTime, _OutputType &output)
{
    double currentTime = startTime, currentTimeStep = 0.0;

    typedef typename _SSASolverType::populationVectorType populationVectorType;
    //typedef typename _SSASolverType::propensitiesVectorType propensitiesVectorType;
    populationVectorType currentPopulation;
    //propensitiesVectorType currentPropensities;

    solver.prepare(realizations, startTime, endTime);

    currentPopulation = solver.getInitialPopulation();

    if (!output.initialize(realizations,startTime,endTime,currentPopulation)) {
        std::cerr << "StochKit ERROR (calibrate): initialization of output object failed, simulation aborted\n";
        return false;
    }

    solver.initialize(startTime);
    currentTime = startTime;
    currentTimeStep = 0.0;

    unsigned nReactionsFired=0;

    //initalize timing variables
#ifdef WIN32
	LARGE_INTEGER begin, end, freq;
	QueryPerformanceFrequency(&freq);
	QueryPerformanceCounter(&begin);
#else
    timeval timer;
    gettimeofday(&timer, NULL);
    double start_time=timer.tv_sec+(timer.tv_usec/1000000.0);
#endif

    //start simulation
    currentTimeStep = solver.selectStepSize();
    currentTime += currentTimeStep;
    while (currentTime<endTime) {
        solver.fireReaction();
        nReactionsFired++;
        currentTimeStep = solver.selectStepSize();
        currentTime += currentTimeStep;
    }

    //calculate total time and calibrator stat
#ifdef WIN32
	QueryPerformanceCounter(&end);
	double totalTime=(double)(end.QuadPart-begin.QuadPart)/freq.QuadPart;
#else
    gettimeofday(&timer,NULL);
    double end_time=timer.tv_sec+(timer.tv_usec/1000000.0);
    double totalTime=(end_time-start_time);//elapsed time, units=seconds
#endif
    double calibratorStat = totalTime / ((double)nReactionsFired);
    output.recordCalibratorStat(calibratorStat);

    return true;
}



}

#endif
