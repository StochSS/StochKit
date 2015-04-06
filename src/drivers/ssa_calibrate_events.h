#ifndef _SSA_CALIBRATE_EVENTS_H_
#define _SSA_CALIBRATE_EVENTS_H_

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
template<typename _SSASolverType, typename _EventHandlerType, typename _OutputType>
bool calibrate_events(_SSASolverType &solver, _EventHandlerType& eventHandler, std::size_t realizations, double startTime, double endTime, _OutputType &output)
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


    eventHandler.initialize(startTime,endTime);
    solver.resetParameters();
    solver.initialize(startTime);
    currentTime = startTime;
    currentTimeStep = 0.0;
    currentPopulation = solver.getCurrentPopulation();

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

    //Check to see if system has already met conditions for state based trigger
    eventHandler.fireStateBasedTriggerEvents(currentTime,currentPopulation);
    //Check to see if any time based events will occur before next reaction fires
    while (eventHandler.nextTriggerTime() < currentTime+currentTimeStep)
    {

        //if so update the time of the next time-trigger
        currentTime=eventHandler.nextTriggerTime();

        //actually fire that time-trigger
        eventHandler.fireTimeBasedTriggerEvents(currentTime,currentPopulation);

        //see if firing that time trigger leads to statebased trigger needing to be fired
        eventHandler.fireStateBasedTriggerEvents(currentTime,currentPopulation);

       //resample the reaction we will fire
        currentTimeStep = solver.selectStepSize();
    }
    currentTime += currentTimeStep;

    //Main While Loop
    while (currentTime<endTime)
    {
        //the local var currentPopulation must be updated to match the solver's
        //currentPopulation after firing a reaction within the scope of the solver
        solver.fireReaction();
        nReactionsFired++;
        currentPopulation = solver.getCurrentPopulation();

        eventHandler.fireStateBasedTriggerEvents(currentTime,currentPopulation);

        currentTimeStep = solver.selectStepSize();

        while (eventHandler.nextTriggerTime() < currentTime+currentTimeStep) {
             currentTime=eventHandler.nextTriggerTime();

            eventHandler.fireTimeBasedTriggerEvents(currentTime,currentPopulation);

            eventHandler.fireStateBasedTriggerEvents(currentTime,currentPopulation);

            currentTimeStep=solver.selectStepSize();
        }

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
