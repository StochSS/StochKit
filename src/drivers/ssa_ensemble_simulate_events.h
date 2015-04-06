#ifndef _SSA_ENSEMBLE_SIMULATE_EVENTS_H_
#define _SSA_ENSEMBLE_SIMULATE_EVENTS_H_

#include <vector>
#include "StandardDriverTypes.h"

namespace STOCHKIT
{

/**
    @brief ensemble simulate function for models with events

    -since the actual currentPopulation of the solvers are protected
    members of the class we must let the eventhandler update a local
    population vector then use the setCurrentPopulation( populationVector )
    method of the solver to update the solver's actual currentPopulation
    -be sure that whenever a trigger fires and thus changes the state
    of the sytems, that we update the propensities of the system and
    resample a new reaction to fire
*/
template<typename _SSASolverType, typename _EventHandlerType, typename _OutputType>
bool ensemble_simulate_events(_SSASolverType &solver, _EventHandlerType& eventHandler, std::size_t realizations, double startTime, double endTime, _OutputType &output)
{
    double currentTime = startTime, currentTimeStep = 0.0;

    typedef typename _SSASolverType::populationVectorType populationVectorType;
    //typedef typename _SSASolverType::propensitiesVectorType propensitiesVectorType;
    populationVectorType currentPopulation;
    //propensitiesVectorType currentPropensities;

    solver.prepare(realizations, startTime, endTime);

    currentPopulation = solver.getInitialPopulation();

    if (!output.initialize(realizations,startTime,endTime,currentPopulation)) {
        std::cerr << "StochKit ERROR (ensembel_simulate): initialization of output object failed, simulation aborted\n";
        return false;
    }

    std::vector<double> outputTimes = output.getOutputTimes();
    std::size_t totalIntervals=outputTimes.size();

    std::size_t currentInterval;

    eventHandler.initialize(startTime,endTime);

    for (std::size_t currentRealization=0; currentRealization!=realizations; ++currentRealization) {
        solver.resetParameters();

        solver.initialize(startTime);
        eventHandler.reset();
        currentTime = startTime;
        currentTimeStep = 0.0;
        currentInterval=0;
        currentPopulation = solver.getCurrentPopulation();
        currentTimeStep = solver.selectStepSize();

        //Check to see if system has already met conditions for state based trigger
        eventHandler.fireStateBasedTriggerEvents(currentTime,currentPopulation);
        //Check to see if any time based events will occur before next reaction fires
        while (eventHandler.nextTriggerTime() < currentTime+currentTimeStep)
        {
            while (currentInterval<totalIntervals && eventHandler.nextTriggerTime() >outputTimes[currentInterval])
            {
                output.record(currentRealization,currentInterval,currentPopulation);
                currentInterval++;
            }

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
            while (currentInterval<totalIntervals && currentTime >=outputTimes[currentInterval])
            {
                output.record(currentRealization,currentInterval,currentPopulation);
                currentInterval++;
            }

            //the local var currentPopulation must be updated to match the solver's
            //currentPopulation after firing a reaction within the scope of the solver
            solver.fireReaction();
            currentPopulation = solver.getCurrentPopulation();

            eventHandler.fireStateBasedTriggerEvents(currentTime,currentPopulation);

            currentTimeStep = solver.selectStepSize();

            while (eventHandler.nextTriggerTime() < currentTime+currentTimeStep) {
                currentTime=eventHandler.nextTriggerTime();

                while (currentInterval<totalIntervals && currentTime>outputTimes[currentInterval])
                {
                    //note only recording if currentTime strictly > output time
                    output.record(currentRealization,currentInterval,currentPopulation);
                    currentInterval++;
                }

                eventHandler.fireTimeBasedTriggerEvents(currentTime,currentPopulation);

                eventHandler.fireStateBasedTriggerEvents(currentTime,currentPopulation);

                //this will catch the case where a time-based trigger event occurs exactly on an output time
                while (currentInterval<totalIntervals && currentTime>=outputTimes[currentInterval])
                {
                    output.record(currentRealization,currentInterval,currentPopulation);
                    currentInterval++;
                }

                currentTimeStep=solver.selectStepSize();
            }

            currentTime += currentTimeStep;
        }

        while (currentInterval<totalIntervals && currentTime>=outputTimes[currentInterval])
        {
            currentPopulation = solver.getCurrentPopulation();
            output.record(currentRealization,currentInterval,currentPopulation);
            currentInterval++;
        }

    }

    return true;
}

}

#endif
