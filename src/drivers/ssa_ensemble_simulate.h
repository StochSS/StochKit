#ifndef _SSA_ENSEMBLE_SIMULATE_H_
#define _SSA_ENSEMBLE_SIMULATE_H_

#include <vector>
#include "StandardDriverTypes.h"

namespace STOCHKIT
{

/**
    @brief function which interacts with solver to run a simulation and interacts with output
    to output the simulation results
*/
template<typename _SSASolverType, typename _OutputType>
bool ensemble_simulate(_SSASolverType &solver, std::size_t realizations, double startTime, double endTime, _OutputType &output)
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

    for (std::size_t currentRealization=0; currentRealization!=realizations; ++currentRealization) {
        solver.initialize(startTime);
        currentTime = startTime;
        currentTimeStep = 0.0;
        currentInterval=0;

        //currentPropensities = solver.getCurrentPropensities();

        currentTimeStep = solver.selectStepSize();
        currentTime += currentTimeStep;
        while (currentTime<endTime) {
            while (currentInterval<totalIntervals && currentTime >=outputTimes[currentInterval]){
                currentPopulation = solver.getCurrentPopulation();
                output.record(currentRealization,currentInterval,currentPopulation);
                currentInterval++;
            }

            solver.fireReaction();
            currentTimeStep = solver.selectStepSize();
            currentTime += currentTimeStep;
        }

        while (currentInterval<totalIntervals && currentTime>=outputTimes[currentInterval]){
            currentPopulation = solver.getCurrentPopulation();
            output.record(currentRealization,currentInterval,currentPopulation);
            currentInterval++;
        }

    }

    return true;
}



}

#endif
