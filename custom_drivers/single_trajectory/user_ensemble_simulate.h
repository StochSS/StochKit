#ifndef _USER_ENSEMBLE_SIMULATE_H_
#define _USER_ENSEMBLE_SIMULATE_H_

#include <vector>
#include <algorithm>
#include "StandardDriverTypes.h"

namespace STOCHKIT
{

/**
    @brief function which interacts with solver to run a simulation and interacts with output
    to output the simulation results
           single trajectory simulation, keep trajectory, max reaction counts = min(1e6,int_limit)
*/
template<typename _SSASolverType, typename _OutputType>
bool ensemble_simulate(_SSASolverType &solver, std::size_t realizations, double startTime, double endTime, _OutputType &output)
{
    double currentTime = startTime, currentTimeStep = 0.0;

    typedef typename _SSASolverType::populationVectorType populationVectorType;
    //typedef typename _SSASolverType::propensitiesVectorType propensitiesVectorType;
    populationVectorType currentPopulation;
    //propensitiesVectorType currentPropensities;

    std::size_t max_steps = std::min((std::size_t)100000,std::numeric_limits<std::size_t>::max());
    std::size_t step_count = 0;
    std::vector<double> outputTime(0);
    std::vector<populationVectorType> outputPop(0);

    if (realizations != 1){
        std::cout << "StochKit MESSAGE (ssa_single_traj_simulate): realizations != 1, may cause segmentation fault or generate huge amount of output\n";
    }

    solver.prepare(realizations, startTime, endTime);

    currentPopulation = solver.getInitialPopulation();


    solver.initialize(startTime);
    currentTime = startTime;
    currentTimeStep = 0.0;

    outputTime.push_back(currentTime);
    outputPop.push_back(currentPopulation);

    currentTimeStep = solver.selectStepSize();
    currentTime += currentTimeStep;
    while (currentTime<endTime && step_count < max_steps) {
        solver.fireReaction();

        currentPopulation = solver.getCurrentPopulation();
        outputTime.push_back(currentTime);
        outputPop.push_back(currentPopulation);

        currentTimeStep = solver.selectStepSize();
        currentTime += currentTimeStep;

        ++step_count;
    }

    if(output.getKeepStats() != false || output.getKeepTrajectories() != true || output.getKeepHistograms() != false){
        std::cout << "StochKit MESSAGE (ssa_single_traj_simulate): recommand output options --no-stats --keep-trajectories\n";
    }

    output.setOutputTimes(outputTime);
    if (!output.initialize(1,startTime,endTime,currentPopulation)) {
        std::cerr << "StochKit ERROR (ssa_single_traj_simulate): initialization of output object failed, simulation aborted\n";
        return false;
    }

    for(std::size_t currentInterval=0;currentInterval<outputTime.size();++currentInterval){
        output.record(0,currentInterval,outputPop[currentInterval]);
    }

    return true;
}



}

#endif
