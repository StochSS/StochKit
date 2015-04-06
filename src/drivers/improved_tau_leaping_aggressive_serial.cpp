/*
 *  improved_tau_leaping_serial.cpp
 *  
 */
#include <iostream>
#include <string>
#include "boost_headers.h"
#include "StandardEventHandler.h"
#include "StandardDriverTypes.h"
#include "SerialIntervalSimulationDriver.h"
#include "improved_tau_leaping_aggressive.h"

using namespace STOCHKIT;

template<typename _solverType, typename _outputType>
bool simulate(_solverType& solver, std::size_t realizations, double startTime, double endTime, _outputType& output)
{
	solver.template simulate<_outputType>(realizations, startTime, endTime, output);
	return true;
};

int main(int ac, char* av[])
{
  typedef improvedTauLeaping_aggressive<StandardDriverTypes::populationType,
    StandardDriverTypes::stoichiometryType, 
    StandardDriverTypes::propensitiesType,
    StandardDriverTypes::graphType> solverType;
  typedef StandardEventHandler<StandardDriverTypes::populationType> eventHandlerType;
  typedef StandardDriverTypes::outputType outputType;

  SerialIntervalSimulationDriver<solverType, eventHandlerType, outputType> driver(ac,av,simulate<solverType, outputType>);

  solverType solver=driver.createMassActionSolver();
  
  //set solver-specific parameters
  solver.setEpsilon(driver.getCommandLine().getEpsilon());
  solver.setThreshold(driver.getCommandLine().getThreshold());
  solver.setSSASteps(driver.getCommandLine().getSSASteps());

  driver.callSimulate(solver);

  driver.writeOutput();

  return 0;
}
