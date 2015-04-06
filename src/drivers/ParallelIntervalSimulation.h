/******************************************************************************
 */

#ifndef _PARALLEL_INTERVAL_SIMULATION_H_
#define _PARALLEL_INTERVAL_SIMULATION_H_

#include "Input_mass_action.h"
#include "CommandLineInterface.h"
#include "StandardDriverTypes.h"
#include "StandardDriverUtilities.h"
#include "Input_mixed_before_compile.h"
#include "Input_events_before_compile.h"
#include "IntervalOutput.h"
#include "StatsOutput.h"
#include "HistogramSingle.h"
#include "boost/thread/thread.hpp"
#include "boost/bind.hpp"
#include "boost/random.hpp"
#include "boost/random/uniform_int.hpp"
#include "boost/filesystem.hpp"
#include <iomanip>
#include <sstream>
#include <ctime>
#include <vector>

#ifdef WIN32
#include <windows.h>
#include <fstream>
#endif

namespace STOCHKIT
{
/**
  @class ParallelIntervalSimulation

  @brief This class is responsible for launching realizations of a simulation in parallel and creating the
  proper output with the help of an output object

  -this class runs executables of ssa_serial.cpp at the system level to minimize the overhead that 
  comes with thread packages such as Boost. This is accomplished the run() function

  @see ssa_serial.cpp
*/
class ParallelIntervalSimulation
{
	
public:
  ParallelIntervalSimulation(int ac, char* av[]);

  void run(std::string executableName);
  
  void runCalibrator(std::string executableName);

  void executable(std::string command);

  std::size_t assignment(std::size_t totalRealizationss, std::size_t threadid);

  /**
    @brief writes the system commands called to launch parallel simulations and saves log files 
    for each thread detailing any errors that occured during simulation

    -this code used to be in mergeouput() but is now a standalone function called in mergeOutput()
    because mergeCalibratorOutput() uses the code as well, and it makes both of the aforementioned
    functions shorter
  */
  void writeCommandLogAndThreadLog();

  void mergeOutput();

  void mergeCalibratorOutput();

  void warnIfLargeOutput();//helper function that prints a warning if simulation will generate a lot of data

  static std::string modifyCmdArgsRealizations(std::string commandLineArguments, std::string subRealizations);
  static std::string modifyCmdArgsSeed(std::string commandLineArguments, std::string newSeed);

protected:
  ParallelIntervalSimulation(int ac, char* av[], int dummy) : commandLine(ac,av), masterProc(0), engine((boost::uint32_t) std::time(0)) { (void)dummy; } //empty constructor; overridden in subclass

  CommandLineInterface commandLine;
  std::size_t numberOfProcesses;

  //records how many realizations were assigned to each processor
  std::vector<std::size_t> processorRealizationAssignments;
  std::size_t masterProc;

  int seedOfSeed;
  boost::mt19937 engine ;
  std::size_t numberOfWorkers;
  boost::thread_group sthreads;

  //initialized in ParallelIntervalSimulation.cpp
  static const size_t nSolvers;
  static const std::string solverNames[3];
};

}

#endif
