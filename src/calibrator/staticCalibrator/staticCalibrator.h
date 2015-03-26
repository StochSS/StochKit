#ifndef _STATIC_CALIBRATOR_H_
#define _STATIC_CALIBRATOR_H_

#include <stdlib.h>
#include <cstdio>
#include <sstream>

#include "StaticCalibratorCommandLineInterface.h"
#include "ModelTag.h"
#include "Input_mixed_before_compile.h"

#include "../GetRunTime.h"

//result associates average time to fire reaction over all
//replications with a string that keeps track of the SSA used
struct result{
  std::string ssaName;
  double averageTimePerStep;
};
bool resultCmpr(result a, result b);

class staticCalibrator
{
 public:
  //constructor
  staticCalibrator(int argc, char* argv[]);

  //destructor deletes commandLineInterface
  ~staticCalibrator();

 private:
  //FUNCTIONS
  
  //performs necessary simulations, save timing data
  //to appropriate variables, set optimal method, and
  //write results to file
  void performSimulations();

  //same as performSimulations but for mixed models
  void performSimulationsMixed();

  //creates output directory to store calibration info for different simulation times
  void createModelCalibrationDirectory();
  void createTimeCalibrationFile();

  void setResultString();

  //VARIABLES
  
  //pointer to commandLineInteface object that's used to parse command line arguments.
  //this is declared in the constructor and deallocated in the destructor of staticCalibrator
  STOCHKIT::StaticCalibratorCommandLineInterface* commandLineInterface;

  //simulation information
  std::string modelFileName;
  std::string modelCalibrationDirectory;

  double simulationTime;
  std::size_t nReplications;

  //average time to fire a reaction for each method
  std::vector<result> results;

  std::string resultString;
};

#endif
