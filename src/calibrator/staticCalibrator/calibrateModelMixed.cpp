#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>

#include "StaticCalibratorCommandLineInterface.h"
#include "../GetRunTimeMixed.h"

//result associates average time to fire reaction over all
//replications with a string that keeps track of the SSA used
struct result{
  std::string ssaName;
  double averageTimePerStep;
};
bool resultCmpr(result a, result b)
{
  return (a.averageTimePerStep < b.averageTimePerStep);
}

//helper functions for writing to file
std::string writeResultsToString(std::vector<result> results);
std::string createModelCalibrationDirectory(std::string modelFileName);
void createTimeCalibrationFile(std::string modelCalibrationDirectory, double simulationTime, std::string resultsString, bool overwritePreviousResults);


int main(int argc, char* argv[])
{
  //get necessary command line options
  STOCHKIT::StaticCalibratorCommandLineInterface commandLineInterface(argc,argv);
  std::string modelFileName = commandLineInterface.getModelFileName();
  double simulationTime = commandLineInterface.getSimulationTime();
  int nReplications = commandLineInterface.getReplications();

  //get total running times for each method over all replications
  //remember getRunTime functions returns average time to fire a reaction
  double odmTime=0.0, nrmTime=0.0, constTime=0.0;
  for(int i=0; i<nReplications; i++)
    {
      odmTime += getODMRunTimeMixed((char*)modelFileName.c_str(), simulationTime);
      nrmTime += getNRMRunTimeMixed((char*)modelFileName.c_str(), simulationTime);
      constTime += getConstantRunTimeMixed((char*)modelFileName.c_str(), simulationTime);
    }
  std::cout << "Simulations Complete!" << std::endl;

  //average total running times for each method over the number
  //of replications to allow for random noise in timing data
  std::vector<result> results;
  results.resize(3);
  results[0].ssaName="ODM";   results[0].averageTimePerStep=odmTime / nReplications;
  results[1].ssaName="NRM";   results[1].averageTimePerStep=nrmTime / nReplications;
  results[2].ssaName="CONST"; results[2].averageTimePerStep=constTime / nReplications;
  std::sort(results.begin(), results.end(), resultCmpr);
  std::cout << "Optimal Method Found!" << std::endl;

  //write out optimal methods in order along with time to string
  std::string resultsString = writeResultsToString(results);

  //create calibration directory if not created and write txt file
  std::string resultsDirectoryName = createModelCalibrationDirectory(modelFileName);
  bool overwritePreviousResults = commandLineInterface.getForce();
  createTimeCalibrationFile(resultsDirectoryName,simulationTime,resultsString,overwritePreviousResults);

  return 0;
}



//HELPER FUNCTIONS

std::string writeResultsToString(std::vector<result> results)
{
  //declare final string to be returned
  std::stringstream ret;
  ret << results[0].ssaName << "\t" << results[0].averageTimePerStep << std::endl;
  ret << results[1].ssaName << "\t" << results[1].averageTimePerStep << std::endl;
  ret << results[2].ssaName << "\t" << results[2].averageTimePerStep << std::endl;

  return ret.str();
}

//create the model directory to save results txt file and return name
//of directory as string
std::string createModelCalibrationDirectory(std::string modelFileName)
{
  std::string modelCalibrationDirectory = modelFileName;

  //truncate ".xml" from model file then append to get output directory
  std::size_t dotXMLSubstringLocation = modelCalibrationDirectory.find_last_of(".");
  modelCalibrationDirectory = modelCalibrationDirectory.substr(0,dotXMLSubstringLocation);
  modelCalibrationDirectory += "_CalibrationInfo";

  //if CalibrationInfo directory has not been created before then create
  try{
    if ( !(boost::filesystem::exists(modelCalibrationDirectory)) ) 
      boost::filesystem::create_directories(modelCalibrationDirectory);
  }
  catch (...) {
    std::cerr << "StochKit ERROR (staticCalibrator::createModelCalibrationDirectory): error creating output directory.\n";
    exit(1);
  }

  return modelCalibrationDirectory;
}

void createTimeCalibrationFile(std::string modelCalibrationDirectory, double simulationTime, std::string resultsString, bool overwritePreviousResults)
{
  //convert simulation time to string as file name
  std::string calibrationFileName = modelCalibrationDirectory;
  calibrationFileName += "/calib";
  std::ostringstream oss;
  oss << simulationTime;
  calibrationFileName += oss.str();
  calibrationFileName += ".txt";

  try{
    //check to see if calibration has already been done on this model
    if (boost::filesystem::exists(calibrationFileName))
      {
	//check to see if the user want to re-calibrate despite
	//already having calibrated the model in the past
	if (!overwritePreviousResults)
	  {
	    std::cout << "StochKit ERROR (staticCalibrator::createTimeCalibrationFile): calibration info already exists for this model and simulation time.\n";
	    std::cout << "Run with --force or -f to overwrite.\n";
	    std::cout << "Calibration terminated.\n";
	    exit(1);
	}
	else
	  {
	    //delete existing file
	    if( remove(calibrationFileName.c_str()) != 0 )
	      std::cerr << "StochKit ERROR (staticCalibrator::createTimeCalibrationFile): could not remove old calibration file\n";
	  }
      }

    //create file
    std::ofstream calibrationFile;
    calibrationFile.open( calibrationFileName.c_str() );
    
    if ( !calibrationFile )
      {
	std::cerr << "StochKit ERROR (staticCalibrator::createTimeCalibrationFile): Unable to open output file. Terminating.\n";
	exit(1);
      }

    calibrationFile << resultsString;
    calibrationFile.close();
    std::cout << "Calibration Data Saved!" << std::endl;
  }
  catch (...) {
    std::cerr << "StochKit ERROR (staticCalibrator::createTimeCalibrationFile): error creating output directory.\n";
    exit(1);
  }
    
}
