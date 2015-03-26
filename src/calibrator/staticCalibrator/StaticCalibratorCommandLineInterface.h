/******************************************************************************
 */

#ifndef _STATIC_CALIB_COMMAND_LINE_INTERFACE_H_
#define _STATIC_CALIB_COMMAND_LINE_INTERFACE_H_

#include <iostream>
#include <string>
#include <vector>
#include <exception>
#include "boost/program_options.hpp"
#include "boost/filesystem.hpp"
#include "Input_tag.h"
#include "ModelTag.h"
#include "DenseVectorSubset.h"

namespace STOCHKIT
{
 class StaticCalibratorCommandLineInterface
 {
	
 public:
  
   //CONSTRUCTOR
   StaticCalibratorCommandLineInterface(int ac, char* av[]);

   //GETTERS FOR COMMAND LINE OPTION VALUES
   std::string getModelFileName() const;
   
   double getSimulationTime() const;

   std::size_t getReplications() const;

   bool getUseSeed() const;
   int getSeed() const;

   bool getForce() const;

   bool getRecompile() const;

   std::string getCmdArgs() const;

   //GETTERS FOR HIDDEN COMMAND LINE OPTIONS
   bool getUseExistingOutputDirs() const;
   
   std::string getStatsDir() const;
   std::string getMeansFileName() const;
   std::string getVariancesFileName() const;
   std::string getStatsInfoFileName() const;
   std::string getTrajectoriesDir() const;
   std::size_t getTrajectoriesOffset() const;
   std::string getHistogramsDir() const;
   std::string getHistogramsInfoFileName() const;
   
   std::size_t getSSASteps() const;

   std::string getGeneratedCodeDir() const;
   
   
   
 protected:
   void parse(int ac, char* av[]);


 private:
    
  boost::program_options::variables_map vm;
  boost::program_options::options_description visible;
  boost::program_options::options_description hidden;
  boost::program_options::options_description combined;
  
  std::string modelFileName;
  double simulationTime;
  std::size_t replications;
  int seed;
  bool useSeed;
  bool recompile;
  std::string cmdArgs;
  bool force;
  std::string outputDir;

  //HIDDEN OPTION VALUES
  bool useExistingOutputDirs;
  std::string statsDir;
  std::string meansFileName;
  std::string variancesFileName;
  std::string statsInfoFileName;
  std::string trajectoriesDir;
  std::size_t trajectoriesOffset;
  std::string histogramsDir;
  std::string histogramsInfoFileName;
  std::string generatedCodeDir;
  std::size_t SSASteps;
  
 };
}

#endif
