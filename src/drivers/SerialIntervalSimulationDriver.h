#ifndef _SERIAL_INTERVAL_SIMULATION_DRIVER_H_
#define _SERIAL_INTERVAL_SIMULATION_DRIVER_H_

#include "Input_mass_action.h"
#include "CommandLineInterface.h"
#include "StandardDriverTypes.h"
#include "StandardDriverUtilities.h"
#include "Input_mixed_before_compile.h"
#ifdef MIXED
  #include "Input_mixed_after_compile.h"
#endif
#ifdef EVENTS
  #include "Input_events.h"
  #include "Input_events_after_compile.h"
#endif
#include "Input_tag.h"
#include "ModelTag.h"
#include <sstream>
#include <cstdio>
#include <string>
#include <cstring>
#include <vector>

#include "CommandPassAux.h"

namespace STOCHKIT
{

  /**
      @class SerialIntervalSimulationDriver

      @brief creates solver, simulates using solver and simulate function, and writes output

      -when creating solver an input object either of type Input_mass_action, Input_events_after_compile
      or Input_mixed_after_compile will be used to parse a model file and output initial populations,
      stoichiometries, and etc. in a usable form that a solver can take in

  */
template<typename _solverType, typename _eventHandlerType, typename _outputType>
class SerialIntervalSimulationDriver
{
	
public:
  typedef typename _solverType::populationVectorType populationVectorType;
  typedef typename _solverType::stoichiometryType stoichiometryType;
  typedef typename _solverType::propensitiesType propensitiesType;
  typedef typename _solverType::dependencyGraphType dependencyGraphType;

#ifdef EVENTS
  SerialIntervalSimulationDriver(int ac, char* av[], bool simulate(_solverType&, _eventHandlerType&, std::size_t, double, double, _outputType&) ):
    commandLine(ac,av),
    output(),
    simulateFunction(simulate)
  {}
#else
  SerialIntervalSimulationDriver(int ac, char* av[], bool simulate(_solverType&, std::size_t, double, double, _outputType&) ):
    commandLine(ac,av),
    output(),
    simulateFunction(simulate)
  {}
#endif

  /**
    @brief returns a newly created mass action solver

    -uses Input_mass_action object to generate data structures such as initialPopulations,
    stochiometry, and etc. for the solver to take in

    @see Input_mass_action 
  */
  _solverType createMassActionSolver() {
    
#ifdef WIN32//for visual studio
	std::string mystring;
	mystring=commandLine.getModelFileName();
	char *const modelFileName = static_cast<char *>(_alloca(mystring.length() + 2048) );
	memset(modelFileName, 0, mystring.length() + 2048);
	strcpy(modelFileName,mystring.c_str());
#else
	char modelFileName[commandLine.getModelFileName().length() + 2048];
	memset(modelFileName, 0, commandLine.getModelFileName().length() + 2048);
        strcpy(modelFileName,commandLine.getModelFileName().c_str());
#endif

    Input_mass_action<populationVectorType, stoichiometryType, propensitiesType, dependencyGraphType> model(modelFileName);

    _solverType solver(model.writeInitialPopulation(),
		       model.writeStoichiometry(),
		       model.writePropensities(),
		       model.writeDependencyGraph()
#ifdef JACOBIAN
				,
				model.writeReactants(),
				model.writeRates()
#endif
#ifdef SPATIAL
				,
				commandLine.getDataFileName()
#endif
				);

    if (commandLine.getUseSeed()) {
      solver.seed(commandLine.getSeed());
    }

    return solver;
  }
  
  //Input_mixed_after_compile won't compile unless CustomPropensityFunctions.h is included (i.e. we're compiling a mixed model)
#ifdef MIXED
  /**
    @brief returns a newly created mixed solver

    -uses Input_mixed_after_compile object to generate data structures such as initialPopulations,
    stochiometry, and etc. for the solver to take in

    @see Input_mixed_after_compile 
  */
  _solverType createMixedSolver() {
    
#ifdef WIN32
	std::string name;
	name=commandLine.getModelFileName();
	char *const modelFileName = static_cast<char *>(_alloca(name.length() + 2048) );
	memset(modelFileName, 0, name.length() + 2048);
        strcpy(modelFileName,name.c_str());
#else
	char modelFileName[commandLine.getModelFileName().length() + 2048];
	memset(modelFileName, 0, commandLine.getModelFileName().length() + 2048);
	strcpy(modelFileName,commandLine.getModelFileName().c_str());
#endif
    
    Input_mixed_after_compile<populationVectorType, stoichiometryType, propensitiesType, dependencyGraphType> model(modelFileName);

    _solverType solver(model.writeInitialPopulation(),
		       model.writeStoichiometry(),
		       model.writePropensities(),
		       model.writeDependencyGraph());

    if (commandLine.getUseSeed()) {
      solver.seed(commandLine.getSeed());
    }

    return solver;
  }

#endif

#ifdef EVENTS
  /**
    @brief returns a newly created mixed solver

    -uses Input_events_after_compile object to generate data structures such as initialPopulations,
    stochiometry, and etc. for the solver to take in

    @see Input_events_after_compile
  */
  _solverType createEventsSolver() {
    
#ifdef WIN32
	std::string name;
	name=commandLine.getModelFileName();
	char *const modelFileName = static_cast<char *>(_alloca(name.length() + 2048) );
	memset(modelFileName, 0, name.length() + 2048);
        strcpy(modelFileName,name.c_str());
#else
	char modelFileName[commandLine.getModelFileName().length() + 2048];
	memset(modelFileName, 0, commandLine.getModelFileName().length() + 2048);
	strcpy(modelFileName,commandLine.getModelFileName().c_str());
#endif
	
  //Arya 7/23/13
	//typedef StandardEventHandler<StandardDriverTypes::populationType> eventsType;
    //Input_events_after_compile<populationVectorType, stoichiometryType, propensitiesType, dependencyGraphType, eventsType, _solverType> model(modelFileName);
  Input_events_after_compile<populationVectorType, stoichiometryType, propensitiesType, dependencyGraphType, _eventHandlerType, _solverType> model(modelFileName);

    _solverType solver(model.writeInitialPopulation(),
		       model.writeStoichiometry(),
		       model.writePropensities(),
		       model.writeDependencyGraph());

    if (commandLine.getUseSeed()) {
      solver.seed(commandLine.getSeed());
    }

    return solver;
  }
#endif

  void callSimulate(_solverType& solver) {
    std::size_t realizations=commandLine.getRealizations();
    double simulationTime=commandLine.getSimulationTime();
    
    std::size_t intervals=commandLine.getIntervals();
    
    //set output options
    output.setOutputTimes(IntervalOutput<populationVectorType>::createUniformOutputTimes(0.0,simulationTime,intervals));
    output.setKeepStats(commandLine.getKeepStats());
    output.setKeepTrajectories(commandLine.getKeepTrajectories());
    output.setKeepHistograms(commandLine.getKeepHistograms());
#ifdef USER_OUTPUT
    output.setKeepUserOutput(commandLine.getKeepUserOutput());
#endif
    output.setHistogramBins(commandLine.getHistogramBins());

    if (commandLine.getSpeciesSubset().size()!=0) {
      output.setSpeciesSubset(commandLine.getSpeciesSubset());
    }

#ifdef EVENTS
    //Arya 7/23/13
    //typedef StandardEventHandler<StandardDriverTypes::populationType> eventsType;
    //unfortunately, we need to create another instance of model here

#ifdef WIN32
	std::string name;
	name=commandLine.getModelFileName();
	char *const modelFileName = static_cast<char *>(_alloca(name.length() + 2048) );
	memset(modelFileName, 0, name.length() + 2048);
        strcpy(modelFileName,name.c_str());
#else
	char modelFileName[commandLine.getModelFileName().length() + 2048];
	memset(modelFileName, 0, commandLine.getModelFileName().length() + 2048);
	strcpy(modelFileName,commandLine.getModelFileName().c_str());
#endif

    //Have the Input_events_after_compile parse the XML file for events then pass that info to the eventsHandler who in turn is passed to the simulateFunction
    //Arya 7/23/13
    Input_events_after_compile<populationVectorType, stoichiometryType, propensitiesType, dependencyGraphType, _eventHandlerType, _solverType> model(modelFileName);
    _eventHandlerType eventsHandler=model.writeEvents(solver);
    
    //solver.template simulateEvents<_outputType>(realizations, 0.0, simulationTime, output, eventsHandler);
    simulateFunction(solver, eventsHandler, realizations, 0.0, simulationTime, output);

#else
    simulateFunction(solver, realizations, 0.0, simulationTime, output);
#endif
  }

  void writeOutput() {
    if (!commandLine.getUseExistingOutputDirs()) {
      StandardDriverUtilities::createOutputDirs(commandLine,false);
    }

    if (commandLine.getKeepStats()) {     
      output.stats.writeMeansToFile(commandLine.getOutputDir()+"/"+commandLine.getStatsDir()+"/"+commandLine.getMeansFileName());
      output.stats.writeVariancesToFile(commandLine.getOutputDir()+"/"+commandLine.getStatsDir()+"/"+commandLine.getVariancesFileName());
      output.stats.writeSimulationInfoFile(commandLine.getOutputDir()+"/"+commandLine.getStatsDir()+"/"+commandLine.getStatsInfoFileName());
    }
    if (commandLine.getKeepTrajectories()) {
      std::size_t trajectoryNumber;
      std::string trajectoryNumberString;
      for (std::size_t i=0; i!=commandLine.getRealizations(); ++i) {
	trajectoryNumber=i+commandLine.getTrajectoriesOffset();
	trajectoryNumberString=StandardDriverUtilities::size_t2string(trajectoryNumber);

	if (commandLine.getLabel()) {
	  IntervalOutput<StandardDriverTypes::populationType>::writeLabelsToFile(commandLine.getOutputDir()+"/"+commandLine.getTrajectoriesDir()+"/trajectory"+trajectoryNumberString+".txt",commandLine.getSpeciesNames());
	  output.trajectories.writeDataToFile(i,commandLine.getOutputDir()+"/"+commandLine.getTrajectoriesDir()+"/trajectory"+trajectoryNumberString+".txt",true,true);
	}
	else {  
	  output.trajectories.writeDataToFile(i,commandLine.getOutputDir()+"/"+commandLine.getTrajectoriesDir()+"/trajectory"+trajectoryNumberString+".txt");
	}
      }
    }
    if (commandLine.getKeepHistograms()) {
      output.histograms.writeHistogramsToFile(commandLine.getOutputDir()+"/"+commandLine.getHistogramsDir()+"/hist",".dat",commandLine.getSpeciesNames());
      //thread 0 will write a histograms info file
      if (commandLine.getHistogramsInfoFileName()!="") {
	std::ofstream outfile;
	outfile.open((commandLine.getOutputDir()+"/"+commandLine.getHistogramsDir()+"/"+commandLine.getHistogramsInfoFileName()).c_str());
	if (!outfile) {
	  std::cerr << "StochKit ERROR (SerialIntervalSimulationDriver::writeOutput): Unable to open histogram info file for writing.\n";
	  exitFunc(1);
	}
	
	outfile << output.histograms.numberOfSpecies() << "\n";
	
	outfile.close();
      }
    }
#ifdef USER_OUTPUT
   if(commandLine.getKeepUserOutput()){
       output.user_output.writeDataToFile(commandLine.getOutputDir()+"/"+commandLine.getUserOutputDir()+"/"+commandLine.getUserOutputFileName());
   }
#endif
  }

#ifdef SPATIAL
		void writeOutput(std::size_t NumberOfSubvolumes) {
			if (!commandLine.getUseExistingOutputDirs()) {
				StandardDriverUtilities::createOutputDirs(commandLine,false);
			}

			if (commandLine.getKeepStats()) {     
				output.stats.writeMeansToFile(commandLine.getOutputDir()+"/"+commandLine.getStatsDir()+"/"+commandLine.getMeansFileName());
				output.stats.writeVariancesToFile(commandLine.getOutputDir()+"/"+commandLine.getStatsDir()+"/"+commandLine.getVariancesFileName());
				output.stats.writeSimulationInfoFile(commandLine.getOutputDir()+"/"+commandLine.getStatsDir()+"/"+commandLine.getStatsInfoFileName());
			}
			if (commandLine.getKeepTrajectories()) {
				std::size_t trajectoryNumber;
				std::string trajectoryNumberString;
				for (std::size_t i=0; i!=commandLine.getRealizations(); ++i) {
					trajectoryNumber=i+commandLine.getTrajectoriesOffset();
					trajectoryNumberString=StandardDriverUtilities::size_t2string(trajectoryNumber);

					if (commandLine.getLabel()) {
						IntervalOutput<StandardDriverTypes::populationType>::writeLabelsToFile(commandLine.getOutputDir()+"/"+commandLine.getTrajectoriesDir()+"/trajectory"+trajectoryNumberString+".txt",commandLine.getSpeciesNames());
						output.trajectories.writeDataToFile(i,commandLine.getOutputDir()+"/"+commandLine.getTrajectoriesDir()+"/trajectory"+trajectoryNumberString+".txt",true,true);
					}
					else {  
						output.trajectories.writeDataToFile(i,commandLine.getOutputDir()+"/"+commandLine.getTrajectoriesDir()+"/trajectory"+trajectoryNumberString+".txt");
					}
				}
			}
			if (commandLine.getKeepHistograms()) {
				std::vector<std::string> species=commandLine.getSpeciesNames();
				std::size_t NumberOfSpecies=species.size();
				std::vector<std::string> labels(NumberOfSpecies*NumberOfSubvolumes);
				std::stringstream ss;
				ss.clear();
				std::size_t k=0;
				for(std::size_t i=0; i<NumberOfSpecies; i++)
				{
					for(std::size_t j=0; j<NumberOfSubvolumes; j++)
					{
						ss<<species[i]<<"_"<<j;
						labels[k]=ss.str();
						ss.str("");
						k++;
					}
				}
				output.histograms.writeHistogramsToFile(commandLine.getOutputDir()+"/"+commandLine.getHistogramsDir()+"/hist",".dat",labels);
				//thread 0 will write a histograms info file
				if (commandLine.getHistogramsInfoFileName()!="") {
					std::ofstream outfile;
					outfile.open((commandLine.getOutputDir()+"/"+commandLine.getHistogramsDir()+"/"+commandLine.getHistogramsInfoFileName()).c_str());
					if (!outfile) {
						std::cerr << "StochKit ERROR (SerialIntervalSimulationDriver::writeOutput): Unable to open histogram info file for writing.\n";
						exit(1);
					}

					outfile << output.histograms.numberOfSpecies() << "\n";

					outfile.close();
				}
			}
		}

		void writeOutput(std::size_t NumberOfRows, std::size_t NumberOfColumns) {
			if (!commandLine.getUseExistingOutputDirs()) {
				StandardDriverUtilities::createOutputDirs(commandLine,false);
			}

			if (commandLine.getKeepStats()) {     
				output.stats.writeMeansToFile(commandLine.getOutputDir()+"/"+commandLine.getStatsDir()+"/"+commandLine.getMeansFileName());
				output.stats.writeVariancesToFile(commandLine.getOutputDir()+"/"+commandLine.getStatsDir()+"/"+commandLine.getVariancesFileName());
				output.stats.writeSimulationInfoFile(commandLine.getOutputDir()+"/"+commandLine.getStatsDir()+"/"+commandLine.getStatsInfoFileName());
			}
			if (commandLine.getKeepTrajectories()) {
				std::size_t trajectoryNumber;
				std::string trajectoryNumberString;
				for (std::size_t i=0; i!=commandLine.getRealizations(); ++i) {
					trajectoryNumber=i+commandLine.getTrajectoriesOffset();
					trajectoryNumberString=StandardDriverUtilities::size_t2string(trajectoryNumber);

					if (commandLine.getLabel()) {
						IntervalOutput<StandardDriverTypes::populationType>::writeLabelsToFile(commandLine.getOutputDir()+"/"+commandLine.getTrajectoriesDir()+"/trajectory"+trajectoryNumberString+".txt",commandLine.getSpeciesNames());
						output.trajectories.writeDataToFile(i,commandLine.getOutputDir()+"/"+commandLine.getTrajectoriesDir()+"/trajectory"+trajectoryNumberString+".txt",true,true);
					}
					else {  
						output.trajectories.writeDataToFile(i,commandLine.getOutputDir()+"/"+commandLine.getTrajectoriesDir()+"/trajectory"+trajectoryNumberString+".txt");
					}
				}
			}
			if (commandLine.getKeepHistograms()) {
				std::vector<std::string> species=commandLine.getSpeciesNames();
				std::size_t NumberOfSpecies=species.size();
				std::vector<std::string> labels(NumberOfSpecies*NumberOfRows*NumberOfColumns);
				std::stringstream ss;
				ss.clear();
				std::size_t l=0;
				for(std::size_t i=0; i<NumberOfSpecies; i++)
				{
					for(std::size_t j=0; j<NumberOfRows; j++)
					{
						for(std::size_t k=0; k<NumberOfColumns; k++)
						{
							ss<<species[i]<<"_"<<j<<k;
							labels[l]=ss.str();
							ss.str("");
							l++;
						}
					}
				}
				output.histograms.writeHistogramsToFile(commandLine.getOutputDir()+"/"+commandLine.getHistogramsDir()+"/hist",".dat",labels);
				//thread 0 will write a histograms info file
				if (commandLine.getHistogramsInfoFileName()!="") {
					std::ofstream outfile;
					outfile.open((commandLine.getOutputDir()+"/"+commandLine.getHistogramsDir()+"/"+commandLine.getHistogramsInfoFileName()).c_str());
					if (!outfile) {
						std::cerr << "StochKit ERROR (SerialIntervalSimulationDriver::writeOutput): Unable to open histogram info file for writing.\n";
						exit(1);
					}

					outfile << output.histograms.numberOfSpecies() << "\n";

					outfile.close();
				}
			}
		}
#endif

  void writeCalibratorOutput()
  {
    if( commandLine.shouldCalibrate() )
      output.writeCalibratorStatToFile( commandLine.getOutputDir()+"/"+commandLine.getCalibratorDir()+"/"+commandLine.getCalibratorFileName() );
  }
 
  bool isMassAction() {
	return commandLine.isMassAction();
  }

  bool isMixed() {
	return commandLine.isMixed();
  }

  bool isEventsEnabled() {
	return commandLine.isEventsEnabled();
  }

  CommandLineInterface getCommandLine() {
	return commandLine;
  }
  
private:
  
	CommandLineInterface commandLine;
	_outputType output;

  #ifdef EVENTS
  bool (*simulateFunction)(_solverType&, _eventHandlerType&, std::size_t, double, double, _outputType&);
  #else
	bool (*simulateFunction)(_solverType&, std::size_t, double, double, _outputType&);
  #endif

};

}

#endif
