/*
 *  FILE:    UserOutput.h
 */

#ifndef _USER_OUTPUT_H_
#define _USER_OUTPUT_H_

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <limits>
#include <iomanip>
#include "boost/foreach.hpp"
#include "boost/tokenizer.hpp"

namespace STOCHKIT
{
 template<typename _populationVectorType>
 class UserOutput
 {	
/*
 * Data structre to record outputs, need to include at least
 * 1. data
 * 2. output time points
 */
  protected:
  // define your own data structure here
  std::vector< std::vector<_populationVectorType> > data;

  // vector of output times
  std::vector<double> outputTimes;
  

/*
 * Interfaces to interact with drivers, need to inlucde at least
 * 1. constructor/destructor
 * 2. virtual bool initialize(std::size_t realizations, double startTime, double endTime, _populationVectorType& samplePopulationVector)
 * 3. virtual void record(std::size_t realization, std::size_t interval,_populationVectorType population)
 * 4. static UserOutput createFromFiles(std::string dataFileName)
 * 5. virtual void merge(UserOutput other)
 * 6. void writeDataToFile(std::string filename)
 * 7. std::vector<double> getOutputTimes()
 * 8. void setOutputTimes(std::vector<double> outputTimes)
 *
 * the function and parameters of each of the functions will be commented below
 *
 */
  public:

/*
 * Function: 
 *    default constructor
 * Parameters:
 *    realizations: the number of realizations
 * Return Value:
 *    none
 */
  UserOutput(std::size_t realizations=0) : data(realizations)
  {}

/*
 * Function: 
 *    default destructor
 * Parameters:
 *    none
 * Return Value:
 *    none
 */
  virtual ~UserOutput() {
  }

/*
 * Function:
 *    initialize output object
 * Parameters:
 *    realizations: number of realizations on that thread
 *    startTime: start time of simulation
 *    endTime: end time of simulation
 *    samplePopulationVector: initial population vector
 * Return value:
 *    true: initialization succeeded
 *    false: initialization failed
 */
  virtual bool initialize(std::size_t realizations, double /*startTime*/, double endTime, _populationVectorType& samplePopulationVector) {
    if (outputTimes.size()==0) {
      std::vector<double> defaultOutputTimes;
      defaultOutputTimes.push_back(endTime);
      setOutputTimes(defaultOutputTimes);
    }
    //need to add checks for consistency. e.g. output times are increasing order, never less than start time or greater than end time

    data.clear();
    resize(realizations,outputTimes.size(),samplePopulationVector.size());
    return true;
  }

/*
 * Function:
 *    return the vector of output times
 * Parameters:
 *    none
 * Return Value:
 *    the vector of output times
 */ 
  std::vector<double> getOutputTimes() {
    return outputTimes;
  }

/*
 * Function:
 *   set the vector of output times
 * Parameters:
 *   outputTimes: the vector of output times
 * Return Value:
 *   none
 */ 
  void setOutputTimes(std::vector<double> outputTimes){
    // in this case we only set one outputTime t_end no matter what -i command option specify
    this->outputTimes=outputTimes;
    if(this->outputTimes.size()!=1){
        this->outputTimes[0]=this->outputTimes.back();
        this->outputTimes.resize(1);
    }
  }

/*
 * Function:
 *   record output at the end of each step
 * Parameters:
 *   realization: the index of realization that is being currently simulated
 *   interval: the index of output time point that is currently being written
 *   population: the current state of simulated system
 * Return Value:
 *   none
 */ 
  virtual void record(std::size_t realization, std::size_t /*interval*/,_populationVectorType population) {
    //to achieve best performance, this function does no consistency checking
    // rewrite the same vector again and again so only end time population will be recorded
    data[realization][0]=population;
  }

/*
 * Function:
 *   create output object from an existing output data file
 * Parameters:
 *   dataFileName: the name of the existing output data file
 * Return Value:
 *   the output object
 */ 
   static UserOutput createFromFiles(std::string dataFileName) {
    std::ifstream dataIn(dataFileName.c_str());
    if (!dataIn) {
      std::cerr << "StochKit ERROR (UserOutput::createFromFiles): Unable to open user output data file "<< dataFileName<<".\n";
      exit(1);
    }
    

    //read lines from output file
    // first line is number of realizations
    std::size_t numberOfRealizations;
    dataIn >> numberOfRealizations;
    // second line is number of species
    std::size_t numberOfSpecies;
    dataIn >> numberOfSpecies;
    // third line the end of simulation time
    double endTime; 
    dataIn >> endTime;

    UserOutput userOutput;

    std::vector<double> defaultOutputTimes;
    defaultOutputTimes.push_back(endTime);
    userOutput.setOutputTimes(defaultOutputTimes);

    userOutput.data.clear();
    userOutput.resize(numberOfRealizations,defaultOutputTimes.size(),numberOfSpecies);
    
    boost::char_separator<char> sep("\t","\n");

    // each line represents a realization
    // all the elements are population of each species at the end of simulation
    std::string line;
    std::size_t tokenCounter;
    _populationVectorType popRow(numberOfSpecies);
    typename _populationVectorType::value_type popVal;

    std::size_t realization = 0;

    // first eat the '\n' left by dataIn >> endTime
    std::getline(dataIn,line);

    while (std::getline(dataIn,line)) {
      boost::tokenizer<boost::char_separator<char> > tokens(line,sep);
      tokenCounter=0;
      BOOST_FOREACH(std::string t, tokens) {
	std::istringstream ss(t);
	ss >> popVal;
	popRow.insert_element(tokenCounter,popVal);
	tokenCounter++;
      }
      //store popRow in data
      userOutput.data[realization][0]=popRow;
      realization++;
    }

    if (realization != numberOfRealizations) {
      std::cerr << "StochKit ERROR (UserOutput::createFromFiles): number of realizations != number of rows of data.\n";
      exit(1);
    }
    return userOutput;

   }

/*
 * Function:
 *   merge another output object to current output object (for combining results from parallel simulation)
 * Parameters:
 *   other: the other output object to be combined with this one
 * Return Value:
 *   none
 */ 
  virtual void merge(UserOutput other) {
    data.insert(data.end(),other.data.begin(),other.data.end());
  }

/*
 * Function:
 *   merge a vector of other output objects to current output object (for combining results from parallel simulation)
 * Parameters:
 *   others: the vector of other output objects to be combined with this one
 * Return Value:
 *   none
 */ 
  virtual void merge(std::vector<UserOutput> others) {
    for (std::size_t i=0; i!=others.size(); ++i) {
      this->merge(others[i]);
    }
  }

/*
 * Function:
 *   write output data to a data file
 * Parameters:
 *   filename: the name of the output data file
 * Return Value:
 *   none
 */ 
  void writeDataToFile(std::string filename, bool printTime=false, bool append=false, bool highPrecision=false) {

    std::ofstream outfile;

    if (append) {
      outfile.open(filename.c_str(),std::ios::out | std::ios::app);
    }
    else {
      outfile.open(filename.c_str());
    }

    if (!outfile) {
      std::cout << "StochKit ERROR (UserOutput::writeDataToFile): Unable to open output file.\n";
      exit(1);
    }
    try {
      // first line is number of realizations
      outfile << data.size() << std::endl;
      // second line is number of species
      if(data.size() != 0 && data[0].size()!=0){
        outfile << data[0][0].size() << std::endl;
      } else {
        outfile << "0" << std::endl;
      }
      // third line the end of simulation time
      if(outputTimes.size()!=0){
        outfile << outputTimes.back() << std::endl;
      } else {
        outfile << "0" << std::endl;
      }
      for(std::size_t realization=0; realization!=data.size(); ++realization){
        for (std::size_t interval=0; interval!=outputTimes.size(); ++interval) {
          if (printTime) {
	    outfile << outputTimes[interval] << "\t";
	  }
	  for (size_t index=0; index!=data[realization][interval].size(); ++index) {
	    if (highPrecision) {
	      outfile << std::setprecision(std::numeric_limits<double>::digits10)<< data[realization][interval][index] << "\t";
	    }
	    else {
	      outfile << std::setprecision(8) << data[realization][interval][index] << "\t";
	    }
	  }
	  outfile << "\n";
        }
      }
      outfile.close();
    }
    catch (...) {
      std::cout << "StochKit ERROR (UserOutput::writeDataToFile): error writing data to output file.\n";
      exit(1);
    }
  }
 
/*
 * Local functions
 */ 
  protected:
  void resize(std::size_t realizations, std::size_t intervals, std::size_t populationVectorSubsetSize) {
    //beware off-by-1 error
    //doesn't change outputTimes
    data.resize(realizations);
    for (std::size_t i=0; i!=realizations; ++i) {
      data[i].resize(intervals, _populationVectorType(populationVectorSubsetSize));
    }
  }

 };
}
#endif
  
