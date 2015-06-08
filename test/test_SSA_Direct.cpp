/*
 *  test_SSA_Direct.cpp
 *  
 */
#include<iostream>
#include<string>
#include <fstream>
#include <sstream>
#include <vector>
#include <boost/numeric/ublas/vector.hpp>
#include<stdio.h>
#include "SSA_Direct.h"
#include "TestModel.h"
#include "StatsAndTrajectoriesOutput.h"
//#include "IntervalOutput.h"

using namespace std;

int main()
{
  typedef boost::numeric::ublas::vector<double> denseVector;
  typedef std::vector<denseVector> denseMatrix;
  typedef std::vector<std::vector<int> > graphType;

  denseVector X=TestModelInitialPopulations<denseVector>();
  
  denseMatrix NU=TestModelStoichiometry<denseMatrix>();
  
  TestModelPropensities<denseVector> PROPENSITIES;  
  graphType DG=TestModelDependencyGraph<graphType>();
 
  typedef SSA_Direct<denseVector, denseMatrix, TestModelPropensities<denseVector>, graphType> solverType;
  solverType ssa(X,NU,PROPENSITIES,DG);

  std::size_t numRuns=3000;
  double t0=0;
  double tf=10;
  std::size_t intervals=0;

  //typedef IntervalOutput<denseVector> outputType;
  typedef StatsAndTrajectoriesOutput<denseVector> outputType;
  outputType output;

  output.setOutputTimes(IntervalOutput<denseVector>::createUniformOutputTimes(t0,tf,intervals));

  /* uncomment to keep a species subset
     std::vector<std::size_t> outputSpecies;
     outputSpecies.push_back(1);
     outputSpecies.push_back(2);
     output.setSpeciesSubset(outputSpecies);//keep data for only species 1 and 2 (not 0)
  */

  output.setKeepTrajectories(false);//set this to true to keep trajectories

  std::cout << "running simulation...\n";
  ssa.simulate<outputType>(numRuns,t0,tf,output);
  std::cout << "finished running simulation...writing to file...\n";

//		gettimeofday(&tv, NULL);
//		elapsed = (tv.tv_sec - startedat_s) * 1000000 + (tv.tv_usec - startedat_us);
//		std::cout << "simulation time : " << elapsed/1000000  << "of " <<numRuns<<" runs" <<"\n";

  //output.writeMeansToFile("dm_means.txt");
  //output.writeVariancesToFile("dm_variances.txt");

  //write trajectories to file, if keeping trajectories data
  if (output.keepTrajectories) {
    for (std::size_t i=0; i!=numRuns; ++i) {
      std::string numSuffix;
      std::stringstream num;
      num << i;
      numSuffix=num.str();
      std::string filename="traj"+numSuffix+".txt";
      output.writeTrajectoryToFile(i,filename);
    }
  }

  std::cout << "done!\n";

  return 0;
}
