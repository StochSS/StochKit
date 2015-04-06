/*
 *  test_SSA_Direct.cpp
 *  
 */
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/thread.hpp>
#include <boost/shared_ptr.hpp>
#include <stdio.h>
#include <iomanip>
#include "SSA_Direct.h"
#include "TestModel.h"
#include "IntervalOutput.h"
#include "StatsOutput.h"
#include "DenseVectorSubset.h"
#include "Random.h"

using namespace std;

extern bool veryClose(double,double,double);

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

  typedef StatsOutput<denseVector> outputType;
  outputType output;

  std::size_t numRuns=100;
  double t0=0;
  double tf=1;
  std::size_t intervals=5;

  output.setOutputTimes(IntervalOutput<denseVector>::createUniformOutputTimes(t0,tf,intervals));

  ssa.simulate<outputType>(numRuns,t0,tf,output);

  output.writeSimulationInfoFile(".createFromFiles_info.txt");
  output.writeMeansToFile("testMeans.txt",true,false,true);
  output.writeVariancesToFile("testVariances.txt",true,false,true);

  //StatsOutput<denseVector>::createFromFiles("models/examples/stochkit_mass_action_output/stats/means.txt","junk",".createFromFiles_info");
  outputType out2=StatsOutput<denseVector>::createFromFiles("testMeans.txt","testVariances.txt",".createFromFiles_info.txt");

  out2.writeMeansToFile("testMeans2.txt");
  out2.writeVariancesToFile("testVariances2.txt",true,false,true);

  //outputType x2=StatsOutput<denseVector>::createFromFiles(x1);

  return 0;
}

bool veryClose(double x, double y,double eps=1E-13) {
  return fabs(x-y)<=eps;
}
