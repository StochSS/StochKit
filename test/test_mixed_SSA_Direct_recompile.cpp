/*
 *  test_mixed_SSA_Direct.cpp
 *  
 */
#include<iostream>
#include<string>
#include <fstream>
#include <vector>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/thread.hpp>
#include <boost/shared_ptr.hpp>
#include<stdio.h>
#include "SSA_Direct.h"
#include "Input_mixed_after_compile.h"
#include "StatsAndTrajectoriesOutput.h"

using namespace std;

int main(int argc, char *argv[])
{
  if(argc!=2){
     cout << "Please use the format test_mixed_Direct_SSA filename" << endl;
     exit(1);
  }

  typedef boost::numeric::ublas::vector<double> denseVector;
  typedef std::vector<denseVector> denseMatrix;
  typedef boost::numeric::ublas::mapped_vector<double> sparseVector;
  typedef std::vector<sparseVector> sparseMatrix; //dense vector of sparse vectors
  typedef std::vector<std::vector<int> > graphType;

  denseVector X;//the population vector should always be dense
  
  //sparseMatrix NU=TestModelStoichiometry<sparseMatrix>();//sparseMatrix is slower for this small problem
  denseMatrix NU;
  
//  TestModelPropensities<denseVector> PROPENSITIES;  
  CustomPropensitySet<denseVector> PROPENSITIES;
  graphType DG;
  

//  Input<denseVector, denseMatrix, CustomPropensitySet<denseVector>, graphType > input("./stochkit.xml");//with dependency graph
  Input_mixed_after_compile<denseVector, denseMatrix, CustomPropensitySet<denseVector>, graphType> input(argv[1]);

  X = input.writeInitialPopulation();
  NU = input.writeStoichiometry();
  DG = input.writeDependencyGraph();
  PROPENSITIES = input.writePropensities();

  typedef SSA_Direct<denseVector, denseMatrix, CustomPropensitySet<denseVector>, graphType> solverType;
  solverType ssa(X,NU,PROPENSITIES,DG);

/*
  //SSA_Direct<denseVector, sparseMatrix, TestModelPropensities<denseVector> > ssa(X,NU,PROPENSITIES);//create solver object without dependency graph
  SSA_Direct<denseVector, denseMatrix, CustomPropensitySet<denseVector>, graphType > ssa(X,NU,PROPENSITIES,DG);//with dependency graph
  //note that the solver will detect that the problem is small and not use the dependency graph even if it exists
*/

  std::size_t numRuns=100000;
  double t0=0;
  double tf=1;
  std::size_t intervals=10;
  
  //typedef 
  StatsAndTrajectoriesOutput<denseVector> output;
  output.setOutputTimes(IntervalOutput<denseVector>::createUniformOutputTimes(t0,tf,intervals));

  /* uncomment to keep a species subset
     std::vector<std::size_t> outputSpecies;
     outputSpecies.push_back(1);
     outputSpecies.push_back(2);
     output.setSpeciesSubset(outputSpecies);//keep data for only species 1 and 2 (not 0)
  */

  output.setKeepTrajectories(false);//set this to true to keep trajectories

  std::cout << "running simulation...\n";
  ssa.simulate<StatsAndTrajectoriesOutput<denseVector> >(numRuns,t0,tf,output);
  std::cout << "finished running simulation...writing to file...\n";
  
  output.writeMeansToFile("mixed_dm_means.txt");
  output.writeVariancesToFile("mixed_dm_variances.txt");

  std::cout << "done!\n";

  return 0;
}

