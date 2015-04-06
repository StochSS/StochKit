/*
 *  test_mass_action_SSA_Direct.cpp
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
#include "Input_mass_action.h"
#include "StatsAndTrajectoriesOutput.h"

using namespace std;

int main(int argc, char *argv[])
{
  if(argc!=2){
     cout << "Please use the format test_mass_action_Direct_SSA filename" << endl;
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
  
  boost::shared_ptr<SimpleMessageHandler> messageHandler(new SimpleMessageHandler());

//  Input<denseVector, denseMatrix, CustomPropensitySet<denseVector>, graphType > input("./stochkit.xml");//with dependency graph
  Input_mass_action<denseVector, denseMatrix, CustomPropensitySet<denseVector>, graphType, SimpleMessageHandler> input(argv[1], messageHandler);//with dependency graph

//  std::cout << "input pass" << std::endl;

  X = input.writeInitialPopulation();
//  cout << "initial population:" << endl;
//  for(unsigned int i=0; i< X.size(); ++i)
//	  cout << " " << X[i];
//  cout << endl;

  NU = input.writeStoichiometry();
//  cout << "stoichiometry matrix:" << endl;
//  for(unsigned int i=0; i<NU.size(); ++i){
//	  for(unsigned int j=0; j<NU[i].size(); ++j)
//		  cout << " " << NU[i][j];
//	  cout << endl;
//  }
  
  DG = input.writeDependencyGraph();
//  cout << "dependency graph:" << endl;
//  for(unsigned int i=0; i<DG.size(); ++i){
//	  for(unsigned int j=0; j<DG[i].size(); ++j)
//		  cout << " " << DG[i][j];
//	  cout << endl;
//  }

  PROPENSITIES = input.writePropensities();

  typedef SSA_Direct<denseVector, denseMatrix, CustomPropensitySet<denseVector>, graphType, SimpleMessageHandler> solverType;
  solverType ssa(X,NU,PROPENSITIES,DG,messageHandler);

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
  StatsAndTrajectoriesOutput<denseVector,SimpleMessageHandler> output(messageHandler);
  output.setOutputTimes(IntervalOutput<denseVector,SimpleMessageHandler>::createUniformOutputTimes(t0,tf,intervals));

  /* uncomment to keep a species subset
     std::vector<std::size_t> outputSpecies;
     outputSpecies.push_back(1);
     outputSpecies.push_back(2);
     output.setSpeciesSubset(outputSpecies);//keep data for only species 1 and 2 (not 0)
  */

  output.setKeepTrajectories(false);//set this to true to keep trajectories

  messageHandler->handleMessage("running simulation...");
  ssa.simulate<StatsAndTrajectoriesOutput<denseVector,SimpleMessageHandler> >(numRuns,t0,tf,output);
  messageHandler->handleMessage("finished running simulation...writing to file...");
  
  output.writeMeansToFile("dm_means.txt");
  output.writeVariancesToFile("dm_variances.txt");

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

  messageHandler->handleMessage("done!");

  return 0;
}

