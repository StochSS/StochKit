/*
 *  test_SSA_ODM.cpp
 *  
 */
#include<iostream>
#include<string>
#include <fstream>
#include <vector>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include<stdio.h>
#include "SSA_ODM.h"
//#include "SSA_Direct.h"
#include "TestModel.h"
#include "Input_mass_action.h"
#include "StandardDriverTypes.h"
#include "StatsAndTrajectoriesOutput.h"
#include "CustomPropensitySet.h"


using namespace std;

void  chooseSpecies(vector<int> &,int , string );
void  chooseTimes(vector<double> &, double,  double, int );
void  chooseTimes(vector<double> &, string, double, double );
void  chooseTimes(vector<double> &, double);


int main()
{

  typedef boost::numeric::ublas::vector<double> denseVector;
  typedef std::vector<denseVector> denseMatrix;
  typedef boost::numeric::ublas::mapped_vector<double> sparseVector;
  typedef std::vector<sparseVector> sparseMatrix; //dense vector of sparse vectors
  typedef std::vector<std::vector<int> > graphType;

  //denseVector X=TestModelInitialPopulations<denseVector>();//the population vector should always be dense
  ////sparseMatrix NU=TestModelStoichiometry<sparseMatrix>();//sparseMatrix is slower for this small problem
  //denseMatrix NU=TestModelStoichiometry<denseMatrix>();
  //TestModelPropensities<denseVector> PROPENSITIES;  
  //graphType DG=TestModelDependencyGraph<graphType>();

  
  std::string modelname="../models/examples/stochkit_mass_action.xml";
  Input_mass_action<denseVector, StandardDriverTypes::stoichiometryType, CustomPropensitySet<denseVector>, graphType> model(const_cast<char*>(modelname.c_str()));
  denseVector X=model.writeInitialPopulation();//the population vector should always be dense
  denseMatrix NU=model.writeStoichiometry();
  CustomPropensitySet<denseVector> PROPENSITIES;
  graphType DG=model.writeDependencyGraph();

  
  //typedef SSA_Direct<denseVector, denseMatrix, TestModelPropensities<denseVector>, graphType > solverType;
  typedef SSA_ODM<denseVector, denseMatrix, CustomPropensitySet<denseVector>, graphType> solverType;
  solverType solver(X,
				    NU,
  					PROPENSITIES,
 					DG);
/*  solverType solver(model.writeInitialPopulation(),
				    model.writeStoichiometry(),
  					model.writePropensities(),
 					model.writeDependencyGraph());

*/
//  solverType ssa(X,NU,PROPENSITIES,DG, 1);//with dependency graph
  
  double t0=0;
  double tf=1;
  int intervals=5;
  int numRuns =1000;//100000;


  StatsAndTrajectoriesOutput<denseVector> output;
  output.setOutputTimes(IntervalOutput<denseVector>::createUniformOutputTimes(t0,tf,intervals));

  output.setKeepTrajectories(false);

  solver.simulate<StatsAndTrajectoriesOutput<denseVector> >(numRuns,t0,tf,output);

  output.writeMeansToFile("odm_means.txt");
  output.writeVariancesToFile("odm_variances.txt");

  return 0;

}


