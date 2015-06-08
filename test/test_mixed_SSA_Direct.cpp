/*
 *  test_mixed_SSA_Direct.cpp
 *  
 */
#include<iostream>
#include<string>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/thread.hpp>
#include <boost/shared_ptr.hpp>
#include<stdio.h>
#include "SSA_Direct.h"
#include "Input_mixed_before_compile.h"
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


//  Input<denseVector, denseMatrix, CustomPropensitySet<denseVector>, graphType > input("./stochkit.xml");//with dependency graph
  Input_mixed_before_compile<denseVector, denseMatrix, CustomPropensitySet<denseVector>, graphType> input(argv[1]);
  input.writeCustomPropensityFunctionFile("./CustomPropensityFunctions.h");

  //compile
  system("make DM_MIXED_RECOMPILE");
  
  //run recompiled version
  string recompiled_command_1("./test_mixed_SSA_Direct_recompile ");
  string recompiled_command_2(argv[1]);
  string recompiled_command = recompiled_command_1 + recompiled_command_2;

  cout << recompiled_command << endl;

  system(recompiled_command.c_str());

  return 0;
}

