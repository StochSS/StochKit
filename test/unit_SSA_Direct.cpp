/*******************************************************************
 *  FILE: unit_SSA_Direct.cpp
 *
 *  AUTHOR: Sheng
 *
 *  CREATED: Apr 6, 2009
 *
 *  SUMMARY: unit tests for the Direct Method of SSA 
 *           (/src/SSA_Direct.h)
 *
 *  NOTES:
 *       - use DimerDecay model (DimerDecay.h)
 *       - generate an executable file unit_SSA_Direct.o by Makefile
 *
 *  TODO:
 *       1. Initialization tests
 *       2. Single step tests
 *
 ******************************************************************/


#include <iostream>
#include <fstream>
#include <vector>
#include <unistd.h>
#include <assert.h>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "SSA_Direct.h"

// DimerDecay model
#include "DimerDecay.h"
const int NumberOfReactions = 4;
const int NumberOfSpecies = 3;

typedef boost::numeric::ublas::vector<double> denseVector;
typedef std::vector<denseVector> denseMatrix;
typedef boost::numeric::ublas::mapped_vector<double> sparseVector;
typedef std::vector<sparseVector> sparseMatrix;
typedef std::vector<std::vector<int> > graphType;

//int initialization_tests(intVector X, intMatrix NU, TestModelPropensities<doubleVector> PROPENSITIES, intMatrix DG, \
//		SSA_Direct<intVector, intMatrix, TestModelPropensities<doubleVector>, intMatrix > ssa);
int initialization_tests(SSA_Direct<denseVector, denseMatrix, TestModelPropensities<denseVector>, graphType>& ssa);
int onestep_tests(SSA_Direct<denseVector, denseMatrix, TestModelPropensities<denseVector>, graphType>& ssa);

int main()
{
     // initialize initial_population, stoichiometric matrix, propensity functions, and dependancy graph
     denseVector X = TestModelInitialPopulations<denseVector>(); // current population unset
     denseMatrix NU = TestModelStoichiometry<denseMatrix>();
     TestModelPropensities<denseVector> PROPENSITIES;
     graphType DG = TestModelDependencyGraph<graphType>();

     // initialize SSA
     SSA_Direct<denseVector, denseMatrix, TestModelPropensities<denseVector>, graphType > ssa(X,NU,PROPENSITIES,DG);//with dependency graph
     //note that the solver will detect that the problem is small and not use the dependency graph even if it exists - worth test

     // initialization tests
     initialization_tests(ssa);
     std::cout << "..Initialization tests passed." << std::endl;
     
     // one step tests
     onestep_tests(ssa);
     std::cout << "..One step tests passed." << std::endl;

     std::cout << "SSA_Direct.h works fine with DimerDecay.h ." << std::endl;

     return 0;
}

//int initialization_tests(intVector X, intMatrix NU, TestModelPropensities<doubleVector> PROPENSITIES, intMatrix DG, \
//		SSA_Direct<intVector, intMatrix, TestModelPropensities<doubleVector>, intMatrix > ssa)
int initialization_tests(SSA_Direct<denseVector, denseMatrix, TestModelPropensities<denseVector>, graphType >& ssa)
{
    // initial population tests
    // X[0] = 10000, X[1] = X[2] = 0
    assert(ssa.initialPopulation.size() == NumberOfSpecies);
    assert(ssa.initialPopulation[0] == 10000);
    assert(ssa.initialPopulation[1] == 0);
    assert(ssa.initialPopulation[2] == 0);
    std::cout << "....Population initialization test passed." << std::endl;

    // stoichiometric matrix tests
    // nu[0][0]=-1; //reaction 0: S0->0
    // nu[1][0]=-2; nu[1][1]=1; //reaction 1: 2S0->S1
    // nu[2][0]=2; nu[2][1]=-1; //reaction 2: S1->2S0
    // nu[3][1]=-1; nu[3][2]=1; //reaction 3: S1->S2
    assert(ssa.stoichiometry.size() == NumberOfReactions);
    assert(ssa.stoichiometry[0].size() == NumberOfSpecies);
    assert(ssa.stoichiometry[0][0] == -1);
    assert(ssa.stoichiometry[2][0] == 2);
    assert(ssa.stoichiometry[3][2] == 1);
    assert(ssa.stoichiometry[1][2] == 0);
    std::cout << "....Stoichiometric matrix initialization test passed." << std::endl;

    // propensity functions tests
    // reaction 0: S0->0 1.0
    // reaction 1: 2S0->S1 0.002
    // reaction 2: S1->2S0 0.5
    // reaction 3: S1->S2 0.04
    // set current population as {1,1,1}, the propensity functions should return {1.0,0.002,0.5,0.04}
    assert(ssa.propensities.NumberOfReactions == NumberOfReactions);
    ssa.currentPopulation = ssa.initialPopulation;
    for( int i = 0; i < NumberOfSpecies; i++){
	ssa.currentPopulation[i] = 2;
    }
    ssa.calculateAllPropensities();
//    for( int i = 0; i < NumberOfReactions; i++){
//	std::cout << ssa.currentPropensities[i] << std::endl;
//    } 
    denseVector correctPropensities(NumberOfReactions);
    correctPropensities[0] = 2.0;
    correctPropensities[1] = 0.002;
    correctPropensities[2] = 1.0;
    correctPropensities[3] = 0.08;
    for( int i = 0; i < NumberOfReactions; i++){
	assert(ssa.currentPropensities[i] == correctPropensities[i]);
    }
    std::cout << "....Propensities initialization test passed." << std::endl;

    // dependency graph tests
    // dg[0] = {0 1}
    // dg[1] = {0 1 2 3}
    // dg[2] = {0 1 2 3}
    // dg[3] = {2 3}
//    for(int i = 0; i < ssa.dependencyGraph.size(); i++){
//	for(int j = 0; j < ssa.dependencyGraph[i].size(); j++){
//	    std::cout << ssa.dependencyGraph[i][j] << ' ';
//	}
//	std::cout << std::endl;
//    }
    assert(ssa.dependencyGraph.size() == 4);
    assert(ssa.dependencyGraph[0].size() == 2);
    assert(ssa.dependencyGraph[2].size() == 4);
    assert(ssa.dependencyGraph[1][2] == 2);
    assert(ssa.dependencyGraph[3][1] == 3);
    std::cout << "....Dependency graph initialization test passed." << std::endl;

    return 0;
}

int onestep_tests(SSA_Direct<denseVector, denseMatrix, TestModelPropensities<denseVector>, graphType>& ssa)
{
    double startTime = 0.0;
    double endTime = 0.1;
    
    // initialize currentTime, currentPopulation, and Propensities
    // check if they are properly initialized
    ssa.initialize(startTime);
    assert(ssa.currentTime == startTime);
    assert(norm_2(ssa.currentPopulation - ssa.initialPopulation) == 0); // there is no == operator for ublas vector, so use norm_2(difference) == 0
//    std::cout << ssa.propensitySum << std::endl;
    assert(ssa.propensitySum == 109990);

    // run one step
    // check if one and only one reaction has been selected and fired
    ssa.fireReaction(ssa.selectReaction());
    denseVector populationDifference = ssa.currentPopulation - ssa.initialPopulation;
    int flag = 0;
    for( int i = 0; i < ssa.stoichiometry.size(); ++i ){
	if( norm_2(populationDifference - ssa.stoichiometry[i]) == 0 ){
	    ++flag;
	}
    }
    assert(flag == 1);

    return 0;
}
