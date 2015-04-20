/*
 *  test_tauleaping.cpp
 *  
 */
#include<iostream>
#include<string>
#include <fstream>
#include <vector>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include<stdio.h>
#include "TauLeaping.h"
#include "TestModel.h"
#include "Output.h"

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

	denseVector X=TestModelInitialPopulations<denseVector>();//the population vector should always be dense
  
	//sparseMatrix NU=TestModelStoichiometry<sparseMatrix>();//sparseMatrix is slower for this small problem
	denseMatrix NU=TestModelStoichiometry<denseMatrix>(); 
	TestModelPropensities<denseVector> PROPENSITIES;  
	graphType DG=TestModelDependencyGraph<graphType>();

//  TauLeaping<denseVector, sparseMatrix, TestModelPropensities<denseVector> > tauleaping(X,NU,PROPENSITIES);
  	TauLeaping<denseVector, denseMatrix, TestModelPropensities<denseVector> > tauleaping(X,NU,PROPENSITIES);
//	note that the solver will detect that the problem is small and not use the dependency graph even if it exists

//	do an ensemble (1000 dimer decay to t=10 runs in about 10s with dense/dense/noDG on Kevin's MacBook, and >1min with sparse/sparse/DG)
//	in a real driver, these variables would be read in as parameters

	int numRuns=10;
	double t0=0;
	double tf=1;
	int intervals=1;
	int speciesChoice=0;  // 0: All species, 1: Species from File
	int timeChoice=2;  // 0:Intervals, 1: Times from file, Other: Endtimes only
	int statChoice=2; //  0:Both data and stat , 1: Stat only, Other: Data only
	int maxNumSpecies=X.size();
	double epsilon=0.02;
	vector<int> selectedSpecies;
	vector<double> selectedTimes;
	string outputFile="dimertau.txt";
	string speciesFile="chosenSpecies.txt";
	string timesFile="times.txt";


	if ( speciesChoice == 1)
		chooseSpecies(selectedSpecies, maxNumSpecies, speciesFile);
	else
	{
		for (int i =0; i < maxNumSpecies; i++)
			selectedSpecies.push_back(i);
	}


	if (timeChoice == 0)
		chooseTimes(selectedTimes, t0, tf, intervals);
	else
	{ 
		if (timeChoice == 1)
			chooseTimes(selectedTimes, timesFile, t0, tf);
		else
			chooseTimes(selectedTimes, tf);
	}


	Output theOutput(selectedSpecies, selectedTimes,numRuns,statChoice, outputFile);
	tauleaping.SetEpsilon(epsilon);
	tauleaping.simulate(numRuns,t0,tf,theOutput);
	theOutput.flushSpecies();
	return 0;
}


void chooseSpecies(vector<int> &whichSpecies,int maxNumSpecies, string filename)
{
    int x;
    ifstream inFile;
       
	   inFile.open(filename.c_str());
                 if (!inFile) {
                    cout << "Unable to open file";
                    exit(1);    // terminate with error
                  }
    
       while (inFile >> x) {
         if (x >= maxNumSpecies){
	   cout << "WARNING: Species number "<< x <<" does not exist"<<endl;
	 }else{
	   whichSpecies.push_back(x);
         }
        
    }
    
    inFile.close();
    
    
    }

void  chooseTimes(vector<double> &times, double t0, double tf, int num)
{
  double timeStep;

   timeStep=(tf-t0)/num;
   for (int i=0; i<=num; i++){
          times.push_back(t0+i*timeStep);
           }
  
}

void  chooseTimes(vector<double> &times, double tf)
{
     times.push_back(tf);
  
}


void  chooseTimes(vector<double> &times, string filename, double t0, double tf)
{
     double t;
     ifstream inFile;
       
	   inFile.open(filename.c_str());
                 if (!inFile) {
                    cout << "Unable to open file";
                    exit(1);    // terminate with error
                  }
    
       while (inFile >> t) {
         if ((t > tf)   || (t < t0)){ 
	   cout << "WARNING: Time  "<< t <<" out of range"<<endl;
	 }else{
	   times.push_back(t);
         }
        
    }
    
    inFile.close();
   
}


