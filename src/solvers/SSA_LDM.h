#ifndef _SSA_LDM_H_
#define _SSA_LDM_H_

#include <iostream>
#include <vector>
#include <deque>
#include <list>
#include <utility>
#include "Random.h"
#include "StandardDriverTypes.h"
#include "LDMTree.h"

#include "SSA_Base.h"

/**
   @file SSA_LDM.h
   
   @brief Header file for SSA_LDM class
*/

namespace STOCHKIT
{
	/**
     @class SSA_LDM

     @brief Gillespie's Stochastic Simulation Algorithm (SSA): logarithmic direct method.
     
     -This algorithm is generally faster than the direct method only when the number of reactions is
	 around 3000-4000 or more. However, it is usually slower than the constant time algorithm
	 and is therefore not recommended for most problems.

     @param _populationVectorType the population vector type, should be dense
	 @param _stoichiometryType
	 @param _propensitiesFunctorType functor takes reaction index and _populationVectorType population
      and returns the propensity for that reaction
	 @param _dependencyGraphType [expand]
	*/

	template<typename _populationVectorType, 
			 typename _stoichiometryType,
			 typename _propensitiesFunctorType,
			 typename _dependencyGraphType>
 	class SSA_LDM
 	: public SSA_Base<_populationVectorType, _stoichiometryType, _propensitiesFunctorType, _dependencyGraphType>
 	{
 	private:
 		/**
	       Default constructor. This is currently an empty function.
       	*/
      	SSA_LDM();

      	typedef SSA_Base<_populationVectorType, _stoichiometryType, _propensitiesFunctorType, _dependencyGraphType> Base;

	public:
		typedef typename _populationVectorType::value_type populationValueType;//for set single species function


		// Base types
    	typedef typename Base::populationVectorType populationVectorType;
    	typedef typename Base::stoichiometryType stoichiometryType;
    	typedef typename Base::propensitiesType propensitiesType;
    	typedef typename Base::dependencyGraphType dependencyGraphType;

		#ifdef MATRIX_STOICHIOMETRY
		  typedef StandardDriverTypes::stoichiometryRow matrixrow;
		#endif

	protected:
		// Base objects
    	using Base::randomGenerator;
    	using Base::initialPopulation;
    	using Base::stoichiometry;
    	using Base::propensities;
    	using Base::dependencyGraph;
    	using Base::NumberOfSpecies;
    	using Base::NumberOfReactions;
    	using Base::currentTime;
    	using Base::currentPopulation;
    	using Base::currentPropensities;
    	using Base::propensitySum;
    	using Base::previousReactionIndex;
    	using Base::stepsSinceCalculateAllPropensities;
    	using Base::defaultMaxStepsCalculateAllPropensities;
    	using Base::maxStepsCalculateAllPropensities;
    	using Base::detectedVerySmallPropensity;

		/**
			@brief creates the post-ordered dependency graph from the normal dependency graph
	 	*/
		void createDependencyGraphPostOrder(const _dependencyGraphType& depGraph);

		/**
			@brief helper function for LDM
		*/
		void checkTree() {
		  //compare the existing tree to a new tree with current propensities
		  std::vector<double> testPropensities(NumberOfReactions);

		  for (std::size_t i=0; i<testPropensities.size(); ++i) {
		    testPropensities[i]=propensities(i,currentPopulation);
		    std::cout << "testPropensities[" << i << "]= " << testPropensities[i] << " ";
		  }
		  std::cout <<"\n";
		  for (std::size_t i=0; i!=currentPopulation.size(); ++i) {
		    std::cout << "x["<<i<<"]=" << currentPopulation[i] << " ";
		  }
		  std::cout << "\n";

		  LDMTree testTree(NumberOfReactions); 
		   testTree.build(testPropensities);

		  for (std::size_t i=0; i<testPropensities.size(); ++i) {
		    std::cout << "propensityTree[" << i << "].sum=" << propensityTree[i].sum << " vs test=" << testTree[i].sum << std::endl;
		  }  
		}

		/**
			@brief helper class for LDM
		*/
        LDMTree propensityTree;

        /**
			@brief helper class for LDM

			-keep a copy of the propensityTree constructed from initial populations
		*/
        LDMTree initialPropensityTree;

        std::vector<std::vector<std::size_t> > dependencyGraphPostOrder;

	public:
		using Base::seed;
    	using Base::validate;
    	using Base::calculateAllPropensities;

    	using Base::setCurrentTime;
    	using Base::getCurrentTime;

    	/**
	     @brief restart the state of the system to its initial state

	     -changes the current state to whatever the initial state is
	     -recall the initial state is stored in the solver class
	     -this function should be called before every realization

	     @param startTime the time we'd like to start the simulation at
	     @warning in 2.1 times will be handled by driver so this must be edited
      	*/
    	void initialize(double startTime = 0.0);

    	/**
	     @brief usual constructor

	     @param initialPop intial population count of each species in the system
	     @param stoich stoichiometry vector or matrix
	     @param propensitiesFunctor 
	     @param depGraph graph that maintains structure of which reactions will affect
	     the propensities of other reactions
	     @param seed the initial integer seed to our simulation which will determine
	     the sequence of random numbers which are generated
      	*/
		SSA_LDM(const _populationVectorType& initialPop,
			    const _stoichiometryType& stoich,
			    const _propensitiesFunctorType& propensitiesFunctor,
			    const _dependencyGraphType& depGraph,
			    int seed=time(NULL))
		: SSA_Base<_populationVectorType, _stoichiometryType, _propensitiesFunctorType, _dependencyGraphType>(initialPop,stoich,propensitiesFunctor,depGraph, seed),
	          propensityTree(depGraph.size()),
        	  initialPropensityTree(depGraph.size()),
		  dependencyGraphPostOrder(depGraph.size())
		{
			randomGenerator.Seed(seed);
			std::vector<double> initialPropensities(depGraph.size());
			for (std::size_t i=0; i!=depGraph.size(); ++i) {
				initialPropensities[i]=propensities(i,initialPopulation);
			}
			initialPropensityTree.build(initialPropensities);
			createDependencyGraphPostOrder(depGraph);
		}

		/**
			@brief destructor
		*/
		virtual ~SSA_LDM() {
		}

		/**
			@brief time until next rxn of any time will happen based on propensity sum

			-returns infinity if propensitySum is less than or equal to 0
			-issues a warning if propensitySum is less than 0

			@return time until the next reaction
		*/
		double selectStepSize();

		/**
	     @brief fire a reaction

		 -tested other implementations that were cleaner, but this one is faster
	     -updates currentPopulation
	     -updates currentPropensities for all affected reactions (determined by dependencyGraph[reactionIndex])
	     -updates currentTime
      	*/
		bool fireReaction();

		/**
			@brief run an ensemble simulation with output recorded at fixed time intervals
				   
			needs work: number of intervals should be a parameter, 
			calls validate before ensemble
			calls initialize before each realization
			
			@param realizations number of simulations in the ensemble
			@param startTime the initial value of currentTime for each realization
			@param endTime the end time of each realization
			@param Output the class that handles storing the output for the simulation
		*/
        template<typename IntervalOutputType>
	  void simulate(std::size_t realizations, double startTime, double endTime, IntervalOutputType& output, bool doValidate=true);
	
	  /**
       @brief set the current population of species[speciesIndex]

       -this is overridden specifically for the NRM
      */
      virtual bool setSingleSpeciesCurrentPopulation(std::size_t speciesIndex, populationValueType newPopulation);

      /**
       @brief get the current population count of each species of the system as a vector

       @param value_type true: value, parameter will be constant afterwards; false: expression, parameter still depends on other parameters
       */
      virtual bool setParameterValue(std::size_t parameterIndex, std::string newParameterValue, bool value_type);

 };//end SSA_LDM class
}

#define _SSA_LDM_IPP_
#include "SSA_LDM.ipp"
#undef _SSA_LDM_IPP_

#endif
