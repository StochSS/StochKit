/*!
\brief the next reaction method of the Stochastic Simulation Algorithm (SSA)
Note that this class behaves differently from the direct method
The fireReaction method also updates currentTime!
To see the time of the next reaction without firing or updating the time, use getNextReactionTime()
*/

#ifndef _SSA_NRM_H_
#define _SSA_NRM_H_

#include <iostream>
#include <vector>
#include <limits>
#include "Random.h"
#include "StandardDriverTypes.h"
#include "BinHeap.h"

#include "SSA_Base.h"

/**
   @file SSA_NRM.h
   
   @brief Header file for SSA_NRM class
*/

namespace STOCHKIT
{
   /**
      @class SSA_NRM
   
      @brief Gibson and Bruck's Next Reaction Method for stochastic simulation.
      
      -This class uses the Next Reaction Method for stochastic simulation.
      -This method draws an exponential random variable for each reaction channel then
      another exponential after each firing of a reaction. This is in contrast with
      the Gillespie (direct) method which requires two random draws per reaction firing.
      
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
    class SSA_NRM
      : public SSA_Base<_populationVectorType, _stoichiometryType, _propensitiesFunctorType, _dependencyGraphType>
    {	
    private:
      /**
	     Default constructor. This is currently an empty function.
      */
      SSA_NRM();
      
      typedef SSA_Base<_populationVectorType, _stoichiometryType, _propensitiesFunctorType, _dependencyGraphType> Base;
   
    public:

      typedef typename _populationVectorType::value_type populationValueType;//for set single species function

      // Base types
      typedef typename Base::populationVectorType populationVectorType;
      typedef typename Base::stoichiometryType stoichiometryType;
      typedef typename Base::propensitiesType propensitiesType;
      typedef typename Base::dependencyGraphType dependencyGraphType;
      
#ifdef MATRIX_STOICHIOMETRY
      using Base::matrixrow;
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
      using Base::previousReactionIndex;
      using Base::stepsSinceCalculateAllPropensities;
      using Base::defaultMaxStepsCalculateAllPropensities;
      using Base::maxStepsCalculateAllPropensities;
      using Base::detectedVerySmallPropensity;

      /**
	     @brief custom created object used to keep track of minimum reaction times

	     @see BinHeap
      */	
      BinHeap rxnHeap;

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
      void initialize(double startTime=0.0);

      /**
	     @brief usual constructor

	     -the constructor must modify the dependency graph slightly to account for a bug
	     where propensities would not change and cause the NRM to crash

	     @param initialPop intial population count of each species in the system
	     @param stoich stoichiometry vector or matrix
	     @param propensitiesFunctor 
	     @param depGraph graph that maintains structure of which reactions will affect
	     the propensities of other reactions
	     @param seed the initial integer seed to our simulation which will determine
	     the sequence of random numbers which are generated

	     @warning Why is constructor in header file?
      */
      SSA_NRM(const _populationVectorType& initialPop,
	      const _stoichiometryType& stoich,
	      const _propensitiesFunctorType& propensitiesFunctor,
	      const _dependencyGraphType& depGraph,
	      int seed=time(NULL)) :
      SSA_Base<_populationVectorType, _stoichiometryType, _propensitiesFunctorType, _dependencyGraphType>(initialPop,stoich,propensitiesFunctor,depGraph, seed),
	rxnHeap((int)NumberOfReactions)
      	{
      	  randomGenerator.Seed(seed);

      	  //fix dependency graph to ensure reaction is in its own dependency graph (previous code would
      	  //not include in A->A+B because propensity is not changed but NRM needs it)
      	  for (std::size_t i=0; i!=dependencyGraph.size(); ++i) {
      	    //insert i into the dependency graph for rxn i
      	    dependencyGraph[i].push_back(i);
      	    //sort and remove duplicates
      	    std::sort(dependencyGraph[i].begin(),dependencyGraph[i].end());
      	    dependencyGraph[i].erase(std::unique(dependencyGraph[i].begin(),dependencyGraph[i].end()),dependencyGraph[i].end());
      	  }
      	}

      //! compiler-generated copy constructor OK

      //! compiler-generated assignment operator OK

      //! destructor
      virtual ~SSA_NRM() {}

      /**
	     @brief run an ensemble simulation with output recorded at fixed time intervals
	
	     -output must have a conforming initialize, getOutputTimes, and record method
	     -outputTimes should be set in output prior to calling simulate
	     -if doValidate=true (the default) calls validate before ensemble
	     -calls initialize before each realization

	     @param realizations number of simulations in the ensemble
	     @param startTime the initial value of currentTime for each realization
	     @param endTime the end time of each realization
	     @param Output the class that handles storing the output for the simulation
      */
      template<typename IntervalOutputType>
	    void simulate(std::size_t realizations, double startTime, double endTime, IntervalOutputType& output, bool doValidate=true);

      /**
	     @brief Queries heap to find when next rxn of ANY type will happen

	     @return the time of the next reaction
      */
      double selectStepSize();

      /**
	     @brief fire a reaction

	     -updates currentPopulation
	     -updates currentPropensities for all affected reactions (determined by dependencyGraph[reactionIndex])
	     -updates currentTime
      */
      bool fireReaction();

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

    };//end SSA_NRM class
}

#define _SSA_NRM_IPP_	
#include "SSA_NRM.ipp"
#undef _SSA_NRM_IPP_

#endif
