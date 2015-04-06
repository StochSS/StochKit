#ifndef _SSA_DIRECT_H_
#define _SSA_DIRECT_H_

#include <iostream>
#include <vector>
#include <limits>
#include "Random.h"
#include "StandardDriverTypes.h"

#include "SSA_Base.h"

/**
   @file SSA_Base.h
   
   @brief Header file for SSA_Direct class
*/

namespace STOCHKIT
{
  /**
     @class SSA_Direct

     @brief Gillespie's Stochastic Simulation Algorithm (SSA): direct method.
     
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
    class SSA_Direct
    : public SSA_Base<_populationVectorType, _stoichiometryType, _propensitiesFunctorType, _dependencyGraphType>
    {
    private:
      /**
	 Default constructor. This is currently an empty function.
       */
      SSA_Direct();
      
      typedef SSA_Base<_populationVectorType, _stoichiometryType, _propensitiesFunctorType, _dependencyGraphType> Base;
      
    public:
      // Base types
      typedef typename Base::populationVectorType populationVectorType;
      typedef typename Base::stoichiometryType stoichiometryType;
      typedef typename Base::propensitiesType propensitiesType;
      typedef typename Base::dependencyGraphType dependencyGraphType;
//      using Base::populationVectorType;
//      using Base::stoichiometryType;
//      using Base::propensitiesType;
//      using Base::dependencyGraphType;
      
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
	@brief selects the index of the next reaction to fire based on currentPropensities
	
	-this is a helper function for fireReaction()
	-calls calculateAllPropensities if stepsSinceCalculateAllPropensities is greater than maxStepsCalculateAllPropensities
	
	@see fireReaction()

	@return -1 if there is an error
	
      */
      int selectReaction();

    public:
      // Base methods
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

      	 @param initialPop intial population count of each species in the system
      	 @param stoich stoichiometry vector or matrix
      	 @param propensitiesFunctor 
      	 @param depGraph graph that maintains structure of which reactions will affect
      	 the propensities of other reactions
      	 @param seed the initial integer seed to our simulation which will determine
      	 the sequence of random numbers which are generated
       */
      SSA_Direct(const _populationVectorType& initialPop,
		 const _stoichiometryType& stoich,
		 const _propensitiesFunctorType& propensitiesFunctor,
	       const _dependencyGraphType& depGraph,
	       int seed=time(NULL)):
      SSA_Base<_populationVectorType, _stoichiometryType, _propensitiesFunctorType, _dependencyGraphType>(initialPop,stoich,propensitiesFunctor,depGraph, seed)
	{	
	}
      
      //! compiler-generated copy constructor OK
      //! compiler-generated assignment operator OK
      
      //! destructor
      virtual ~SSA_Direct() {
      }
      
      /**
      	@brief run an ensemble simulation with output recorded at fixed time intervals
      	
      	-THIS FUNCTION IS ONLY IN STOCHKIT 2.1 AS OF NOW FOR DEVELOPMENT PURPOSES AND WILL SOON BE REMOVED
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
      	@brief selects the time step based on the propensitySum
      	
      	-one the time step is selected one should usually fire the
      	reaction directory after

      	@return infinity if propensitySum is less than or equal to 0 or
      	issues a warning if propensitySum is less than 0
      */
      double selectStepSize();
      
      /**
      	@brief fire a reaction
      	
      	-updates currentPopulation
      	-updates currentPropensities for all affected reactions (determined by dependencyGraph[reactionIndex])
      	-updates propensitiesSum
      	-increments stepsSinceCalculateAllPropensities
      	
      	@return 1 if reaction was fired successfully or 0 otherwise
      */
      bool fireReaction();

    };//end SSA_Direct class
}

#define _SSA_DIRECT_IPP_	
#include "SSA_Direct.ipp"
#undef _SSA_DIRECT_IPP_

#endif
