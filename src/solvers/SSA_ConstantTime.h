/*!
\brief the amortized constant-complexity (composition-rejection) method of the Stochastic Simulation Algorithm (SSA)
*/

#ifndef _SSA_CONSTANTTIME_H_
#define _SSA_CONSTANTTIME_H_

#include <iostream>
#include <vector>
#include <deque>
#include <list>
#include <utility>
#include <math.h>
#include "Random.h"
#include "StandardDriverTypes.h"
#include "ConstantTimeGroup.h"
#include "ConstantTimeGroupCollection.h"

#include "SSA_Base.h"

/**
  @file SSA_ConstantTime.h

  @brief Header file for the constant complexity (in the number of reaction channels) method
*/

namespace STOCHKIT
{
  /**
     @class SSA_ConstantTime

     @brief This is an exact SSA that scales constantly in the number of reaction channels.

     -This algorithm is generally faster than the direct method only when the number of reactions is
     around 3000-4000 or more.  For a description of the algorithm see
     A. Slepoy, A.P. Thompson, and S.J. Plimpton. J Chem Phys 128(20):205101 2008. or
     S. Mauch, M. Stalzer "Efficient formulations for exact stochastic simulation of chemical systems"
     IEEE/ACM Trans. on Comp. Bio. and Bioinformatics, 30 April 2009.
     
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
    class SSA_ConstantTime
    : public SSA_Base<_populationVectorType, _stoichiometryType, _propensitiesFunctorType, _dependencyGraphType>
    {	
    private:
      /**
	 Default constructor. This is currently an empty function.
      */
      SSA_ConstantTime();
      
      typedef SSA_Base<_populationVectorType, _stoichiometryType, _propensitiesFunctorType, _dependencyGraphType> Base;
      
    public:
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

      /**
	     @brief groups of propensities
      */
      ConstantTimeGroupCollection groups;
      
      /**
    	 @brief propensities based on initial population

    	 -this object is unique to the constantTime SSA and is required because of a
    	 implementation detail related to the "groups" in this solver

    	 @warning type should probably be changed to std::vector<double>
      */
      boost::numeric::ublas::vector<double> initialPropensities;
      
      /**
	     @brief groups based on initial conditions
      */
      ConstantTimeGroupCollection initialGroups;

      /**
      	@brief update the "groups" data structure
      	
      	@param affectedReactionIndex the index of the propensity that changed
      */
      void updateGroups(int affectedReactionIndex, double oldPropensity);

      /**
        @brief selects the index of the next reaction to fire based on currentPropensities

        @return -1 if there is an error
  
      */
      int selectReaction();
      
    public:

      typedef typename _populationVectorType::value_type populationValueType;//for set single species function


      //Base methods
      using Base::seed;
      using Base::validate;

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

    	 @warning why is the definition of the constructor in the header file?
       */
      SSA_ConstantTime(const _populationVectorType& initialPop,
		       const _stoichiometryType& stoich,
		       const _propensitiesFunctorType& propensitiesFunctor,
		       const _dependencyGraphType& depGraph,
		       int seed=time(NULL)) :
      SSA_Base<_populationVectorType, _stoichiometryType, _propensitiesFunctorType, _dependencyGraphType>(initialPop,stoich,propensitiesFunctor,depGraph, seed),
            	 groups(NumberOfReactions),
            	 initialPropensities(NumberOfReactions),
            	 initialGroups(NumberOfReactions)
    	{
    	  randomGenerator.Seed(seed);
    	  for (std::size_t i=0; i!=NumberOfReactions; ++i){
    	    initialPropensities[i]=propensities(i,initialPopulation);
    	  }
    	  groups.build(initialPropensities);
    	  initialGroups=groups;
    	}
      
      //! destructor
      virtual ~SSA_ConstantTime() {
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

      	@bug simulate segfaults when trying to fireReaction sometimes
      */
      template<typename IntervalOutputType>
	void simulate(std::size_t realizations, double startTime, double endTime, IntervalOutputType& output, bool doValidate=true);
      
      /**
       @brief selects the step size based on the propensitySum
       
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
      

    };//end SSA_ConstantTime class
}

#define _SSA_CONSTANTTIME_IPP_
#include "SSA_ConstantTime.ipp"
#undef _SSA_CONSTANTTIME_IPP_

#endif
