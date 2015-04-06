/*!
\brief the Base class of the Stochastic Simulation Algorithm (SSA)
*/

#ifndef _SSA_BASE_H_
#define _SSA_BASE_H_

#include <iostream>
#include <sstream>
#include <vector>
#include <limits>
#include "Random.h"
#include "StandardDriverTypes.h"

/**
   @file SSA_Base.h
   
   @brief Header file for SSA_Base class which serves as a parent class to all SSA's
   by containing all members and functions which the different SSA's have in common
 */

namespace STOCHKIT
{
  /**
      @class SSA_Base

      @brief Serves as a parent class to all SSA's by containing all members and functions which
      the different SSA's have in common

      @param _populationVectorType the population vector type, should be dense
      @param _stoichiometryType the stoichiometry matrix type, could be vector<vector> or boost::matrx
      @param _propensitiesFunctorType functor takes reaction index and _populationVectorType population
      and returns the propensity for that reaction
      @param _dependencyGraphType [expand]
  */
  template<typename _populationVectorType, 
    typename _stoichiometryType,
    typename _propensitiesFunctorType,
    typename _dependencyGraphType>
    class SSA_Base
    {
    public:	
      typedef _populationVectorType populationVectorType;
      typedef typename _populationVectorType::value_type populationValueType;
      typedef _stoichiometryType stoichiometryType;
      typedef _propensitiesFunctorType propensitiesType;
      typedef boost::numeric::ublas::vector<double> propensitiesVectorType;
      typedef _dependencyGraphType dependencyGraphType;
      
#ifdef MATRIX_STOICHIOMETRY
      typedef StandardDriverTypes::stoichiometryRow matrixrow;
#endif
      
    protected:
      /**
	  @brief Class member to generate all random numbers
	  
	  -Note that STOCHKIT::RandomGenerator is a wrapper class to a boost functions and classes
	  -NOTE:change the RandomGenerator class to swap generators
       */
      STOCHKIT::RandomGenerator randomGenerator;
      
      /**
	  @brief vector describing the initial populations of the various species

	  -currentPopulation should be set to initialPopulation at the beginning of eachrealization in initialize()
      */
      _populationVectorType initialPopulation;
      
      /**
	 @brief the stoichiometric matrix which tells us how each reaction will affect our system
    
	 -should actually be a dense vector of (usually sparse) vectors of dimension NumberOfReactions x NumberOfSpecies
	 so that currentPopulation+=stoichiometry[reactionNumber] modifies the population based on
	 the stoichiometry of the given reaction number
      */
      _stoichiometryType stoichiometry;
      
      /**
	@brief the propensities functor which stores the propensity of each each reaction
	
	propensities(rxn, pop) returns the propensity of reaction number rxn based on population pop
	after a simulation step, currentPropensities[rxn] should equal propensities(rxn, currentPopulation)
	but since propensities() is a function call, accessing current propensities should be done with
	currentPropensities for efficiencies sake
	
	@see currentPropensities
      */
      _propensitiesFunctorType propensities;

      /**
	 @brief the dependency graph which describes the propensities that are affected by each reaction

	should be a dense vector (size=NumberOfReactions) of variable length vectors
	where dependencyGraph[rxn] returns the variable length vector of reaction indices that are affected by reaction rxn
	e.g. if dependencyGraph[4][0]=2 and dependencyGraph[4][1]=4 and dependencyGraph[4][2]=5 then dependencyGraph[4].size()=3 and
	reaction 4 affects reactions 2, 4, and 5
	
	@see fireReaction
      */
      _dependencyGraphType dependencyGraph;
      
      
      /**
	 @brief total number of species in the system
      */
      std::size_t NumberOfSpecies;
      
      /**
	 @brief number of reactions in the system
      */
      std::size_t NumberOfReactions;
      
      /**
	 @brief current time of the simulation, should be incremented at each simulation time step
      */
      double currentTime;
      
      /**
	 @brief current population of the simulation
      */
      _populationVectorType currentPopulation;
      
      /**
	 @brief current propensities of the simulation

	 @warning type should probably be changed to std::vector<double>
       */
      propensitiesVectorType currentPropensities;

      /**
	 @brief current sum of propensities of the simulation, used to determine time step
      */ 
      double propensitySum;
      
      /**
	 @brief index of the last reaction that fired
	 
	 default and error value is -1
      */
      int previousReactionIndex;
      
      /**
      @brief counter for simulation steps taken since the last time calculateAllPropensities was called

      is used to ensure that roundoff errors in propensities do not accumulate
      selectReaction() uses a simple strategy to recalculate all propensities using maxStepsCalculateAllPropensities
      for conservative strategy, see S. Mauch, M. Stalzer "Efficient formulations for
      exact stochastic simulation of chemical systems" IEEE/ACM Trans. on Comp. Bio. and Bioinformatics, 30 April 2009
      
      @see calculateAllPropensities
      @see defaultMaxStepsCalculateAllPropensities
      @see maxStepsCalculateAllPropensities
      @see selectReaction()
      */
      std::size_t stepsSinceCalculateAllPropensities;
      
      /**
	 @brief default value for maxStepsCalculateAllPropensities
      */
      static const std::size_t defaultMaxStepsCalculateAllPropensities=10000;
      
      /**
	 @brief maximum number of steps allowed before calling calculateAllPropensities
	 
	 @see stepsSinceCalculateAllPropensities
      */
      std::size_t maxStepsCalculateAllPropensities;
      
    private:
      /**
	 @brief default constructor
	 
	 -not implemented
      */
      SSA_Base();
      
    public:

      /**
	 @brief usual constructor
	 
	 compiler-generated copy constructor OK
	 compiler-generated assignment operator OK
      */
      SSA_Base(const _populationVectorType& initialPop,
	       const _stoichiometryType& stoich,
	       const _propensitiesFunctorType& propensitiesFunctor,
	       const _dependencyGraphType& depGraph,
	       int seed=time(NULL));
      
      /**
	 @brief destructor
      */
      virtual ~SSA_Base() {
      }
      
      /**
	 @brief seeds the random number generator
      */
      void seed(int seed);
      
      /**
	 @brief empty prepare() function for all the solvers that do not need to prepare
      */
      virtual bool prepare(std::size_t realizations, double startTime, double endTime);
      
      /**
	 @brief consistency checks to validate that the class is set up properly for a simulation, should be
	 called before an ensemble as in simulate(...,doValidate=true)
      */
      bool validate(double startTime, double endTime);
      
      /**
	@brief update all the propensities

	-updates currentPropensities by calling the propensities functor for each reaction using currentPopulation
	-updates propensitySum
	-resets stepsSinceCalculateAllPropensities to 0
      */
      void calculateAllPropensities();
      
      /**
	 @brief returns current time of the system

	 @warning we eventually plan to have the driver keep track of time so this should be removed
       */
      double getCurrentTime();
      
      /**
	 @brief set current time of the system

	 @warning we eventually plan to have the driver keep track of time so this should be removed
       */
      bool setCurrentTime(double newCurrentTime);
      
      /**
	 @brief get the current population count of each species of the system as a vector
       */
      _populationVectorType getCurrentPopulation();

      /**
   @brief get the Initial population count of each species of the system as a vector
       */
      _populationVectorType getInitialPopulation();
      
      /**
	 @brief set the current population count of every species of the system as a vector
       */
      bool setCurrentPopulation(_populationVectorType& newPopulation);

      /**
	 @brief set the current population of species[speciesIndex]

   -this is virtual because certain solvers such as the NRM solver need to override this function
       */
      virtual bool setSingleSpeciesCurrentPopulation(std::size_t speciesIndex, populationValueType newPopulation);

      /**
	 @brief get the current population count of each species of the system as a vector
       */
      propensitiesVectorType getCurrentPropensities();

      /**
	 @brief get the current population count of each species of the system as a vector

	 @param value_type true: value, parameter will be constant afterwards; false: expression, parameter still depends on other parameters
       */
      virtual bool setParameterValue(std::size_t parameterIndex, std::string newParameterValue, bool value_type);

      /**
	 @brief get previous reaction index for ssa, currently undefined behavior for tau-leaping
       */
      int getPreviousReactionIndex();
      
      /**
	 @brief checks and warns the user about very small propensities which can cause round off error
       */
      bool detectedVerySmallPropensity;

      ParameterSet& referenceToParametersList();

      /**
        @brief resets parameters to initial values

        -used particularly when events change a parameter value and the value must be reset in the
        next realization
      */
      void resetParameters();

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE

#ifndef BOOST_RANDOM_NO_STREAM_OPERATORS
      /**
	 @brief get current Random Number Generator state
       */
      std::string getCurrentRNGState();

      /**
	 @brief set current Random Number Generator state
       */
      bool setCurrentRNGState(std::string RNGStateString);
#endif

#endif
 
    };//end SSA_Base class
}

#define _SSA_BASE_IPP_	
#include "SSA_Base.ipp"
#undef _SSA_BASE_IPP_

#endif
