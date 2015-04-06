#ifndef _SSA_ODM_H_
#define _SSA_ODM_H_

#include<iostream>
#include <vector>
#include <limits>
#include <time.h>
#ifndef WIN32
#include <sys/time.h>
#endif
#include "SSA_Direct.h"
#include "Random.h"
#include "StandardDriverTypes.h"
#include "CustomPropensitySet.h"

/**
   @file SSA_ODM.h
   
   @brief header file for SSA_ODM. Note this solver does NOT have an accompanying .ipp file
*/
namespace STOCHKIT
{
  /**
     @class SSA_ODM

     @brief Cao, Li and Petzold's Optimized Direct Method of the Stochastic Simulation Algorithm (SSA).

     -runs a single trajectory to determine the most commonly fired reactions and
     then orders these reactions in descending order so that when selecting a reaction
     the search only must go through the beginning of the list of reactions usually

     @param _populationVectorType the population vector type, should be compatible with
     _stoichiometryType (see below), and as input to _propensitiesFunctorType, must have a .size() function
     @param _stoichiometryType should be compatible with _populationVectorType--that is,
     when a reaction fires, we will take _populationVectorType+=_stoichiometryType[reactionIndex].
     size should be equal to number of reactions, must have a .size() method
     @param _propensitiesFunctorType functor takes reaction index and _populationVectorType population
     and returns the propensity for that reaction
     @param _dependencyGraphType [expand]
     
   */

  template<typename _populationVectorType, 
    typename _stoichiometryType,
    typename _propensitiesFunctorType,
    typename _dependencyGraphType>
    class SSA_ODM
    : public SSA_Direct <_populationVectorType, _stoichiometryType, _propensitiesFunctorType, _dependencyGraphType>
    {	
    private:
      /**
	 Default constructor. This is currently an empty function.
      */
      SSA_ODM();
    
      typedef SSA_Direct<_populationVectorType, _stoichiometryType, _propensitiesFunctorType, _dependencyGraphType> Base;

    protected:
      using Base::NumberOfReactions;
      using Base::propensities;
      using Base::currentTime;

      /**
	 @class FireFrequencyReport
	 
	 @brief to store the index and the frequency of the reactions in the first trial run
      */
      class FireFrenquencyRecord{
      public:
	int reactionIndex;//the index of the reaction
	int fireFrequency;//firing frequency of the reaction with index reactionIndex
      };
      
      /**
	 @class ValueGreater
	 
	 @brief sets the criterion for comparing two firing frequency
	 
	 -takes two values of the same type than those contained in the range, 
	 returns true if the first argument goes before the second argument, 
	 and false otherwise. 
      */
      class ValueGreater :
      public std::binary_function<FireFrenquencyRecord, FireFrenquencyRecord, bool> {
      public:
	bool operator()(const FireFrenquencyRecord& x, const FireFrenquencyRecord& y) const {
	  return x.fireFrequency > y.fireFrequency;
	}
      };
      
      /**
	 the vector to record the reactions' firing frequency
      */
      std::vector<FireFrenquencyRecord> FrequencyVector;
      std::vector<std::size_t> newReactionOrder;
      
      /**
	 the vector to record the old index for current index
      */
      std::vector<int> IndexMap;
    
      /**
	 the propensities fucntion type
      */
      typedef double (_propensitiesFunctorType::* PropensityMember) (_populationVectorType&);
     
      /**
	 @brief Set the initial value for the firing frequency of the reactions, initial are all 0s
      */
      void InitializeFrequencyVector()
      {
	FrequencyVector.resize(NumberOfReactions);	
	newReactionOrder.resize(NumberOfReactions);
	IndexMap.resize(NumberOfReactions);	
	//printf("number of Reactions %d \n", numberofReactions);
	for(std::size_t i=0; i < NumberOfReactions; i++){
	  FrequencyVector[i].reactionIndex = i;
	  FrequencyVector[i].fireFrequency = 0;
	  newReactionOrder[i] = i;
	  IndexMap[i] = i;
	  
#ifdef DEBUG
	  printf("FrequencyVector[%d].reactionIndex = %d\n", i, FrequencyVector[i].reactionIndex);
	  printf("FrequencyVector[%d].fireFrequency = %d\n", i, FrequencyVector[i].fireFrequency);
#endif
	}
      }
      
      /**
	 @brief firing frequncy increases 1 for fired reaction
      */
      void RecordFireReaction(int previousIndex)
      {
	if (previousIndex!=-1) {
	  FrequencyVector[previousIndex].fireFrequency++;
	}
      }
      
      /**
	@brief Get the reaction fring requency
      
	run the whole simulation once, recording the firing frequency for each reaction
      */
      void getReactionsFireFrequency(double startTime, double endTime)
      {
	currentTime = startTime;
	Base::initialize(startTime);
	currentTime+=Base::selectStepSize();
	while (currentTime<endTime) {
	  Base::fireReaction();
	  RecordFireReaction(Base::previousReactionIndex);
	  currentTime+=Base::selectStepSize();
	}
      }
    
      /**
	 @brief reorder the vectors
	 
	 reordering the vectos according to the firing frequency saved in the Frequency vector
	 Please remember the _reorderVectorType has to be a vector
      */
      template<typename _reorderVectorType>
	_reorderVectorType ReorderVectors(_reorderVectorType& orgVector) {
	int numberofReactions = NumberOfReactions;
	_reorderVectorType rV(numberofReactions);
	for(std::size_t i=0; i <  NumberOfReactions; i++){
	  int rI = FrequencyVector[i].reactionIndex;
	  rV[i] = orgVector[rI];
	}
	return rV;
      }
      
      /**
	 @brief Generate new DG
      */
      _dependencyGraphType GenerateDG(_dependencyGraphType oriDg)
      {
	//obtain the new index based on old index
	for(std::size_t i=0; i <  NumberOfReactions; i++){
	  int rI = FrequencyVector[i].reactionIndex;
	  IndexMap[rI] = i;
	}
	_dependencyGraphType dg(NumberOfReactions);
	for(std::size_t i=0; i <  NumberOfReactions; i++){
	  int rI = FrequencyVector[i].reactionIndex;
	  int dgSize = oriDg[rI].size();
	  dg[i].resize(dgSize);
	  for(int j=0; j < dgSize; j++){
	    int tmpIndex = oriDg[rI][j];
	    int newI = IndexMap[tmpIndex];
	    dg[i][j] = newI;
	  }
	}
	return dg;
      }
      
      /**
	 @brief Reorder the reactions
      */
      void reorderReactions()
      {
	sort(FrequencyVector.begin(), FrequencyVector.end(), ValueGreater());
#ifdef DEBUG
	printf("New order:\n");
	for(int i=0; i < NumberOfReactions; i++){
	  printf("FrequencyVector[%d].reactionIndex = %d\n", i, FrequencyVector[i].reactionIndex);
	  printf("FrequencyVector[%d].fireFrequency = %d\n", i, FrequencyVector[i].fireFrequency);
	}
#endif
        for(std::size_t i=0; i <  NumberOfReactions; i++){
          newReactionOrder[i] = FrequencyVector[i].reactionIndex;
        }

	Base::stoichiometry = ReorderVectors<_stoichiometryType>(Base::stoichiometry);
	//Base::propensities._propensityFunctions = ReorderVectors<vector<PropensityMember> >(Base::propensities._propensityFunctions);
	//Base::propensities.propensities = ReorderVectors<_propensitiesFunctorType::tempType >(Base::propensities.propensities);
	//typedef typename CustomPropensitySet<_populationVectorType>::tempType myType;
	//Base::propensities.propensities = ReorderVectors<myType>(Base::propensities.propensities);
	Base::propensities.reOrderPropensities(newReactionOrder);
	Base::dependencyGraph = GenerateDG(Base::dependencyGraph);
	
      }
    
      /**
	 @brief Presimulating the model once with start time and end time, 
	 collecting the firing requencies of reactions and
	 reorder the stoichiometry matrix and penpesities
	 
	 @param startTime the initial value of currentTime for each realization
	 @param endTime the end time of each realization
      */
      void presimulation(double startTime, double endTime)
      {
	InitializeFrequencyVector();
	getReactionsFireFrequency(startTime, endTime);
	reorderReactions();
      }


    public:	
      /**
	 @brief Prepare all the simulation by calling presimulation() 
      */
      virtual bool prepare(std::size_t /*realizations*/, double startTime, double endTime)
      {     
#ifdef DEBUG
          printf("ori stoichiometry:\n");
          for(int i=0; i < NumberOfReactions; i++){
            for(int j=0; j< Base::NumberOfSpecies; j++){
              printf("%f ", Base::stoichiometry[i][j]);
            }
            printf("\n");
          }
          printf("Start presimulation\n");
#endif
          presimulation(startTime, endTime);
#ifdef DEBUG
          printf("Finish presimulation\n");
          printf("\n new stoichiometry:\n");
          for(int i=0; i < NumberOfReactions; i++){
            for(int j=0; j< Base::NumberOfSpecies; j++){
              printf("%f ", Base::stoichiometry[i][j]);
            }
            printf("\n");
          }
#endif
	 return true;
       }
 
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
    SSA_ODM(const _populationVectorType& initialPop,
	    const _stoichiometryType& stoich,
	    const _propensitiesFunctorType& propensitiesFunctor,
	    const _dependencyGraphType& depGraph, 
	    int seed=time(NULL)):
      SSA_Direct<_populationVectorType, _stoichiometryType, _propensitiesFunctorType, _dependencyGraphType>(initialPop,stoich,propensitiesFunctor,depGraph, seed),
	FrequencyVector(Base::NumberOfReactions)
	  {
	  }
      
      //! destructor
      virtual ~SSA_ODM() {
      }
      
      /**
	 @brief run an ensemble simulation with output recorded at fixed time intervals
	 
	 @warning needs work:
	 the reactions need reordered or not should be a parameter,
	 number of the intervals should be a parameter,
	 Output should be specifically a base "IntervalOutputObject" (or similar) class with a default
	 or Output should be a template parameter
	 need to specify the functionality that Output must provide
	 there is an error in the main loop--it always records output one reaction after it should
	 calls validate before ensemble
	 calls initialize before each realization
	 
	 @param realizations number of simulations in the ensemble
	 @param startTime the initial value of currentTime for each realization
	 @param endTime the end time of each realization
	 @param Output the class that handles storing the output for the simulation
      */
      template<typename IntervalOutputType>
	void simulate(std::size_t realizations, double startTime, double endTime, IntervalOutputType& output, bool doValidate=true)
	{
	  
#ifdef DEBUG
	  printf("ori stoichiometry:\n");
	  for(int i=0; i < NumberOfReactions; i++){
	    for(int j=0; j< Base::NumberOfSpecies; j++){
	      printf("%f ", Base::stoichiometry[i][j]);
	    }
	    printf("\n");
	  }
	  printf("Start presimulation\n");
#endif
	  presimulation(startTime, endTime);
#ifdef DEBUG
	  printf("Finish presimulation\n");
	  printf("\n new stoichiometry:\n");
	  for(int i=0; i < NumberOfReactions; i++){
	    for(int j=0; j< Base::NumberOfSpecies; j++){
	      printf("%f ", Base::stoichiometry[i][j]);
	    }
	    printf("\n");
	  }
#endif
	  Base::template simulate<IntervalOutputType>(realizations, startTime, endTime, output, doValidate);
	}

    };//end class
}

#endif
