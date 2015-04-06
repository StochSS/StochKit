/*!
	\brief the direct method of the Stochastic Simulation Algorithm (SSA)
*/

#if !defined(_SSA_NRM_IPP_)
#error This file is the implementation of SSA_NRM
#endif
#include "CommandPassAux.h"

namespace STOCHKIT
{

template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
void
SSA_NRM<_populationVectorType, 
		   _stoichiometryType,
		   _propensitiesFunctorType,
		   _dependencyGraphType>::
initialize(double startTime) {
	previousReactionIndex=-1;
	currentTime=startTime;
	currentPopulation=initialPopulation;
	calculateAllPropensities();

	//draw random time for each reaction and store in intialTimes
	int nRxns = (int)NumberOfReactions;
    double* initialTimes = (double*)malloc( nRxns * sizeof(double) );
    for(int i=0; i<nRxns; i++ )
    	{
            if(currentPropensities[i] == 0)
	      		initialTimes[i] = std::numeric_limits<double>::infinity();
	    	else
	      		initialTimes[i] = randomGenerator.Exponential(1.0/currentPropensities[i]);
	  	}
     
	//put in to BinHeap and order BinHeap
	rxnHeap.initializeHeap( initialTimes );

	free(initialTimes);
}

template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
inline
bool
SSA_NRM<_populationVectorType, 
		   _stoichiometryType,
		   _propensitiesFunctorType,
		   _dependencyGraphType>::
fireReaction() {

	//get the index and the time it took to fire the rxn
	std::size_t firedReactionIndex = rxnHeap.getNextRxnNumber();
	double timeUntilFiredReaction = rxnHeap.getNextRxnTime();

	std::vector<int> tempPop;
	tempPop.resize(currentPopulation.size());
	for(std::size_t l=0; l<currentPopulation.size(); l++)
	{
		tempPop[l] = currentPopulation[l];
	}

	//change populations according to reaction
	currentPopulation+=stoichiometry[firedReactionIndex]; 

	bool isT = false;
	for(std::size_t l=0; l<currentPopulation.size(); l++)
	{
		if( currentPopulation[l] - tempPop[l] > 0 )
		{
			isT = true;
		}
	}

	//subtract time that passed while last rxn occurred
	double oldNextRxnTime;

	std::vector<double> newNextRxnTime;
	newNextRxnTime.resize( NumberOfReactions );

	for (std::size_t i=0; i<NumberOfReactions; ++i)
  	{
  		oldNextRxnTime = rxnHeap.getRxnTime(i);
  		newNextRxnTime[i] = oldNextRxnTime - timeUntilFiredReaction;   
  	}

  	//modulate the next rxn times of the rxns whose propensities were Affected
  	//i.e. apply formula to newNextRxnTimes to account for changed propensities
  	for (std::size_t i=0; i!=dependencyGraph[firedReactionIndex].size(); ++i)
  	{
  		std::size_t dependentRxnIndex = dependencyGraph[firedReactionIndex][i];

  		double oldPropensity, newPropensity;
  			
	    oldPropensity=currentPropensities[dependentRxnIndex];

	    currentPropensities[dependentRxnIndex]=propensities(dependentRxnIndex,currentPopulation);
	    newPropensity = currentPropensities[dependentRxnIndex];

	    //handle reactions with positive propensity normally
	    if( newPropensity != 0.0 )
	    {
		    if( dependentRxnIndex != firedReactionIndex && oldPropensity != 0.0 )
		     	newNextRxnTime[dependentRxnIndex] = (oldPropensity/newPropensity) * newNextRxnTime[dependentRxnIndex];

		    //draw new random number for the actual rxn fired and when adjustment
		    //doesn't work due to an oldPropensity value of zero
		    else
		     	newNextRxnTime[dependentRxnIndex] = randomGenerator.Exponential( 1.0/newPropensity );
	    }

	    //handle reactions with 0 propensities by setting the next rxn time to infinity
	    else
	    {
	      	newNextRxnTime[dependentRxnIndex] = std::numeric_limits<double>::infinity();
	    }
  	}

  	//update rxn time in rxn heap
  	for (std::size_t i=0; i<NumberOfReactions; ++i)
  	{
  		rxnHeap.setNewRxnTime(i, newNextRxnTime[i]);
  	}
    
	previousReactionIndex = firedReactionIndex;

	stepsSinceCalculateAllPropensities++;
	return true;
}

template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
template<typename IntervalOutputType>
void
SSA_NRM<_populationVectorType, 
		   _stoichiometryType,
		   _propensitiesFunctorType,
		   _dependencyGraphType>::
simulate(std::size_t realizations, double startTime, double endTime, IntervalOutputType& output, bool doValidate) {

	if (doValidate) {
		if (!validate(startTime,endTime)) {
			std::cerr << "StochKit ERROR (SSA_NRM::simulate): validate() failed, simulation aborted\n";
			exitFunc(1);
		}		
	}

	if (!output.initialize(realizations,startTime,endTime,initialPopulation)) {
		std::cerr << "StochKit ERROR (SSA_NRM::simulate): initialization of output object failed, simulation aborted\n";
		exitFunc(1);
	}
	
	std::vector<double> outputTimes = output.getOutputTimes();
	std::size_t totalIntervals=outputTimes.size();

	std::size_t currentInterval;

	for (std::size_t currentRealization=0; currentRealization!=realizations; ++currentRealization) {
	  initialize(startTime);

	    currentInterval=0;    
	    //updateCurrentTime();
	    //while (currentTime<endTime) {
		while (rxnHeap.getNextRxnTime()<endTime) {
			//while (currentInterval<totalIntervals && currentTime >=outputTimes[currentInterval]) {
			while (currentInterval<totalIntervals && rxnHeap.getNextRxnTime() >=outputTimes[currentInterval]) {
				output.record(currentRealization,currentInterval,currentPopulation);
				currentInterval++;
			}

			fireReaction();
			//updateCurrentTime();//fireReaction now updates
		}
	    while (currentInterval<totalIntervals && rxnHeap.getNextRxnTime()>=outputTimes[currentInterval]) {
			output.record(currentRealization,currentInterval,currentPopulation);
			currentInterval++;
		}	
	}//end main for
}//end simulate	

template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
inline
double
SSA_NRM<_populationVectorType, 
		   _stoichiometryType,
		   _propensitiesFunctorType,
		   _dependencyGraphType>::
selectStepSize(){
  return rxnHeap.getNextRxnTime();
}

template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
inline
bool
SSA_NRM<_populationVectorType, 
		   _stoichiometryType,
		   _propensitiesFunctorType,
		   _dependencyGraphType>::
setSingleSpeciesCurrentPopulation(std::size_t speciesIndex, populationValueType newPopulation)
{
	//set new population
	currentPopulation[speciesIndex]=newPopulation;

	std::vector<double> nextRxnTime;
	nextRxnTime.resize( NumberOfReactions );
	for (std::size_t i=0; i<NumberOfReactions; ++i)
  	{
  		nextRxnTime[i] = rxnHeap.getRxnTime(i);
  	}

  	//modulate the next rxn times for all reactions
  	//i.e. apply formula to newNextRxnTimes to account for changed propensities
  	for (std::size_t i=0; i!=NumberOfReactions; ++i)
  	{
  		double oldPropensity, newPropensity;
  			
	    oldPropensity=currentPropensities[i];

	    currentPropensities[i]=propensities(i,currentPopulation);
	    newPropensity = currentPropensities[i];

	    //handle reactions with positive propensity normally
	    if( newPropensity != 0.0 )
	    {
		    if( oldPropensity != 0.0 )
		     	nextRxnTime[i] = (oldPropensity/newPropensity) * nextRxnTime[i];

		    //draw new random number when adjustment
		    //doesn't work due to an oldPropensity value of zero
		    else
		     	nextRxnTime[i] = randomGenerator.Exponential( 1.0/newPropensity );
	    }

	    //handle reactions with 0 propensities by setting the next rxn time to infinity
	    else
	    {
	      	nextRxnTime[i] = std::numeric_limits<double>::infinity();
	    }
  	}

  	//update rxn time in rxn heap
  	for (std::size_t i=0; i<NumberOfReactions; ++i)
  	{
  		rxnHeap.setNewRxnTime(i, nextRxnTime[i]);
  	}

    return true;
}

template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
inline
bool
SSA_NRM<_populationVectorType, 
		   _stoichiometryType,
		   _propensitiesFunctorType,
		   _dependencyGraphType>::
setParameterValue(std::size_t parameterIndex, std::string newParameterValue, bool value_type)
{
	bool updateStatus;
	updateStatus = propensities.ParametersList.updateParameter(parameterIndex, newParameterValue, value_type);
	if (!updateStatus){
		return false; 
	}

	std::vector<double> nextRxnTime;
	nextRxnTime.resize( NumberOfReactions );
	for (std::size_t i=0; i<NumberOfReactions; ++i)
  	{
  		nextRxnTime[i] = rxnHeap.getRxnTime(i);
  	}

  	//modulate the next rxn times for all reactions
  	//i.e. apply formula to newNextRxnTimes to account for changed propensities
  	for (std::size_t i=0; i!=NumberOfReactions; ++i)
  	{
  		double oldPropensity, newPropensity;
  			
	    oldPropensity=currentPropensities[i];

	    currentPropensities[i]=propensities(i,currentPopulation);
	    newPropensity = currentPropensities[i];

	    //handle reactions with positive propensity normally
	    if( newPropensity != 0.0 )
	    {
		    if( oldPropensity != 0.0 )
		     	nextRxnTime[i] = (oldPropensity/newPropensity) * nextRxnTime[i];

		    //draw new random number when adjustment
		    //doesn't work due to an oldPropensity value of zero
		    else
		     	nextRxnTime[i] = randomGenerator.Exponential( 1.0/newPropensity );
	    }

	    //handle reactions with 0 propensities by setting the next rxn time to infinity
	    else
	    {
	      	nextRxnTime[i] = std::numeric_limits<double>::infinity();
	    }
  	}

  	//update rxn time in rxn heap
  	for (std::size_t i=0; i<NumberOfReactions; ++i)
  	{
  		rxnHeap.setNewRxnTime(i, nextRxnTime[i]);
  	}

	return true;
}


}//namespace Stochkit
