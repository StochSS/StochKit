/*!
	\brief the amortized constant-complexity (composition-rejection) method of the Stochastic Simulation Algorithm (SSA)
*/

#if !defined(_SSA_CONSTANTTIME_IPP_)
#error This file is the implementation of SSA_ConstantTime
#endif
#include "CommandPassAux.h"

namespace STOCHKIT
{

template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
inline
double
SSA_ConstantTime<_populationVectorType, 
		   _stoichiometryType,
		   _propensitiesFunctorType,
		   _dependencyGraphType>::
selectStepSize() {
	if (groups.getPropensitySum()<0.0) {
	        //if propensitySum negative, recalculate all propensities
	        groups.recalculatePropensitySum();
		//if still negative, give warning and return infinity
		if (groups.getPropensitySum()<0.0) {
			std::cerr << "StochKit WARNING (SSA_ConstantTime::selectStepSize): propensitySum<0, returning step size=infinity\n";
			return std::numeric_limits<double>::infinity();
		}
	}
	return randomGenerator.Exponential(1.0/groups.getPropensitySum());
}

template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
void
SSA_ConstantTime<_populationVectorType, 
		   _stoichiometryType,
		   _propensitiesFunctorType,
		   _dependencyGraphType>::
initialize(double startTime) {
	previousReactionIndex=-1;
	currentTime=startTime;
	currentPopulation=initialPopulation;
	currentPropensities=initialPropensities;
	groups=initialGroups;
}

template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
int
SSA_ConstantTime<_populationVectorType, 
		   _stoichiometryType,
		   _propensitiesFunctorType,
		   _dependencyGraphType>::
selectReaction() {
	return groups.selectReaction(randomGenerator);
}//end selectReaction

template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
bool
SSA_ConstantTime<_populationVectorType, 
		   _stoichiometryType,
		   _propensitiesFunctorType,
		   _dependencyGraphType>::
fireReaction() {
        int reactionIndex = selectReaction();

	previousReactionIndex=reactionIndex;

	if (reactionIndex==-1) {
		return false;
	}

	#ifdef MATRIX_STOICHIOMETRY
		matrixrow dX(stoichiometry, reactionIndex);
		typename matrixrow::iterator it;
		for(it=dX.begin();it!=dX.end();it++) {
			currentPopulation[it.index()]+=*it;
		}
	#else
		currentPopulation+=stoichiometry[reactionIndex];
	#endif


	int affectedReactionIndex;
	double oldPropensity;

	//update affected reactions
	std::size_t numAffectedReactions=dependencyGraph[previousReactionIndex].size();
	for (std::size_t i=0; i!=numAffectedReactions; ++i) {
		affectedReactionIndex=dependencyGraph[previousReactionIndex][i];

		oldPropensity=currentPropensities[affectedReactionIndex];
		currentPropensities[affectedReactionIndex]=propensities(affectedReactionIndex,currentPopulation);

		//now we need to update groups data structures
		groups.update(affectedReactionIndex, oldPropensity, currentPropensities[affectedReactionIndex]);
	}

	return true;
}//end fireReaction

template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
template<typename IntervalOutputType>
void
SSA_ConstantTime<_populationVectorType, 
		   _stoichiometryType,
		   _propensitiesFunctorType,
		   _dependencyGraphType>::
simulate(std::size_t realizations, double startTime, double endTime, IntervalOutputType& output, bool doValidate) {

	if (doValidate) {
		if (!validate(startTime,endTime)) {
			std::cerr << "StochKit ERROR (SSA_ConstantTime::simulate): validate() failed, simulation aborted\n";
			exitFunc(1);
		}
	}

	if (!output.initialize(realizations,startTime,endTime,initialPopulation)) {
		std::cerr << "StochKit ERROR (SSA_ConstantTime::simulate): initialization of output object failed, simulation aborted\n";
		exitFunc(1);
	}
	
	std::vector<double> outputTimes = output.getOutputTimes();
	std::size_t totalIntervals=outputTimes.size();

	std::size_t currentInterval;

	for (std::size_t currentRealization=0; currentRealization!=realizations; ++currentRealization) {
		initialize(startTime);
		currentInterval=0;
		
		currentTime+=selectStepSize();
		while (currentTime<endTime) {	
			while (currentInterval<totalIntervals && currentTime>=outputTimes[currentInterval]){
				output.record(currentRealization,currentInterval,currentPopulation);
				currentInterval++;
			}
			fireReaction();
			currentTime+=selectStepSize();
		}
		while (currentInterval<totalIntervals && currentTime >=outputTimes[currentInterval]){
			output.record(currentRealization,currentInterval,currentPopulation);
			currentInterval++;
		}
	}
}//end simulate	

template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
inline
bool
SSA_ConstantTime<_populationVectorType, 
		   _stoichiometryType,
		   _propensitiesFunctorType,
		   _dependencyGraphType>::
setSingleSpeciesCurrentPopulation(std::size_t speciesIndex, populationValueType newPopulation)
{
	//set new population
	currentPopulation[speciesIndex]=newPopulation;

	double oldPropensity;

	//update affected reactions
	for (std::size_t i=0; i!=NumberOfReactions; ++i) {

		oldPropensity=currentPropensities[i];
		currentPropensities[i]=propensities(i,currentPopulation);

		//now we need to update groups data structures
		groups.update(i, oldPropensity, currentPropensities[i]);
	}

    return true;
}

template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
inline
bool
SSA_ConstantTime<_populationVectorType, 
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

	double oldPropensity;

	//update affected reactions
	for (std::size_t i=0; i!=NumberOfReactions; ++i) {

		oldPropensity=currentPropensities[i];
		currentPropensities[i]=propensities(i,currentPopulation);

		//now we need to update groups data structures
		groups.update(i, oldPropensity, currentPropensities[i]);
	}

	return true;
}



}
