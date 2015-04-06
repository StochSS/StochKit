/*!
	\brief the logarithmic direct method of the Stochastic Simulation Algorithm (SSA)
*/

#if !defined(_SSA_LDM_IPP_)
#error This file is the implementation of SSA_LDM
#endif
#include "CommandPassAux.h"

namespace STOCHKIT
{

template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
void
SSA_LDM<_populationVectorType, 
		   _stoichiometryType,
		   _propensitiesFunctorType,
		   _dependencyGraphType>::
createDependencyGraphPostOrder(const _dependencyGraphType& depGraph) {
	for (std::size_t i=0; i!=depGraph.size(); ++i) {
		dependencyGraphPostOrder[i]=initialPropensityTree.sortPostOrder(depGraph[i]);
	}
}

template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
void
SSA_LDM<_populationVectorType, 
		   _stoichiometryType,
		   _propensitiesFunctorType,
		   _dependencyGraphType>::
initialize(double startTime) {
	  currentTime=startTime;
	  currentPopulation=initialPopulation;
	  propensityTree=initialPropensityTree;
}

template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
inline
double
SSA_LDM<_populationVectorType, 
		   _stoichiometryType,
		   _propensitiesFunctorType,
		   _dependencyGraphType>::
selectStepSize() {
	if (propensityTree.tree[0].sum<0.0) {
		std::cerr << "StochKit WARNING (SSA_LDM::selectStepSize): propensitySum<0, returning step size=infinity\n";
		return std::numeric_limits<double>::infinity();
	}
	return randomGenerator.Exponential(1.0/propensityTree.tree[0].sum);
}

template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
bool
SSA_LDM<_populationVectorType, 
		   _stoichiometryType,
		   _propensitiesFunctorType,
		   _dependencyGraphType>::
fireReaction() {

	if (propensityTree.tree[0].sum<=0.0) {
		return -1;
	}
	previousReactionIndex=propensityTree.selectReactionIndex(randomGenerator.ContinuousOpen(0.0, 1.0));

	int reactionIndex=previousReactionIndex;

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

	std::vector<double> updatedPropensities(dependencyGraphPostOrder[reactionIndex].size());
	for (std::size_t i=0; i!=updatedPropensities.size(); ++i) {
		updatedPropensities[i]=propensities(dependencyGraphPostOrder[reactionIndex][i],currentPopulation);
	}
	propensityTree.updateTree(dependencyGraphPostOrder[reactionIndex],updatedPropensities);

	return true;
}//end fireReaction

template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
template<typename IntervalOutputType>
void
SSA_LDM<_populationVectorType, 
		   _stoichiometryType,
		   _propensitiesFunctorType,
		   _dependencyGraphType>::
simulate(std::size_t realizations, double startTime, double endTime, IntervalOutputType& output, bool doValidate) {
	

	if (doValidate) {
		if (!validate(startTime,endTime)) {
			std::cerr << "StochKit ERROR (SSA_Direct::simulate): validate() failed, simulation aborted\n";
			exitFunc(1);
		}		
	}

	if (!output.initialize(realizations,startTime,endTime,initialPopulation)) {
		std::cerr << "StochKit ERROR (SSA_Direct::simulate): initialization of output object failed, simulation aborted\n";
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
			
			while (currentInterval<totalIntervals && currentTime >=outputTimes[currentInterval]){
				output.record(currentRealization,currentInterval,currentPopulation);
				currentInterval++;
			}

			fireReaction();
			currentTime+=selectStepSize();
		}
		while (currentInterval<totalIntervals && currentTime>=outputTimes[currentInterval]){
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
SSA_LDM<_populationVectorType, 
		   _stoichiometryType,
		   _propensitiesFunctorType,
		   _dependencyGraphType>::
setSingleSpeciesCurrentPopulation(std::size_t speciesIndex, populationValueType newPopulation)
{
	//set new population
	currentPopulation[speciesIndex]=newPopulation;

	calculateAllPropensities();
	std::vector<double> updatedPropensities;
	updatedPropensities.resize( currentPropensities.size() );
	for(std::size_t i=0; i<NumberOfReactions; i++)
	{
		updatedPropensities[i] = currentPropensities[i];
	}
	propensityTree.updateTree(updatedPropensities);

    return true;
}

template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
inline
bool
SSA_LDM<_populationVectorType, 
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

	calculateAllPropensities();
	std::vector<double> updatedPropensities;
	updatedPropensities.resize( currentPropensities.size() );
	for(std::size_t i=0; i<NumberOfReactions; i++)
	{
		updatedPropensities[i] = currentPropensities[i];
	}
	propensityTree.updateTree(updatedPropensities);

	return true;
}


}
