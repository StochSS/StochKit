#if !defined(_SSA_DIRECT_IPP_)
#error This file is the implementation of SSA_Direct
#endif

#include "CommandPassAux.h"

namespace STOCHKIT
{

template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
void
SSA_Direct<_populationVectorType, 
		   _stoichiometryType,
		   _propensitiesFunctorType,
		   _dependencyGraphType>::
initialize(double startTime) {
	previousReactionIndex=-1;
	currentTime=startTime;
	currentPopulation=initialPopulation;
	calculateAllPropensities();
}

template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
inline
double
SSA_Direct<_populationVectorType, 
		   _stoichiometryType,
		   _propensitiesFunctorType,
		   _dependencyGraphType>::
selectStepSize() {

	if (propensitySum<0.0) {
			#ifdef DEBUG
				std::cout << "StochKit DEBUG (SSA_Direct::selectStepSize): detected negative propensitySum, recalculating\n";
			#endif
	        //if propensitySum negative, recalculate all propensities
	        calculateAllPropensities();
		//if still negative, give warning and return infinity
		if (propensitySum<0.0) {
			std::cerr << "StochKit WARNING (SSA_Direct::selectStepSize): propensitySum<0, returning step size=infinity\n";
			return std::numeric_limits<double>::infinity();
		}
	}

	return randomGenerator.Exponential(1.0/propensitySum);
}

template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
inline
int
SSA_Direct<_populationVectorType, 
		   _stoichiometryType,
		   _propensitiesFunctorType,
		   _dependencyGraphType>::
selectReaction() {

	previousReactionIndex=-1;
	if (stepsSinceCalculateAllPropensities>maxStepsCalculateAllPropensities) {
		calculateAllPropensities();
	}

	//generate a uniform random number between (0,propensitySum)
	double r=0;
	while (r==0) {
		r=randomGenerator.ContinuousOpen(0,1)*propensitySum;
	}
	double jsum=0;
	while (jsum < r) {
		++previousReactionIndex;
		//test that we don't run off end of array
		if (previousReactionIndex==(int)NumberOfReactions) {
			#ifdef DEBUG
				std::cout << "StochKit DEBUG (SSA_Direct::selectReaction): detected numerical error in propensities, recalculating\n";
			#endif
			calculateAllPropensities();
			return selectReaction();
		}
		else {
			#ifdef DEBUG
				if (currentPropensities[previousReactionIndex]<0.0) {
					std::cout << "StochKit DEBUG (SSA_Direct::selectReaction): detected negative propensity ("<<currentPropensities[previousReactionIndex]<<") for reaction index "<<previousReactionIndex<<"\n";
					std::cout << "currentPopulation was:\n";
					for (std::size_t i=0; i!=currentPopulation.size(); ++i) {
						std::cout << "currentPopulation["<<i<<"]="<<currentPopulation[i]<<"\n";
					}
					std::cout << "Terminating.\n";
					exitFunc(1);
				}
			#endif			
			jsum+=currentPropensities[previousReactionIndex];
		}
	}
  
	return previousReactionIndex;
}

template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
inline
bool
SSA_Direct<_populationVectorType, 
		   _stoichiometryType,
		   _propensitiesFunctorType,
		   _dependencyGraphType>::
fireReaction() {
        int reactionIndex = selectReaction();

	if (reactionIndex==-1) {
		#ifdef DEBUG
			std::cout << "StochKit DEBUG (SSA_Direct::fireReaction): attempt to fire reaction index = -1. Terminating.\n";
			exitFunc(1);
		#endif
		return false;
	}
	else {
		//if not -1, assumes valid reactionIndex
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
		for (std::size_t i=0; i!=dependencyGraph[reactionIndex].size(); ++i) {
			affectedReactionIndex=dependencyGraph[reactionIndex][i];
			oldPropensity=currentPropensities[affectedReactionIndex];
			currentPropensities[affectedReactionIndex]=propensities(affectedReactionIndex,currentPopulation);
			#ifdef DEBUG
				if (currentPropensities[affectedReactionIndex]!=currentPropensities[affectedReactionIndex]) {
					std::cout << "StochKit DEBUG (SSA_Direct::fireReaction): detected 'NaN (not a number)' propensity for reaction index "<<affectedReactionIndex<<" while updating after firing reaction index "<<reactionIndex<<".\n";
				}
				if (currentPropensities[affectedReactionIndex]==std::numeric_limits<double>::infinity()) {
					std::cout << "StochKit DEBUG (SSA_Direct::fireReaction): detected 'infinity' propensity for reaction index "<<affectedReactionIndex<<" while updating after firing reaction index "<<reactionIndex<<".\n";
				}
				if (currentPropensities[affectedReactionIndex]<0.0) {
					std::cout << "StochKit DEBUG (SSA_Direct::fireReaction): detected negative propensity for reaction index "<<affectedReactionIndex<<" while updating after firing reaction index "<<reactionIndex<<".\n";
					std::cout << "updated currentPopulation is:\n";
					for (std::size_t i=0; i!=currentPopulation.size(); ++i) {
						std::cout << "currentPopulation["<<i<<"]="<<currentPopulation[i]<<"\n";
					}
					std::cout << "population before firing reaction was:\n";
					#ifdef MATRIX_STOICHIOMETRY
			                        matrixrow dX(stoichiometry, reactionIndex);
			                        typename matrixrow::iterator it;
			                        for(it=dX.begin();it!=dX.end();it++) {
                        			        currentPopulation[it.index()]-=*it;
			                        }
						_populationVectorType oldPop = currentPopulation;
					#else
						_populationVectorType oldPop=currentPopulation-=stoichiometry[reactionIndex];
					#endif
					for (std::size_t i=0; i!=oldPop.size(); ++i) {
						std::cout << "oldPopulation["<<i<<"]="<<oldPop[i]<<"\n";
					}
					std::cout << "Terminating.\n";
					exitFunc(1);
				}
			#endif			
			propensitySum+=currentPropensities[affectedReactionIndex]-oldPropensity;
		}
		stepsSinceCalculateAllPropensities++;
		return true;
	}
}

template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
template<typename IntervalOutputType>
void
SSA_Direct<_populationVectorType, 
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
		
}

