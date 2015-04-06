/*!
	\brief the direct method of the Stochastic Simulation Algorithm (SSA)
*/

#if !defined(_SSA_BASE_IPP_)
#error This file is the implementation of SSA_Base
#endif
#include "CommandPassAux.h"

namespace STOCHKIT
{

template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
SSA_Base<_populationVectorType, 
		   _stoichiometryType,
		   _propensitiesFunctorType,
		   _dependencyGraphType>::
SSA_Base(const _populationVectorType& initialPop,
	const _stoichiometryType& stoich,
	const _propensitiesFunctorType& propensitiesFunctor,
	const _dependencyGraphType& depGraph,
	int seed) :
		initialPopulation(initialPop),
		stoichiometry(stoich),
		propensities(propensitiesFunctor),
		dependencyGraph(depGraph),
	        NumberOfSpecies(initialPop.size()),
		#ifdef MATRIX_STOICHIOMETRY
	        	NumberOfReactions(stoich.size1()),
			currentPropensities(stoichiometry.size1()),
		#else
	        	NumberOfReactions(stoich.size()),
			currentPropensities(stoichiometry.size()),
		#endif
		previousReactionIndex(-1),
		maxStepsCalculateAllPropensities(defaultMaxStepsCalculateAllPropensities),
		detectedVerySmallPropensity(false)
{
	randomGenerator.Seed(seed);
}

template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
void
SSA_Base<_populationVectorType, 
		   _stoichiometryType,
		   _propensitiesFunctorType,
		   _dependencyGraphType>::
seed(int seed) {
	randomGenerator.Seed(seed);
}

template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
bool
SSA_Base<_populationVectorType, 
		   _stoichiometryType,
		   _propensitiesFunctorType,
		   _dependencyGraphType>::
prepare(std::size_t /*realizations*/, double /*startTime*/, double /*endTime*/) {
	return true;
}

template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
bool
SSA_Base<_populationVectorType, 
		   _stoichiometryType,
		   _propensitiesFunctorType,
		   _dependencyGraphType>::
validate(double startTime, double endTime) {
	if (startTime>=endTime) {
		std::cout << "StochKit ERROR (SSA_Base::validate): startTime not before endTime\n";
		return false;
	}

	std::size_t N=initialPopulation.size();
	#ifdef MATRIX_STOICHIOMETRY
		std::size_t M=stoichiometry.size1();
	#else
		std::size_t M=stoichiometry.size();
	#endif
	if (N==0) {
		std::cout << "StochKit ERROR (SSA_Base::validate): initial population size=0\n";
		return false;
	}
	if (N!=NumberOfSpecies) {
		std::cout << "StochKit ERROR (SSA_Base::validate): Number of species does not equal initial population size\n";
		return false;
	}
	if (M!=NumberOfReactions) {
		std::cout << "StochKit ERROR (SSA_Base::validate): Number of reactions does not equal stoichiometry size\n";
		return false;
	}
	if (M!=propensities.size()) {
		std::cout << "StochKit ERROR (SSA_Base::validate): Number of reactions does not equal propensities size\n";
		return false;
	}

	//check initial populations are all non-negative
	for (std::size_t i=0; i!=NumberOfSpecies; ++i) {
		if (initialPopulation[i]<0) {
			std::cout << "StochKit ERROR (SSA_Base::validate): negative value detected in initial population\n";
			return false;
		}
	}

	//check that propensities, evaluated with initial population, are all non-negative
	for (std::size_t i=0; i!=NumberOfReactions; ++i) {
		if (propensities(i,initialPopulation)<0.0) {
			std::cout << "StochKit ERROR (SSA_Base::validate): negative propensity detected based on initial population\n";
			return false;
		}
	}
			
	return true;	  
}

template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
void
SSA_Base<_populationVectorType, 
		   _stoichiometryType,
		   _propensitiesFunctorType,
		   _dependencyGraphType>::
calculateAllPropensities() {
	double smallestNonzeroPropensity=std::numeric_limits<double>::max();
	
	propensitySum=0.0;
	for (std::size_t i=0; i!=NumberOfReactions; ++i) {
		currentPropensities[i]=propensities(i,currentPopulation);
		#ifdef DEBUG
				if (currentPropensities[i]!=currentPropensities[i]) {
					std::cout << "StochKit DEBUG (SSA_Base::calculateAllPropensities): detected 'NaN (not a number)' propensity for reaction index "<<i<<".\n";
					std::cout << "currentPopulation was:\n";
					for (std::size_t j=0; j!=currentPopulation.size(); ++j) {
						std::cout << "currentPopulation["<<j<<"]="<<currentPopulation[j]<<"\n";
					}
					std::cout << "Terminating.\n";
					exitFunc(1);
				}
				if (currentPropensities[i]==std::numeric_limits<double>::infinity()) {
					std::cout << "StochKit DEBUG (SSA_Base::calculateAllPropensities): detected 'infinity' propensity for reaction index "<<i<<".\n";
					std::cout << "currentPopulation was:\n";
					for (std::size_t j=0; j!=currentPopulation.size(); ++j) {
						std::cout << "currentPopulation["<<j<<"]="<<currentPopulation[j]<<"\n";
					}
					std::cout << "Terminating.\n";
					exitFunc(1);
				}
			if (currentPropensities[i]<0.0) {
				std::cout << "StochKit DEBUG (SSA_Base::calculateAllPropensities): detected negative propensity ("<<currentPropensities[i]<<") for reaction index "<<i<<"\n";
				std::cout << "currentPopulation was:\n";
				for (std::size_t j=0; j!=currentPopulation.size(); ++j) {
					std::cout << "currentPopulation["<<j<<"]="<<currentPopulation[j]<<"\n";
				}
				std::cout << "Terminating.\n";
				exitFunc(1);
			}
		#endif
		propensitySum+=currentPropensities[i];
		if (currentPropensities[i]>0.0 && currentPropensities[i]<smallestNonzeroPropensity) {
			smallestNonzeroPropensity=currentPropensities[i];
		}
	}
	stepsSinceCalculateAllPropensities=0;
	
	if (propensitySum>0.0 && smallestNonzeroPropensity/propensitySum<2E-10) { //per S.Mauch, M.Stalzer. (2009) "Efficient Formulations for Exact..."
		if (detectedVerySmallPropensity==false) {
			detectedVerySmallPropensity=true;
			std::cout << "StochKit WARNING (SSA_Base::calculateAllPropensities): detected very small propensity value, biased sampling of small propensity reactions may occur\n";
		}
	}
	#ifdef DEBUG
		//a reasonable place to check for possible time step inaccuracy
		if (propensitySum>0.0 && currentTime>2E21/propensitySum) { //per Mauch, Stalzer. (2009) "Efficient Formulations..."
			std::cout << "StochKit DEBUG (SSA_Base::calculateAllPropensities): ratio of average time step size to simulation currentTime is very small, may lead to step size inaccuracies\n";
		}
	#endif
}

template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
inline
double
SSA_Base<_populationVectorType, 
		   _stoichiometryType,
		   _propensitiesFunctorType,
		   _dependencyGraphType>::
getCurrentTime(){
	return currentTime;
}


template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
inline
bool
SSA_Base<_populationVectorType, 
		   _stoichiometryType,
		   _propensitiesFunctorType,
		   _dependencyGraphType>::
setCurrentTime(double newCurrentTime){
	currentTime = newCurrentTime;
	return true;
}

template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
inline
_populationVectorType
SSA_Base<_populationVectorType, 
		   _stoichiometryType,
		   _propensitiesFunctorType,
		   _dependencyGraphType>::
getCurrentPopulation(){
	return currentPopulation;
}

template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
inline
_populationVectorType
SSA_Base<_populationVectorType, 
		   _stoichiometryType,
		   _propensitiesFunctorType,
		   _dependencyGraphType>::
getInitialPopulation(){
	return initialPopulation;
}

template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
inline
bool
SSA_Base<_populationVectorType, 
		   _stoichiometryType,
		   _propensitiesFunctorType,
		   _dependencyGraphType>::
setCurrentPopulation(_populationVectorType& newPopulation){
#ifdef DEBUG
	if (newPopulation.size()!=NumberOfSpecies) {
		std::cerr<< "StochKit ERROR (SSA_Base::setCurrentPopulation): setting currentPopulation with greater or fewer species than NumberOfSpecies, simulation failure likely." << std::endl;
	}
	//loop over population to ensure non-negativity
	for (std::size_t i=0; i!=NumberOfSpecies; ++i) {
		if (newPopulation[i]<0) {
			std::cerr << "StochKit ERROR (SSA_Base::setCurrentPopulation): setting currentPopulation with one or more species with negative population, aborting." << std::endl;
			exitFunc(1);
  		}
	}
#endif
	currentPopulation=newPopulation;
	calculateAllPropensities();
	return true;
}

template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
inline
bool
SSA_Base<_populationVectorType, 
		   _stoichiometryType,
		   _propensitiesFunctorType,
		   _dependencyGraphType>::
setSingleSpeciesCurrentPopulation(std::size_t speciesIndex, populationValueType newPopulation)
{
    #ifdef DEBUG
		if (speciesIndex>=NumberOfSpecies) {
			std::cerr<<"StochKit ERROR (SSA_Base::setSingleSpeciesCurrentPopulation): attempt to setsetting currentPopulation with more species than NumberOfSpecies, will likely cause simulation failure." <<std::endl;
		}
		if (newPopulation<0) {
			std::cerr << "StochKit ERROR (SSA_Base::setSingleSpeciesCurrentPopulation): attempt to set species species population with a negative value, aborting." << std::endl;
			exitFunc(1);
		}
    #endif

    currentPopulation[speciesIndex]=newPopulation;

    calculateAllPropensities();
    return true;
}

		
template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
inline
typename SSA_Base<_populationVectorType, 
		   _stoichiometryType,
		   _propensitiesFunctorType,
		   _dependencyGraphType>::propensitiesVectorType
SSA_Base<_populationVectorType, 
		   _stoichiometryType,
		   _propensitiesFunctorType,
		   _dependencyGraphType>::
getCurrentPropensities(){
	return currentPropensities;
}

template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
inline
bool
SSA_Base<_populationVectorType, 
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
	return true;
}


template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
inline
int
SSA_Base<_populationVectorType, 
		   _stoichiometryType,
		   _propensitiesFunctorType,
		   _dependencyGraphType>::
getPreviousReactionIndex(){
	return previousReactionIndex;
}

template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
inline
ParameterSet&
SSA_Base<_populationVectorType, 
		   _stoichiometryType,
		   _propensitiesFunctorType,
		   _dependencyGraphType>::
referenceToParametersList(){
    return propensities.getReferenceToParametersList();
}

template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
inline
void
SSA_Base<_populationVectorType, 
		   _stoichiometryType,
		   _propensitiesFunctorType,
		   _dependencyGraphType>::
resetParameters(){
    propensities.reset();
}

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE

#ifndef BOOST_RANDOM_NO_STREAM_OPERATORS
template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
inline
std::string
SSA_Base<_populationVectorType, 
		   _stoichiometryType,
		   _propensitiesFunctorType,
		   _dependencyGraphType>::
getCurrentRNGState(){
	std::stringstream ss(std::stringstream::in | std::stringstream::out);
	ss.str("");
	ss.clear();
	ss << randomGenerator;
	return ss.str();
}


template<typename _populationVectorType, 
		 typename _stoichiometryType,
		 typename _propensitiesFunctorType,
		 typename _dependencyGraphType>
inline
bool
SSA_Base<_populationVectorType, 
		   _stoichiometryType,
		   _propensitiesFunctorType,
		   _dependencyGraphType>::
setCurrentRNGState(std::string RNGStateString){
	std::stringstream ss(RNGStateString, std::stringstream::in | std::stringstream::out);
	ss >> randomGenerator;
	return true;
}

#endif

#endif


}

