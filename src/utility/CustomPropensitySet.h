#if !defined(__CustomPropensitySet_h__)
#define __CustomPropensitySet_h__

#include <iostream>
#include <vector>
#include <algorithm>

#include "Parameter.h"
#include "CustomPropensity.h"
#include "CustomSimplePropensity.h"

namespace STOCHKIT
{
 template<typename _populationVectorType>
 class CustomPropensitySet
 {	
 protected:
  std::vector<CustomPropensity<_populationVectorType> > customPropensities;
  std::vector<CustomSimplePropensity<_populationVectorType> > simplePropensities;
  std::vector<CustomPropensity<_populationVectorType> *> propensities; // sorted propensities vector, currently only different than the originalPropensities in ODM
  std::vector<CustomPropensity<_populationVectorType> *> originalPropensities; // the original propensities vector read from input file
  std::vector<std::pair<unsigned int, unsigned int> > originalPropensities_index; // index pair (i,j) : i=0: simple; i=1: custom; j: position in corresponding vector
  std::vector<std::size_t> propensities_order;

  std::vector<std::string> simplePropensityRateConstantExpressions;

  //ParameterSet ParametersList;
  ParameterSet initialParametersList;

  //! class to handle calculation of simple math expression strings
  StringCalculator simpleCalculator;



 private:
        //! not allowing default constructor
	CustomPropensitySet()
	{
	}

        // calculate rate based on the value stored in parameterslist
        double rateCalculation(std::string equation)
        {
                std::vector<unsigned int> ParametersAffectRate;
                std::vector<unsigned int>::iterator para_it; // iterator of parameters in link graph

                ParametersAffectRate = ParametersList.analyzeParameterExpression(equation);

                bool calculationStatus = false;

                for( para_it = ParametersAffectRate.begin(); para_it < ParametersAffectRate.end(); ++para_it ){
                        if( ParametersList[*para_it].CalculateFlag == -1 ){
                                calculationStatus = ParametersList.calculateParameter(*para_it);
                                if(!calculationStatus){
                                        std::cerr << "StochKit ERROR (CustomPropensitySet::rateCalculation): while calculating rate " << equation << std::endl;
                                        return BADRESULT;
                                }
                        }
                }

                std::string substitutedEquation = ParametersList.parameterSubstitution(equation);
                if( substitutedEquation.empty() ){
                        std::cerr << "StochKit ERROR (CustomPropensitySet::rateCalculation): while calculating rate " << equation << std::endl;
                        return BADRESULT;
                }

                return simpleCalculator.calculateString(substitutedEquation);
        }


 public:

 	//Arya 7/23/13
 	//added because reference to private member is not allowed
 	ParameterSet ParametersList;

	CustomPropensitySet(std::vector<Parameter> existingParametersList):
    		initialParametersList(existingParametersList),
		simpleCalculator()
	{
		ParametersList = initialParametersList;
#ifdef DEBUG
		std::cout << "CustomPropensitySet Parameter Values:\n";
		for( std::size_t i = 0; i < ParametersList.size(); ++i){      
			std::cout << ParametersList[i].Value << " ";
		}
		std::cout << std::endl;
#endif
	}

	CustomPropensitySet(ParameterSet &existingParametersList):
    		initialParametersList(existingParametersList),
		simpleCalculator()
	{
		ParametersList = initialParametersList;
#ifdef DEBUG
		std::cout << "CustomPropensitySet Parameter Values:\n";
		for( std::size_t i = 0; i < ParametersList.size(); ++i){      
			std::cout << ParametersList[i].Value << " ";
		}
		std::cout << std::endl;
#endif
	}

	double operator()(const int n, _populationVectorType& populations) {
	  return (*propensities[n])(populations, ParametersList);
	}
	
	std::size_t size() {
		return propensities.size();
	}

	//! default destructor ok
//	~CustomPropensitySet() {
//	}

	//! copy-constructor
	CustomPropensitySet(const CustomPropensitySet& other)
	{
		ParametersList = other.ParametersList;
		initialParametersList = other.initialParametersList;
		customPropensities = other.customPropensities;
		simplePropensities = other.simplePropensities;
		originalPropensities_index = other.originalPropensities_index;
		propensities_order = other.propensities_order;

		//added 12/17/13 to fix bug where events models with type 0 reactions would not run- Arya
		simplePropensityRateConstantExpressions = other.simplePropensityRateConstantExpressions;

		originalPropensities.clear();
		// re-direct pointers
		for(unsigned int i=0; i < originalPropensities_index.size(); ++i){
			if(originalPropensities_index[i].first == 0){
				originalPropensities.push_back(&simplePropensities[originalPropensities_index[i].second]);
			} else {
				originalPropensities.push_back(&customPropensities[originalPropensities_index[i].second]);
			}
		}

		propensities.clear();
                for(std::size_t i=0;i<propensities_order.size();++i){
                        propensities.push_back(originalPropensities[propensities_order[i]]);
                }
	}

	//! assignment operator
	CustomPropensitySet& operator=(const CustomPropensitySet& other)
	{
		if(this != &other){ // protect against invalid self-assignment
			ParametersList = other.ParametersList;
			initialParametersList = other.initialParametersList;
			customPropensities = other.customPropensities;
			simplePropensities = other.simplePropensities;
			originalPropensities_index = other.originalPropensities_index;
			propensities_order = other.propensities_order;

			originalPropensities.clear();
			// re-direct pointers
			for(unsigned int i=0; i < originalPropensities_index.size(); ++i){
				if(originalPropensities_index[i].first == 0){
					originalPropensities.push_back(&simplePropensities[originalPropensities_index[i].second]);
				} else {
					originalPropensities.push_back(&customPropensities[originalPropensities_index[i].second]);
				}
			}			

			propensities.clear();
        	        for(std::size_t i=0;i<propensities_order.size();++i){
                	        propensities.push_back(originalPropensities[propensities_order[i]]);
	                }
		}
		return *this;
	}

	bool reset()
	{
		ParametersList = initialParametersList;
		for(std::size_t i=0;i<originalPropensities_index.size();++i){
			if(updatePropensity(i) == false){
				std::cerr << "StochKit ERROR (CustomPropensitySet::reset): while updating propensty values while resetting parameters list. " << std::endl;
				return false;
			}
		}
		return true;
	}

	bool resetOrder()
	{
		std::sort(propensities_order.begin(),propensities_order.end());
		propensities = originalPropensities;
		return true;
	}

	bool pushSimplePropensity(std::string rateExpression)
	{
		propensities_order.push_back(originalPropensities_index.size());
		originalPropensities_index.push_back(std::pair<unsigned int, unsigned int>(0,simplePropensities.size()));
		simplePropensityRateConstantExpressions.push_back(rateExpression);
		double rate =  rateCalculation(rateExpression);
		if( rate == BADRESULT ){
			std::cerr << "StochKit ERROR (CustomPropensitySet::pushSimplePropensity): while calculating rate expression " << rateExpression << std::endl;
			return false;
		}
		simplePropensities.push_back(CustomSimplePropensity<_populationVectorType>(rate));
		originalPropensities.push_back(&simplePropensities.back());
		propensities.push_back(&simplePropensities.back());
#ifdef DEBUG
		std::cout << "Reaction #" << originalPropensities_index.size() << " pushed." << std::endl;
#endif

		return true;
	}

	bool pushSimplePropensity(std::string rateExpression, int reactant1)
	{
		propensities_order.push_back(originalPropensities_index.size());
		originalPropensities_index.push_back(std::pair<unsigned int, unsigned int>(0,simplePropensities.size()));
		simplePropensityRateConstantExpressions.push_back(rateExpression);

		double rate =  rateCalculation(rateExpression);
		if( rate == BADRESULT ){
			std::cerr << "StochKit ERROR (CustomPropensitySet::pushSimplePropensity): while calculating rate expression " << rateExpression << std::endl;
			return false;
		}
		simplePropensities.push_back(CustomSimplePropensity<_populationVectorType>(rate, reactant1));
		originalPropensities.push_back(&simplePropensities.back());
		propensities.push_back(&simplePropensities.back());
#ifdef DEBUG
		std::cout << "Reaction #" << originalPropensities_index.size() << " pushed." << std::endl;
#endif

		return true;
	}

	bool pushSimplePropensity(std::string rateExpression, int reactant1, int reactant2)
	{
		propensities_order.push_back(originalPropensities_index.size());
		originalPropensities_index.push_back(std::pair<unsigned int, unsigned int>(0,simplePropensities.size()));
		simplePropensityRateConstantExpressions.push_back(rateExpression);
		double rate =  rateCalculation(rateExpression);
		if( rate == BADRESULT ){
			std::cerr << "StochKit ERROR (CustomPropensitySet::pushSimplePropensity): while calculating rate expression " << rateExpression << std::endl;
			return false;
		}
		simplePropensities.push_back(CustomSimplePropensity<_populationVectorType>(rate, reactant1, reactant2));
		originalPropensities.push_back(&simplePropensities.back());
		propensities.push_back(&simplePropensities.back());
#ifdef DEBUG
		std::cout << "Reaction #" << originalPropensities_index.size() << " pushed." << std::endl;
#endif

		return true;
	}

	bool pushSimplePropensity(std::string rateExpression, int reactant1, int reactant2, int reactant3)
	{
		propensities_order.push_back(originalPropensities_index.size());
		originalPropensities_index.push_back(std::pair<unsigned int, unsigned int>(0,simplePropensities.size()));
		simplePropensityRateConstantExpressions.push_back(rateExpression);
		double rate =  rateCalculation(rateExpression);
		if( rate == BADRESULT ){
			std::cerr << "StochKit ERROR (CustomPropensitySet::pushSimplePropensity): while calculating rate expression " << rateExpression << std::endl;
			return false;
		}
		simplePropensities.push_back(CustomSimplePropensity<_populationVectorType>(rate, reactant1, reactant2, reactant3));
		originalPropensities.push_back(&simplePropensities.back());
		propensities.push_back(&simplePropensities.back());
#ifdef DEBUG
		std::cout << "Reaction #" << originalPropensities_index.size() << " pushed." << std::endl;
#endif

		return true;
	}

	bool pushCustomPropensity(double (*CustomPropensityFunc)(_populationVectorType&, ParameterSet&))
	{
		propensities_order.push_back(originalPropensities_index.size());
		originalPropensities_index.push_back(std::pair<unsigned int, unsigned int>(1,customPropensities.size()));
		customPropensities.push_back(CustomPropensity<_populationVectorType>(CustomPropensityFunc));
		originalPropensities.push_back(&customPropensities.back());
		propensities.push_back(&customPropensities.back());
#ifdef DEBUG
		std::cout << "Reaction #" << originalPropensities_index.size() << " pushed." << std::endl;
#endif

		return true;
	}

	// update propensity of reaction n in the original order
	bool updatePropensity(unsigned int n)
	{
		if(originalPropensities_index[n].first==0) // simple propensity
		{
			std::string rateExpression = simplePropensityRateConstantExpressions[originalPropensities_index[n].second];
			double rate =  rateCalculation(rateExpression);
			if( rate == BADRESULT ){
				std::cerr << "StochKit ERROR (CustomPropensitySet::updateSimplePropensity): while calculating rate expression " << rateExpression << std::endl;
				return false;
			}
			simplePropensities[originalPropensities_index[n].second].updateRateConstant(rate);
		}
		return true;
	}

	bool updateParameter(unsigned int index, std::string Expression, bool value_type)
	{
		if(ParametersList.isLinkedWithReactions() == false){
	                std::cerr << "StochKit ERROR (CustomPropensitySet::updateParameter): ParameterSet was not successfully linked with reactions, thus it is not able to update propensities. Please linke parameters with reactions properly." << std::endl;
			return false;
		}

		ParametersList.updateParameter(index, Expression, value_type);
		for(unsigned int i=0;i<ParametersList[index].AffectReactions.size();++i){
			updatePropensity(ParametersList[index].AffectReactions[i]);
		}
		return true;
	}

	bool updateParameter(std::string parameterId, std::string Expression, bool value_type)
	{
	        std::size_t i=0;
	        while(ParametersList[i].Id.compare(parameterId)!=0){
	                ++i;
	        }
	        if(i==ParametersList.size()){
	                std::cerr << "StochKit ERROR (CustomPropensitySet::updateParameter): parameter with Id \"" << parameterId << "\" cannot be found. " << std::endl;
	                return false;
	        } else {
	                return updateParameter(i,Expression,value_type);
		}
        }

	double getParameterValue(unsigned int index)
	{
		if(ParametersList.calculateParameter(index) == false){
	                std::cerr << "StochKit ERROR (CustomPropensitySet::getParameterValue): while calculating the value of Parameter[" << index << ". " << std::endl;
			return 0.0; 
		} else {
			return ParametersList[index].Value;
		}
	}

	double getParameterValue(std::string parameterId)
	{
	        std::size_t i=0;
	        while(ParametersList[i].Id.compare(parameterId)!=0){
	                ++i;
	        }
	        if(i==ParametersList.size()){
	                std::cerr << "StochKit ERROR (CustomPropensitySet::getParameterValue): parameter with Id \"" << parameterId << "\" cannot be found. " << std::endl;
	                return false;
	        } else {
			ParametersList.calculateParameter(i);
	                return ParametersList[i].Value;
		}
        }


    ParameterSet& getReferenceToParametersList()
    {
    	return ParametersList;
    }

	/*
	 *  reOrder Propensities according to newOrder
	 *  e.g.: newOrder: {2,1,3}
	 *  then newPropensities: {originalReaction2, originalReaction1, originalReaction3}
	 *
	 */
	bool reOrderPropensities(const std::vector<std::size_t>& newOrder)
	{
		if(propensities_order.size() != newOrder.size()){
	                std::cerr << "StochKit ERROR (CustomPropensitySet::reOrderPropensities): the size of new order vector is different from the size of propensities vector. " << std::endl;
	                return false;
		}
#ifdef DEBUG
		std::vector<std::size_t> sortedNewOrder(newOrder);
		std::sort(sortedNewOrder.begin(),sortedNewOrder.end());
		for(std::size_t i=0;i<sortedNewOrder.size();++i){
			if(sortedNewOrder[i]!=i){
	                	std::cerr << "StochKit ERROR (CustomPropensitySet::reOrderPropensities): please make sure the newOrder vector is a permutation of 0,1,...,m-1 , assuming there are m reactions. " << std::endl;
		                return false;
			}
		}
#endif
		propensities_order = newOrder;
		for(std::size_t i=0;i<propensities_order.size();++i){
			propensities[i] = originalPropensities[propensities_order[i]];
		}

		return true;
	}
 };
}

#endif
