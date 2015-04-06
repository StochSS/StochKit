/*!
	\brief Parameter class definition
*/

#ifndef _PARAMETER_H_
#define _PARAMETER_H_

#ifndef MAXPARAMETERNUM
#define MAXPARAMETERNUM 200 // maximum parameter number there could be
#endif

#ifndef BADRESULT
#define BADRESULT -32678 // indicate bad result such as divisor to be 0 or something in calculation
#endif

#include <iostream>
#include <sstream>
#include <algorithm>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include <boost/shared_ptr.hpp>
#include <vector>
#include <limits>
#include "VectorManipulation.h"
#include "SupportedFunctions.h"
#include "StringCalculator.h"

namespace STOCHKIT
{
 //! class to store parameters and related information
 class Parameter{
	public:
		std::string Id;  // Id of the parameter
		int Type; // 0 = constant, 1 = changeable during simulation
		std::string Expression;  // expression of the parameter
		int CalculateFlag; // -1 = not calculated yet, 0 = calculating, 1 = already calculated
		double Value; // value of the parameter
		std::vector<unsigned int> AffectParameters;
		std::vector<unsigned int> ParametersAffectThis;
		std::vector<unsigned int> AffectReactions;

		Parameter():
			CalculateFlag(-1)
		{
		}

		//! default copy constructor ok

		//! default destructor ok
 };

 //! class to store parameter lists and related operation
 class ParameterSet{
	protected:
		//! class to handle calculation of simple math expression strings
		StringCalculator simpleCalculator;

		SupportedFunctions knownFunctions;

		bool isLinkedWithReactionsFlag;

		std::vector<Parameter> ListOfParameters;

	public:
		ParameterSet():
			isLinkedWithReactionsFlag(false)
		{
			ListOfParameters.clear();
		}

		//! constructor with a list as arguement
		ParameterSet(std::vector<Parameter> existingListOfParameters):
			isLinkedWithReactionsFlag(false)
		{
			ListOfParameters = existingListOfParameters;
		}

		//! default copy constructor ok

		//! default destructor ok

		Parameter& operator[](unsigned int i){
			return ListOfParameters[i];
		}

		unsigned int size(){
			return ListOfParameters.size();
		}

		void push_back(const Parameter& x){
			ListOfParameters.push_back(x);
		}

		Parameter& back(){
			return ListOfParameters.back();
		}

		//! caluclate the ListOfParameters[index]
		bool calculateParameter(unsigned int index);
			
		bool calculateParameters();

		// parameter substitution
		// for example: 
		// if P2 = 3, then  2*P2 -> 2*3
		// returns NULL if parameter not found or there are too many parameters or there are some expression dead loops in parameters list
		std::string parameterSubstitution(std::string equation);

		//! Analyze parameter's expression and return a sorted vector of paramter orders that affect this expression
		std::vector<unsigned int> analyzeParameterExpression(std::string Expression);

		//! link parameters with each other
		bool linkParameters();

		//! update ListOfParameters[index] to a new expression and calculate all related parameters
		// value_type = true: value, will be constant afterwards
		// value_type = false: expression, still depend on other parameters
		bool updateParameter(unsigned int index, std::string Expression, bool value_type);

		//! update parameter named "parameterId" to a new expression and calculate all related parameters
		// value_type = true: value, will be constant afterwards
		// value_type = false: expression, still depend on other parameters
		// return false if parameter not found
		bool updateParameter(std::string parameterId, std::string Expression, bool value_type);

		//! mark ListOfParameters[index] and all the parameters affected by this parameter as needed to be calculated
		bool markParameterCalculationFlag(unsigned int index);

		//! to set isLinkedWithReactionsFlag to true
		bool finishLinkingWithReactions();

		//! to determine if all the parameters have been linked with reactions
		bool isLinkedWithReactions();

		//! return the position of the parameter which matches given Id
		//  -1 if not found
		int findParameterWithId(std::string givenId);
 }; // end of ParameterSet class
}

#endif

