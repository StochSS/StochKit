/*!
	\brief Parameter class definition
*/
#include "Parameter.h"

namespace STOCHKIT
{

bool ParameterSet::calculateParameter(unsigned int index)
{
#ifdef DEBUG
	if( index >= ListOfParameters.size() ){
		std::cerr << "StochKit ERROR (Parameter::calculateParameter): calculate non-existing parameter" << std::endl;
		return false;
	}
#endif

	if(ListOfParameters[index].CalculateFlag == 1){
		return true;
	} else if(ListOfParameters[index].CalculateFlag == -1){
		ListOfParameters[index].CalculateFlag = 0;
		std::vector<unsigned int>::iterator para_it; // iterator of parameters in link graph

		for( para_it = ListOfParameters[index].ParametersAffectThis.begin(); para_it < ListOfParameters[index].ParametersAffectThis.end(); ++para_it ){
			if( ListOfParameters[*para_it].CalculateFlag != 1 ){
				if(!calculateParameter(*para_it)){
					return false;
				}
			}
		}

		std::string parameterExpression = parameterSubstitution(ListOfParameters[index].Expression);
		if( parameterExpression.empty() ){
			std::cerr << "StochKit ERROR (Parameter::calculateParameter): while calculating parameter " + ListOfParameters[index].Id << std::endl;
			return false;
		}

		ListOfParameters[index].Value = simpleCalculator.calculateString(parameterExpression);
		if(ListOfParameters[index].Value == BADRESULT){
			std::cerr << "StochKit ERROR (Parameter::calculateParameter): while calculating parameter " + ListOfParameters[index].Id << std::endl;
			return false;
		}

		ListOfParameters[index].CalculateFlag = 1;

		return true;
	} else { // CalculateFlag == 0
		std::cerr << "StochKit ERROR (Parameter::calculateParameter): there is a loop in parameters link graph" << std::endl;
		return false;
	}
}

bool ParameterSet::calculateParameters()
{
        // parameter substitution and calculate parameters
	for( unsigned int i = 0; i < ListOfParameters.size(); ++i ){
		if(!calculateParameter(i)){
			return false;
		}
	}

	return true;
}

// parameter substitution
// for example: 
// if P2 = 3, then  2*P2 -> 2*3
// returns NULL if parameter not found or there are too many parameters or there are some expression dead loops in parameters list
std::string ParameterSet::parameterSubstitution(std::string equation)
{
	// locate a parameter name in equation
	unsigned int begin = 0, end = 0;
	std::string substitutedEquation = equation;

	while( begin < substitutedEquation.length() ){
		if( (isalpha(substitutedEquation.at(begin)) || substitutedEquation.at(begin) == '_') && ( (begin==0) || ((substitutedEquation.at(begin) != 'e') && (substitutedEquation.at(begin) != 'E')) || !(isalnum(substitutedEquation.at(begin-1)) || substitutedEquation.at(begin-1) == '_') )){
			end = begin+1;
			while( (end<substitutedEquation.length()) && (isalnum(substitutedEquation.at(end)) || substitutedEquation.at(end) == '_') )
				++end;
			
			std::string parameterName = substitutedEquation.substr(begin,end-begin);
			
			// search for the parameter in ListOfParameters
			unsigned int i = 0;
			while( (i < ListOfParameters.size()) && (parameterName.compare(ListOfParameters[i].Id)!=0) )
				++i;
			
			if( i == ListOfParameters.size() ){
				unsigned int j = 0;
 				while( (j < knownFunctions.functions.size()) && (parameterName.compare(knownFunctions.functions[j].first)!=0) )
 					++j;
 				if( j == knownFunctions.functions.size() ){
 					std::cerr << "StochKit ERROR (Parameter::parameterSubstitution): parameter " + parameterName + " not found in parameters list\n";
 					substitutedEquation.clear();
 					return substitutedEquation;
 				}
 				else{
 					begin += knownFunctions.functions[j].first.length();
  				}
			}
			else{
				// substitute parameter with its value
				std::ostringstream parameterValue;
				parameterValue << ListOfParameters[i].Value;
				substitutedEquation.replace(begin, end-begin, parameterValue.str());
				begin += parameterValue.str().length();
			}
		}
		else{
			++begin;
		}
	}
	
	return substitutedEquation;
}
	
bool ParameterSet::linkParameters()
{
	std::vector<unsigned int>::iterator vec_it;

	for( unsigned int i=0; i<ListOfParameters.size(); ++i ){
		ListOfParameters[i].ParametersAffectThis = analyzeParameterExpression(ListOfParameters[i].Expression);
		
		for(vec_it = ListOfParameters[i].ParametersAffectThis.begin(); vec_it < ListOfParameters[i].ParametersAffectThis.end(); ++vec_it){
			insertToSortedArray<std::vector<unsigned int>, unsigned int>(ListOfParameters[*vec_it].AffectParameters, i);
		}
	}
	
	return true;
}

//! update ListOfParameters[index] to a new expression and calculate all related parameters
bool ParameterSet::updateParameter(unsigned int index, std::string Expression, bool value_type)
{
	std::vector<unsigned int> new_ParametersAffectThis;
	std::vector<unsigned int> newlyAffectingParameters;
	std::vector<unsigned int> stopAffectingParameters;
	std::vector<unsigned int>::iterator vec_it;

	// update desired parameter
	if( value_type == true){
		new_ParametersAffectThis.clear();
	} else{
		new_ParametersAffectThis = analyzeParameterExpression(Expression);
	}

	// update link graph
	stopAffectingParameters = vectorDifference(ListOfParameters[index].ParametersAffectThis, new_ParametersAffectThis);
	newlyAffectingParameters = vectorDifference(new_ParametersAffectThis, ListOfParameters[index].ParametersAffectThis);
	for(vec_it = stopAffectingParameters.begin(); vec_it < stopAffectingParameters.end(); ++vec_it){
		delFromSortedArray<std::vector<unsigned int>, unsigned int>(ListOfParameters[*vec_it].AffectParameters, index);
	}
	for(vec_it = newlyAffectingParameters.begin(); vec_it < newlyAffectingParameters.end(); ++vec_it){
		insertToSortedArray<std::vector<unsigned int>, unsigned int>(ListOfParameters[*vec_it].AffectParameters, index);
	}

	if( value_type == true ){
		std::string parameterExpression = parameterSubstitution(Expression);
		if( parameterExpression.empty() ){
			std::cerr << "StochKit ERROR (Parameter::updateParameter): while updating parameter " + ListOfParameters[index].Id << std::endl;
			return false;
		}

		ListOfParameters[index].Value = simpleCalculator.calculateString(parameterExpression);
		if(ListOfParameters[index].Value == BADRESULT){
			std::cerr << "StochKit ERROR (Parameter::updateParameter): while updating parameter " << ListOfParameters[index].Id << std::endl;
			return false;
		}

		std::ostringstream parameterValue;
		parameterValue << ListOfParameters[index].Value;

		ListOfParameters[index].Expression = parameterValue.str();
	} else{
		ListOfParameters[index].Expression = Expression;
	}
	
	ListOfParameters[index].ParametersAffectThis = new_ParametersAffectThis;

	// mark all parameters need to be calculated
	markParameterCalculationFlag(index);

	// calculate all parameters
	if(!calculateParameters()){
		std::cerr << "StochKit ERROR (Parameter::updateParameter): while updating parameter " << ListOfParameters[index].Id << std::endl;
		return false;
	}

	return true;
}

//! update parameter named "parameterId" to a new expression and calculate all related parameters
//return false if parameter name not found
bool ParameterSet::updateParameter(std::string parameterId, std::string Expression, bool value_type)
{
	std::size_t i=0;
	while(ListOfParameters[i].Id.compare(parameterId)!=0){
		++i;
	}
	if(i==ListOfParameters.size()){
		std::cerr << "StochKit ERROR (Parameter::updateParameter): parameter with Id \"" << parameterId << "\" cannot be found. " << std::endl;
		return false;
	} else {
		return updateParameter(i,Expression,value_type);
	}
}

//! mark ListOfParameters[index] and all the parameters affected by this parameter as needed to be calculated
bool ParameterSet::markParameterCalculationFlag(unsigned int index)
{
	ListOfParameters[index].CalculateFlag = -1;
	std::vector<unsigned int>::iterator para_it; // iterator of parameters in link graph

	for( para_it = ListOfParameters[index].AffectParameters.begin(); para_it < ListOfParameters[index].AffectParameters.end(); ++para_it ){
		markParameterCalculationFlag(*para_it);
	}

	return true;
}

//! Analyze parameter's expression and return a sorted vector of paramter orders that affect this expression
std::vector<unsigned int> ParameterSet::analyzeParameterExpression(std::string Expression)
{
	std::vector<unsigned int> PATE; // Parameters Affect This Expression

	// locate a parameter name in equation
	unsigned int begin = 0, end = 0;

	while( begin < Expression.length() ){
		if( (isalpha(Expression.at(begin)) || Expression.at(begin) == '_') && ( (begin==0) || ((Expression.at(begin) != 'e') && (Expression.at(begin) != 'E')) || !(isalnum(Expression.at(begin-1)) || Expression.at(begin-1) == '_' ))){
			end = begin+1;
			while( (end<Expression.length()) && (isalnum(Expression.at(end)) || Expression.at(end) == '_') )
				++end;
			
			std::string parameterName = Expression.substr(begin,end-begin);
			
			// search for the parameter in ListOfParameters
			unsigned int i = 0;
			while( (i < ListOfParameters.size()) && (parameterName.compare(ListOfParameters[i].Id)!=0) )
				++i;
			
			if( i == ListOfParameters.size() ){
				begin = end;				
			}
			else{
				// record parameter affecting this expression
				PATE.push_back(i);
				begin = end;
			}
		}
		else{
			++begin;
		}
	}

	sort(PATE.begin(),PATE.end());

	std::vector<unsigned int>::iterator PATE_it;

#ifdef WIN32
	if(!PATE.empty())
	{
		for( PATE_it = (PATE.begin()+1 ); PATE_it < PATE.end(); ++PATE_it ){
			if( *PATE_it == *(PATE_it-1) ){
				PATE.erase(PATE_it);
				--PATE_it;
			}
		}
	}
#else
	for( PATE_it = (PATE.begin() + 1); PATE_it < PATE.end(); ++PATE_it ){
		if( *PATE_it == *(PATE_it-1) ){
			PATE.erase(PATE_it);
			--PATE_it;
		}
	}
#endif

	return PATE;
}

//! to set isLinkedWithReactionsFlag to true
bool ParameterSet::finishLinkingWithReactions()
{
	isLinkedWithReactionsFlag = true;
	return isLinkedWithReactionsFlag;
}

//! to determine if all the parameters have been linked with reactions
bool ParameterSet::isLinkedWithReactions()
{
	return isLinkedWithReactionsFlag;
}

//! return the position of the parameter which matches given Id
//  -1 if not found
int ParameterSet::findParameterWithId(std::string givenId)
{
	int found = -1;
        for( unsigned int i = 0; i < ListOfParameters.size(); ++i ){
                if(ListOfParameters[i].Id.compare(givenId) == 0){
                        found = (int)i;
			break;
                }
        }
	return found;
}

}
