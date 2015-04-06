/*!
	\brief Text file input handler

	ALL IDS MUST START WITH AN ALPHABETIC LETTER, FOLLOWED BY LETTERS OR DIGITS
//read file
//parse by libxml2
//some check
//write data structure
//write to a file if need compile
//compile
*/

#ifndef _INPUT_H_
#define _INPUT_H_

#ifndef MAXPARAMETERNUM
#define MAXPARAMETERNUM 200 // maximum parameter number there could be
#endif

#define BADRESULT -32678 // indicate bad result such as divisor to be 0 or something in calculation

#include <iostream>
#include <sstream>
#include <algorithm>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include "boost/shared_ptr.hpp"
#include <vector>
#include <limits>
#include "Parameter.h"
#include "StringCalculator.h"
#include "VectorManipulation.h"
#include "CommandPassAux.h"

namespace STOCHKIT
{
 template<typename _populationVectorType,
	typename _stoichiometryType,
	typename _propensitiesFunctorType,
	typename _dependencyGraphType>
 class Input
 {
 protected:
	typedef typename _populationVectorType::value_type _populationValueType;

	bool WhiteSpaceinId;

	int NumberOfReactions, NumberOfSpecies;

	int NumberOfParameters;

	//! class to store parameters and related information
	ParameterSet ParametersList;
	
	//! class to store species and related information
	class Species{
		public:
			std::string Id; // Id of the species
			std::string InitialPopulation;  // initial population, only valid in SpeciesList
			std::vector<int> AffectReactions; // only valid in SpeciesList, record what reactions does one Species affect
	};
	std::vector<Species> SpeciesList;

	class SpeciesReference{
		public:
			std::string Id;
			int Index; // this species is SpeciesList[Index]
			int Stoichiometry; // used in Reactions
	};
	
	//! class to store reactions and related information
	class Reaction{
		public:
			std::string Id;  // Id of the reaction
			int Type; // type of the reaction, 0 = mass-action, 1 = michaelis-menten, 2 = customized
			std::string Rate; // for mass-action
			std::string Vmax;  // for michealis-menten
			std::string Km;  // for michealis-menten
			std::string Customized; // for customized propensity
			std::vector<SpeciesReference> Reactants;  // reactant list
			std::vector<SpeciesReference> Products;  // product list
	};
	std::vector<Reaction> ReactionsList;

	//! class to handle calculation of simple math expression strings
	StringCalculator simpleCalculator;

 public:

	Input(char *xmlFilename):
		simpleCalculator()
	{
		NumberOfReactions = 0;
		NumberOfSpecies = 0;
		NumberOfParameters = 0;
		WhiteSpaceinId = false;
		SpeciesList.clear();
		ReactionsList.clear();
		
		if(!parseXmlFile(xmlFilename)){
			exitFunc(1);
		}

                if(!ParametersList.linkParameters()){
			exitFunc(1);
		}

		if(!ParametersList.calculateParameters()){
			exitFunc(1);
		}

		if(!linkParametersAndReactions()){
			exitFunc(1);
		}

		if(!linkSpeciesAndReactions()){
			exitFunc(1);
		}

		if(!checkUniqueID()){
			exitFunc(1);
		}
	}

 protected:
	bool parseXmlFile(char *xmlFilename);

	bool recordNumberOfReactions(xmlNodePtr cur);

	bool recordNumberOfSpecies(xmlNodePtr cur);

	bool recordParametersList(xmlNodePtr cur);

	bool recordReactionsList(xmlNodePtr cur);

	bool recordSpeciesList(xmlNodePtr cur);

	bool linkParametersAndReactions();

	// record species reference order of reactants and products in reactions, meanwhile write affect reactions in species list
	bool linkSpeciesAndReactions();

	_populationValueType populationCalculation(std::string equation);

	//! check to see if all IDs are unique, return true if they are, return false if there are duplicate IDs
	bool checkUniqueID();

public:
	_populationVectorType writeInitialPopulation();

	_stoichiometryType writeStoichiometry();
	
	_dependencyGraphType writeDependencyGraph();
	
#ifdef JACOBIAN
	std::vector<std::vector<int> > writeReactants();
	std::vector<double> writeRates();
#endif
 }; // end of Input class
}

#define _INPUT_IPP_
#include "Input.ipp"
#undef _INPUT_IPP_

#endif

