/*
*  analyze model and choose which ssa method to use
*/

#include "AutomaticSelection.h"

namespace STOCHKIT
{

/*
 * Automatically select SSA method to use based on NumberOfReactions, NumberOfRealizations, and NumberOfProcesses
 *
 * Return value:
 * -1: ERROR
 *  0: SSA_Direct  (1:SSA_Direct with dense stoichiometry, not determined in this function)
 * 10: SSA_ODM     (11:SSA_ODM with dense stoichiometry, not determined in this function)
 * 30: SSA_ConstantTime
 * 40: SSA_NRM
 */
int SSA_AutomaticSelection(std::size_t numberOfReactions, std::size_t numberOfRealizations, int numberOfProcesses)
{

	//tests needed to determine Cutoffs
	std::size_t constantOverODMCutoff=6000;//if number of reactions is >=5000, use constnat time instead of ODM
	std::size_t constantOverDirectCutoff=2000;
	std::size_t realizationsODMoverDirect=5;//single processor

#ifdef WIN32
	int methodChoice = -1;
#else
	std::size_t methodChoice = -1;
#endif

	//number of realizations vs processors determines if ODM or Direct is better
	//if default number of processes, 0, is chosen, automatically determine
	if (numberOfProcesses==0) {
		numberOfProcesses=boost::thread::hardware_concurrency();
		if (numberOfProcesses==0) {
			numberOfProcesses=1;
		}
	}
	//changed next line to workaround bug in odm (odm crashes if model has just one reaction)
	if (numberOfRealizations>=realizationsODMoverDirect*numberOfProcesses && numberOfReactions>1) {
		methodChoice = 10;
	}
	else {
		methodChoice = 0;
	}

	//use number of reactions and realizations to determine if we should use constant time
	if (methodChoice==10) {
		if ( (numberOfReactions/constantOverDirectCutoff+numberOfReactions*numberOfRealizations/constantOverODMCutoff) > numberOfRealizations ) {
			methodChoice = 30;
		}
	}
	else {
		if (numberOfReactions>=constantOverDirectCutoff) {
			methodChoice = 30;
		}
	}

	return methodChoice;
}

/*
 * Automatically select dense or sparse stoichiometry representation
 *
 * Return value:
 * True: use dense stoichiometry
 * False: use sparse stoichiometry
 */
bool denseStoichiometryOrNot(std::size_t numberOfSpecies)
{
        std::size_t denseStoichiometryCutoff=64;//dense/sparse cutoff roughly this--doesn't make significant difference until much larger
	bool denseFlag = false;
	if (numberOfSpecies<=denseStoichiometryCutoff) {
		denseFlag = true;
	}

	return denseFlag;
}

}

