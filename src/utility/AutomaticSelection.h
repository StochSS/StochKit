#ifndef _AUTOMATIC_SELECTION_H_
#define _AUTOMATIC_SELECTION_H_

/*
*  analyze model and choose which ssa method to use
*/

#include <iostream>
#include "boost_headers.h"

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
int SSA_AutomaticSelection(std::size_t numberOfReactions, std::size_t numberOfRealizations, int numberOfProcesses);

/*
 * Automatically select dense or sparse stoichiometry representation
 *
 * Return value:
 * True: use dense stoichiometry
 * False: use sparse stoichiometry
 */
bool denseStoichiometryOrNot(std::size_t numberOfSpecies);

}

#endif

