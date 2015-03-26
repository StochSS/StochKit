//getRunTime functions for mixed models. This will get compiled only at run time
//Note: only difference is the type of the model (Input_mass_action vs. Input_mixed_after_compile)


#ifndef _GET_RUN_TIME_MIXED_H_
#define _GET_RUN_TIME_MIXED_H_

#include <iostream>
#include <string>

#include "StandardDriverTypes.h"

#include "SSA_ODM.h"
#include "SSA_Direct.h"
#include "SSA_ConstantTime.h"
#include "SSA_NRM.h"

#include "Input_mixed_after_compile.h"

#include "SimulateSingleTrajectory.h"

typedef SSA_ODM<StandardDriverTypes::populationType,
		   StandardDriverTypes::stoichiometryType,
		   StandardDriverTypes::propensitiesType,
		   StandardDriverTypes::graphType> ssa_odm;

typedef SSA_Direct<StandardDriverTypes::populationType,
		   StandardDriverTypes::stoichiometryType,
		   StandardDriverTypes::propensitiesType,
		   StandardDriverTypes::graphType> ssa_dir;

typedef SSA_ConstantTime<StandardDriverTypes::populationType,
		StandardDriverTypes::stoichiometryType,
		StandardDriverTypes::propensitiesType,
		StandardDriverTypes::graphType> ssa_cst;

typedef SSA_NRM<StandardDriverTypes::populationType,
		StandardDriverTypes::stoichiometryType,
		StandardDriverTypes::propensitiesType,
		StandardDriverTypes::graphType> ssa_nrm;

//return average time to fire a single reaction over all reactions fired
//NOTE:separate functions are needed for the different solvers because the
//interface to each respective solver is different and thus the method for 
//firing a single trajectory varies from solver to solver 

double getODMRunTimeMixed(char* modelFileName, double simulationTime);
double getDirectRunTimeMixed(char* modelFileName, double simulationTime);
double getConstantRunTimeMixed(char* modelFileName, double simulationTime);
double getNRMRunTimeMixed(char* modelFileName, double simulationTime);

#endif

