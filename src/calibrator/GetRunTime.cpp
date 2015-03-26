#include "GetRunTime.h"

double getODMRunTime(char* modelFileName, double simulationTime)
{

  Input_mass_action
    <StandardDriverTypes::populationType,
     StandardDriverTypes::stoichiometryType,
     StandardDriverTypes::propensitiesType, 
     StandardDriverTypes::graphType>
    model(modelFileName);
  
  ssa_odm solver(model.writeInitialPopulation(),
		 model.writeStoichiometry(),
		 model.writePropensities(),
		 model.writeDependencyGraph());
  
  return simulateSingleTrajectoryODM(solver, simulationTime);

}

double getDirectRunTime(char* modelFileName, double simulationTime)
{
  Input_mass_action
    <StandardDriverTypes::populationType,
     StandardDriverTypes::stoichiometryType,
     StandardDriverTypes::propensitiesType, 
     StandardDriverTypes::graphType>
    model(modelFileName);
  
  ssa_dir solver(model.writeInitialPopulation(),
	     model.writeStoichiometry(),
	     model.writePropensities(),
	     model.writeDependencyGraph());
  
  return simulateSingleTrajectoryDirect(solver, simulationTime);

}

double getConstantRunTime(char* modelFileName, double simulationTime)
{
  Input_mass_action
    <StandardDriverTypes::populationType,
     StandardDriverTypes::stoichiometryType,
     StandardDriverTypes::propensitiesType, 
     StandardDriverTypes::graphType>
    model(modelFileName);
  
  ssa_cst solver(model.writeInitialPopulation(),
	     model.writeStoichiometry(),
	     model.writePropensities(),
	     model.writeDependencyGraph());
  
  return simulateSingleTrajectoryDirect(solver, simulationTime);

}

double getNRMRunTime(char* modelFileName, double simulationTime)
{
  Input_mass_action
    <StandardDriverTypes::populationType,
     StandardDriverTypes::stoichiometryType,
     StandardDriverTypes::propensitiesType, 
     StandardDriverTypes::graphType>
    model(modelFileName);

  ssa_nrm solver(model.writeInitialPopulation(),
	     model.writeStoichiometry(),
	     model.writePropensities(),
	     model.writeDependencyGraph());
  
  return simulateSingleTrajectoryNRM(solver, simulationTime);

}
