#include "GetRunTimeMixed.h"

double getODMRunTimeMixed(char* modelFileName, double simulationTime)
{

  Input_mixed_after_compile
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

double getDirectRunTimeMixed(char* modelFileName, double simulationTime)
{
  Input_mixed_after_compile
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

double getConstantRunTimeMixed(char* modelFileName, double simulationTime)
{
  Input_mixed_after_compile
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

double getNRMRunTimeMixed(char* modelFileName, double simulationTime)
{
  Input_mixed_after_compile
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
