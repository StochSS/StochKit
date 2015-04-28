/**
  @file ssa_serial.cpp

  @brief main() function launched by ParallelIntervalSimulation for running Stochkit on a single processor

  -relies on function driver_handler to launch a simulation and write output using the SerialIntervalSimulationDriver class

  @see ParallelIntervalSimulation
  @see SerialIntervalSimulationDriver
 */
#include "boost_headers.h"

#include <iostream>
#include <string>
#include "StandardEventHandler.h"
#include "StandardDriverTypes.h"
#include "SerialIntervalSimulationDriver.h"
#include "SSA_Headers.h"
#ifdef USER_DRIVER
  #include "user_ensemble_simulate.h"
#else
  #include "ssa_ensemble_simulate.h"
  #include "ssa_ensemble_simulate_events.h"
#endif
#include "ssa_calibrate.h"
#include "ssa_calibrate_events.h"

#define STANDARD_SOLVER_MEMBER_TYPE_SPARSE StandardDriverTypes::populationType, StandardDriverTypes::stoichiometryType, StandardDriverTypes::propensitiesType, StandardDriverTypes::graphType
#define STANDARD_SOLVER_MEMBER_TYPE_DENSE StandardDriverTypes::populationType, StandardDriverTypes::denseStoichiometryType, StandardDriverTypes::propensitiesType, StandardDriverTypes::graphType

using namespace STOCHKIT;

template<typename _solverType, typename _eventHandlerType, typename _outputType> 
bool driver_handler(int ac, char* av[], bool calibrating);

int main(int ac, char* av[])
{
  typedef StandardEventHandler<StandardDriverTypes::populationType> eventHandlerType;
  typedef StandardDriverTypes::outputType outputType;
 
  CommandLineInterface commandLine(ac,av);

  int methodId = commandLine.getMethod();

  bool calibrating = commandLine.shouldCalibrate();

  bool simulationStatus=false;
// -1: ERROR
//  0: SSA_Direct  (1:SSA_Direct with dense stoichiometry)
// 10: SSA_ODM     (11:SSA_ODM with dense stoichiometry)
// 30: SSA_ConstantTime
// 40: SSA_NRM
  switch(methodId)
  {
     case 0:
       {
	typedef SSA_Direct<STANDARD_SOLVER_MEMBER_TYPE_SPARSE> solverType;
	simulationStatus = driver_handler<solverType, eventHandlerType, outputType>(ac, av, calibrating);
	break;
       }
     case 1:
       {
	typedef SSA_Direct<STANDARD_SOLVER_MEMBER_TYPE_DENSE> solverType;
	simulationStatus = driver_handler<solverType, eventHandlerType, outputType>(ac, av, calibrating);
	break;
       }
     case 10:
       {
	typedef SSA_ODM<STANDARD_SOLVER_MEMBER_TYPE_SPARSE> solverType;
	simulationStatus = driver_handler<solverType, eventHandlerType, outputType>(ac, av, calibrating);
	break;
       }
     case 11:
       {
	typedef SSA_ODM<STANDARD_SOLVER_MEMBER_TYPE_DENSE> solverType;
	simulationStatus = driver_handler<solverType, eventHandlerType, outputType>(ac, av, calibrating);
	break;
       }
    case 30:
       {
	typedef SSA_ConstantTime<STANDARD_SOLVER_MEMBER_TYPE_SPARSE> solverType;
  simulationStatus = driver_handler<solverType, eventHandlerType, outputType>(ac, av, calibrating);
	break;
       }
     case 40:
       {
	typedef SSA_NRM<STANDARD_SOLVER_MEMBER_TYPE_SPARSE> solverType;
	simulationStatus = driver_handler<solverType, eventHandlerType, outputType>(ac, av, calibrating);
	break;
       }
     default:
       {
          std::cerr << "StochKit ERROR: Solver type not recognized." << std::endl;    
          return -1;
       }
  }

  if(simulationStatus == false){
    std::cerr << "StochKit ERROR: Simulation was not successfully carried out." << std::endl;    
    return -1;
  }

  return 0;
}

/**
  -switch simply allows the selection of either ensemble_simulate function or calibrate function
  -cases do not differ after declaration of driver object except for writeOutput(), but we must have redundant code because 
  otherwise the driver object would not be at the correct scope
*/
template<typename _solverType, typename _eventHandlerType, typename _outputType>
bool driver_handler(int ac, char* av[], bool calibrating)
{
  switch( calibrating == true )
  {
    case 1:
    {
        #ifdef EVENTS
        SerialIntervalSimulationDriver<_solverType, _eventHandlerType, _outputType> driver(ac,av, calibrate_events<_solverType, _eventHandlerType, _outputType>);
        #else
        SerialIntervalSimulationDriver<_solverType, _eventHandlerType, _outputType> driver(ac,av, calibrate<_solverType, _outputType>);
        #endif

        _solverType *solverPtr=NULL;
        if(driver.isMassAction()){
      	solverPtr=new _solverType(driver.createMassActionSolver());
        } else if(driver.isMixed()){
      #ifdef MIXED
      	solverPtr=new _solverType(driver.createMixedSolver());
      #else
      	if(solverPtr!=NULL){
      		delete solverPtr;
      	}
      	std::cerr << "StochKit ERROR: Mixed solver needs to be compiled first." << std::endl;
      	return false;
      #endif
        } else if(driver.isEventsEnabled()) {
      #ifdef EVENTS
        solverPtr=new _solverType(driver.createEventsSolver());
      #else
        if(solverPtr!=NULL){
          delete solverPtr;
        }
        std::cerr << "StochKit ERROR: Mixed solver needs to be compiled first." << std::endl;
        return false;
      #endif

        } else {
      	if(solverPtr!=NULL){
      		delete solverPtr;
      	}
      	std::cerr << "StochKit ERROR: Model is other than mass-action, mixed, or events_enabled type. Please make sure there is no mistake in the model file." << std::endl;
      	return false;
        }
      	
        
        //set solver-specific parameters
        //none for ssa_direct

        driver.callSimulate(*solverPtr);
        driver.writeCalibratorOutput();

        if(solverPtr!=NULL){
      	delete solverPtr;
        }
    }

    default:
    {
        #ifdef EVENTS
        SerialIntervalSimulationDriver<_solverType, _eventHandlerType, _outputType> driver(ac,av, ensemble_simulate_events<_solverType, _eventHandlerType, _outputType>);
        #else
        SerialIntervalSimulationDriver<_solverType, _eventHandlerType, _outputType> driver(ac,av, ensemble_simulate<_solverType, _outputType>);
        #endif

        _solverType *solverPtr=NULL;
        if(driver.isMassAction()){
        solverPtr=new _solverType(driver.createMassActionSolver());
        } else if(driver.isMixed()){
      #ifdef MIXED
        solverPtr=new _solverType(driver.createMixedSolver());
      #else
        if(solverPtr!=NULL){
          delete solverPtr;
        }
        std::cerr << "StochKit ERROR: Mixed solver needs to be compiled first." << std::endl;
        return false;
      #endif
        } else if(driver.isEventsEnabled()) {
      #ifdef EVENTS
        solverPtr=new _solverType(driver.createEventsSolver());
      #else
        if(solverPtr!=NULL){
          delete solverPtr;
        }
        std::cerr << "StochKit ERROR: Mixed solver needs to be compiled first." << std::endl;
        return false;
      #endif

        } else {
        if(solverPtr!=NULL){
          delete solverPtr;
        }
        std::cerr << "StochKit ERROR: Model is other than mass-action, mixed, or events_enabled type. Please make sure there is no mistake in the model file." << std::endl;
        return false;
        }
        
        
        //set solver-specific parameters
        //none for ssa_direct

        driver.callSimulate(*solverPtr);

        driver.writeOutput();

        if(solverPtr!=NULL){
        delete solverPtr;
        }
    }
  }//end switch

  return true;
}//end driver_handler
