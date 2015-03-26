#ifndef _SIMULATE_SINGLE_TRAJECTORY_H_
#define _SIMULATE_SINGLE_TRAJECTORY_H_

#include <vector>
#include <utility>
#include <sys/time.h>
#include "StandardDriverTypes.h"

using namespace STOCHKIT;


/**
   Simulate single trajectory for ODM method.
*/

template <class SSA_ODMSolver_T>
double simulateSingleTrajectoryODM(SSA_ODMSolver_T& solver, double simulationTime) {	
  unsigned nReactionsFired=0;

  solver.presimulation(0.0, simulationTime);
  solver.initialize();
  solver.setCurrentTime(solver.getCurrentTime()+solver.selectStepSize());

  timeval timer;
  gettimeofday(&timer, NULL);
  double start_time=timer.tv_sec+(timer.tv_usec/1000000.0);
  
  while (solver.getCurrentTime() < simulationTime)
    {
      
      int rxnIndex=solver.selectReaction();
      solver.fireReaction(rxnIndex);
      nReactionsFired++;
      
      solver.setCurrentTime(solver.getCurrentTime()+solver.selectStepSize());
    }
  
  gettimeofday(&timer,NULL);
  double end_time=timer.tv_sec+(timer.tv_usec/1000000.0);
  double total_time=(end_time-start_time);//elapsed time, units=seconds
  
  return total_time/nReactionsFired;
}


/**
   Simulate single trajectory for certain direct method variants. This includes
   Constant, and Direct. This is the only single trajectory that can
   be used on the direct method variants.
 */
template <class SSA_Solver_T>
double simulateSingleTrajectoryDirect(SSA_Solver_T& solver, double simulationTime) {	
  unsigned nReactionsFired=0;

  solver.initialize();
  
  solver.setCurrentTime(solver.getCurrentTime()+solver.selectStepSize());
  
  timeval timer;
  gettimeofday(&timer, NULL);
  double start_time=timer.tv_sec+(timer.tv_usec/1000000.0);

  while (solver.getCurrentTime() < simulationTime)
    {
      
      int rxnIndex=solver.selectReaction();
      solver.fireReaction(rxnIndex);
      nReactionsFired++;
      
      solver.setCurrentTime(solver.getCurrentTime()+solver.selectStepSize());
    }
  
  gettimeofday(&timer,NULL);
  double end_time=timer.tv_sec+(timer.tv_usec/1000000.0);
  double total_time=(end_time-start_time);//elapsed time, units=seconds
  
  return total_time/nReactionsFired;
}

/**
   Simulate single trajectory for all NRM derived methods including RTC.
   Only these methods can use this function. This is because these methods
   have significantly different interfaces than the direct method derivatives.
 */
template <class SSA_Solver_T>
double simulateSingleTrajectoryNRM(SSA_Solver_T& solver, double simulationTime) {	
  unsigned nReactionsFired=0;

  solver.initialize();

  //timing
  timeval timer;
  gettimeofday(&timer, NULL);
  double start_time=timer.tv_sec+(timer.tv_usec/1000000.0);
  //timing
  
  while (solver.getCurrentTime() < simulationTime)
    { 
      solver.fireReaction();
      nReactionsFired++;
    }
  
  //timing
  gettimeofday(&timer,NULL);
  double end_time=timer.tv_sec+(timer.tv_usec/1000000.0);
  double total_time=(end_time-start_time);//elapsed time, units=seconds
  //timing

  return total_time/nReactionsFired;
}

#endif
