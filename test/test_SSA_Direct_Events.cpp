/*
 *  test_SSA_Direct.cpp
 *  
 */
#include<iostream>
#include<string>
#include <fstream>
#include <sstream>
#include <vector>
#include <boost/numeric/ublas/vector.hpp>
#include<stdio.h>
#include <boost/shared_ptr.hpp>


#include "SSA_Direct_Events.h"
#include "TestModel.h"
#include "StatsAndTrajectoriesOutput.h"

#include "TimeBasedTrigger.h"
#include "CustomStateBasedTrigger.h"
#include "CustomStateBasedTriggerFunctions.h"
#include "CustomChangeSingleSpeciesFunctions.h"
#include "ChangeSingleSpeciesEventAction.h"
#include "ChangeParameterEventAction.h"
#include "FixedValueActionFunction.h"
#include "StateBasedTriggerEvent.h"
#include "StandardEventHandler.h"

using namespace std;

int main()
{

  typedef boost::numeric::ublas::vector<double> denseVector;
  typedef std::vector<denseVector> denseMatrix;
  typedef std::vector<std::vector<int> > graphType;

  denseVector X=TestModelInitialPopulations<denseVector>();
  
  denseMatrix NU=TestModelStoichiometry<denseMatrix>();
  
  TestModelPropensities<denseVector> PROPENSITIES;  
  graphType DG=TestModelDependencyGraph<graphType>();
  

  std::size_t numRuns=3;
  double t0=0;
  double tf=10;
  std::size_t intervals=10;

  typedef StatsAndTrajectoriesOutput<denseVector> outputType;

  outputType output;
  output.setKeepTrajectories(false);
  output.setOutputTimes(IntervalOutput<denseVector>::createUniformOutputTimes(t0,tf,intervals));


  boost::shared_ptr<SimpleMessageHandler> messageHandler(new SimpleMessageHandler());
  //messageHandler->setOutputOption(SimpleMessageHandler::DISPLAY_MESSAGES_ALL);

  typedef SSA_Direct_Events<denseVector, denseMatrix, TestModelPropensities<denseVector>, graphType, SimpleMessageHandler> solverType;
  typedef boost::function<void (double, denseVector&)> EventAction;
  solverType ssa(X,NU,PROPENSITIES,DG,messageHandler);

  StandardEventHandler<denseVector> eventHandler;

  //create two time-based triggers
  TimeBasedTrigger timeTrigger0(1.0);//trigger fires at t=1
  TimeBasedTrigger timeTrigger1(5.0);//trigger fires at t=5
  TimeBasedTrigger timeTrigger2(7.0);//trigger fires at t=7

  //create two (custom) state-based triggers
  CustomStateBasedTrigger<denseVector> stateTrigger0(&_customStateTrigger0<denseVector>);//trigger returns true if species 1 >=10, false otherwise
  CustomStateBasedTrigger<denseVector> stateTrigger1(&_customStateTrigger1<denseVector>);//trigger returns true if species 0<=9990, false otherwise

  //create an action that sets x0=10
  FixedValueActionFunction<denseVector> valFunc1(10.0);
  ChangeSingleSpeciesEventAction<solverType,FixedValueActionFunction<denseVector> > action1(0,valFunc1,ssa);//action sets species 0 to value of 10

  //create an action that sets x1=200
  ChangeSingleSpeciesEventAction<solverType,FixedValueActionFunction<denseVector> > action2(1,200.0,ssa);//action sets species 1 to value of 200

  //create an action that sets parameter index 0 to value of 50.0
  ChangeParameterEventAction<solverType,FixedValueActionFunction<denseVector> > action3(0,50.0,ssa);//action sets parameter 0 to value of 50

  //create an action that sets parameter index 1 to value of 10.0
  ChangeParameterEventAction<solverType,FixedValueActionFunction<denseVector> > action4(1,valFunc1,ssa);//action sets parameter 1 to value of 10

  //create a custom action that divides the population of x0 in half
  typedef double (*customActionFunction)(double,denseVector&);
  ChangeSingleSpeciesEventAction<solverType,customActionFunction> action5(0,&_customChangeSingleSpeciesFunction1,ssa);//action sets species 0 population to half of its current value

  //add events to eventHandler
  eventHandler.insertTimeEvent(new TimeBasedTriggerEvent<denseVector>(timeTrigger1,action1));//at t=5, set species 0 to 10

  //create a time event with multiple actions
  std::vector<EventAction> actions;
  actions.push_back(action1);
  actions.push_back(action3);


  StateBasedTriggerEvent<denseVector> testEvent(stateTrigger1,action4);

  eventHandler.insertTimeEvent(new TimeBasedTriggerEvent<denseVector>(timeTrigger2,actions));//at t=7, set species 0 to 10 and parameter 0 to 50

  eventHandler.insertStateEvent(new StateBasedTriggerEvent<denseVector>(stateTrigger1,action4));//when species 0 is <= 9990, set parameter 1 to 10

  actions.push_back(action2);
  eventHandler.insertStateEvent(new StateBasedTriggerEvent<denseVector>(stateTrigger0,actions,false));//when species 1>=10 (for the first time only), do actions 1,3 and 2

  eventHandler.insertTimeEvent(new TimeBasedTriggerEvent<denseVector>(timeTrigger0,action5));

  ssa.simulateEvents<outputType, StandardEventHandler<denseVector> >(numRuns,t0,tf,output,eventHandler);

  return 0;
}
