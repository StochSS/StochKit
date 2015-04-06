/*
 *  test_SSA_Direct.cpp
 *  
 */
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <iomanip>
#include "ConstantTimeGroup.h"
#include "Random.h"
#include "ConstantTimeGroupCollection.h"
#include "StandardDriverTypes.h"
#include "SSA_ConstantTime.h"
#include "Input_mass_action.h"

using namespace std;


int main()
{
  STOCHKIT::RandomGenerator randomGenerator;

  /*
  double myPropensityValue=6.5;
  int exponent=ConstantTimeGroup::calculateGroupExponent(myPropensityValue);

  ConstantTimeGroup group1(exponent);

  std::vector<double> propensities;
  propensities.push_back(6.5);
  propensities.push_back(4.01);
  propensities.push_back(8);
  propensities.push_back(7.3);

  std::vector<int> wGI(4,-1);//wGI=withinGroupIndexes

  group1.insert(propensities[0], 0, wGI);

  std::cout << "group1.getGroupSum()="<<group1.getGroupSum()<<"\n";
  std::cout << "group1.getGroupExponent()="<<group1.getGroupExponent()<<"\n";
  
  STOCHKIT::RandomGenerator randomGenerator;

  for (std::size_t i=0; i!=20; ++i) {
    std::cout << "selectReaction returned: "<<group1.selectReactionIndex(randomGenerator)<<"\n";
  }

  group1.insert(propensities[1], 1, wGI);

  std::vector<std::size_t> counts(2,0);
  for (std::size_t i=0; i!=10510000; ++i) {
    ++counts[group1.selectReactionIndex(randomGenerator)];
  }
  for (std::size_t i=0; i!=counts.size(); ++i) {
    std::cout << "selected rxn "<<i<<": "<<counts[i]<<"\n";
  }

  group1.insert(propensities[2], 2, wGI);

  std::vector<std::size_t> counts3(3,0);
  for (std::size_t i=0; i!=18510000; ++i) {
    ++counts3[group1.selectReactionIndex(randomGenerator)];
  }
  for (std::size_t i=0; i!=counts3.size(); ++i) {
    std::cout << "selected rxn "<<i<<": "<<counts3[i]<<"\n";
  }

  std::cout << "group1.getGroupSum()="<<group1.getGroupSum()<<"\n";

  std::cout << "wGI[0]="<<wGI[0]<<"\n";
  std::cout << "wGI[1]="<<wGI[1]<<"\n";
  std::cout << "wGI[2]="<<wGI[2]<<"\n";
  std::cout << "wGI[3]="<<wGI[3]<<"\n";

  std::cout << "removing reaction 0...\n";
  group1.remove(0,wGI[0],wGI);

  std::cout << "group1.getGroupSum()="<<group1.getGroupSum()<<"\n";

  std::cout << "wGI[0]="<<wGI[0]<<"\n";
  std::cout << "wGI[1]="<<wGI[1]<<"\n";
  std::cout << "wGI[2]="<<wGI[2]<<"\n";
  std::cout << "wGI[3]="<<wGI[3]<<"\n";

  std::cout << "removing reaction 1...\n";
  group1.remove(1,wGI[1],wGI);

  std::cout << "group1.getGroupSum()="<<group1.getGroupSum()<<"\n";

  std::cout << "wGI[0]="<<wGI[0]<<"\n";
  std::cout << "wGI[1]="<<wGI[1]<<"\n";
  std::cout << "wGI[2]="<<wGI[2]<<"\n";
  std::cout << "wGI[3]="<<wGI[3]<<"\n";

  std::cout << "inserting index 0...\n";
  group1.insert(propensities[0], 0, wGI);  

  std::cout << "group1.getGroupSum()="<<group1.getGroupSum()<<"\n";

  std::cout << "wGI[0]="<<wGI[0]<<"\n";
  std::cout << "wGI[1]="<<wGI[1]<<"\n";
  std::cout << "wGI[2]="<<wGI[2]<<"\n";
  std::cout << "wGI[3]="<<wGI[3]<<"\n";

  //std::cout << "exponent if passed 0: "<<ConstantTimeGroup::calculateGroupExponent(0.0)<<"\n";//error

  std::cout << "creating groups...\n";
  ConstantTimeGroupCollection groups(propensities.size());
  std::cout << "building...\n";
  groups.build(propensities);

  groups.printGroups();

  std::vector<int> counts4(5,0);
  for (std::size_t i=0; i!=25810000; ++i) {
    int reaction=groups.selectReaction(randomGenerator);
    ++counts4[reaction+1];//add 1 so index 0 corresponds to returning rxn -1
  }
  std::cout << "selected rxn -1: "<<counts4[0]<<"\n";

  for (std::size_t i=1; i!=counts4.size(); ++i) {
    std::cout << "selected rxn "<<i-1<<": "<<counts4[i]<<"\n";
  }
  */

  /*
  std::vector<double> propensities2(6,0.0);
  propensities2[1]=23;
  propensities2[3]=66;
  propensities2[4]=257;
  propensities2[5]=67;
  
  ConstantTimeGroupCollection groups2(propensities2.size());
  std::cout << "building...\n";
  groups2.build(propensities2);
  
  groups2.printGroups();
  

  std::vector<int> groupCounts(7,0);
  std::size_t num=(std::size_t)(groups2.getPropensitySum()*10000.0);
  std::cout<<"num="<<num<<"\n";
  for (std::size_t i=0; i!=num; ++i) {
    //    for (std::size_t i=0; i!=5; ++i) {
    int reaction=groups2.selectGroupIndex(randomGenerator);
    ++groupCounts[reaction+1];//add 1 so index 0 corresponds to returning rxn -1
  }
  std::cout << "selected group -1: "<<groupCounts[0]<<"\n";
  for (std::size_t i=1; i!=groupCounts.size(); ++i) {
    std::cout << "selected group "<<i-1<<": "<<groupCounts[i]<<"\n";
  }

  std::vector<int> counts6(7,0);
  for (std::size_t i=0; i!=(std::size_t)(groups2.getPropensitySum()*10000); ++i) {
    //for (std::size_t i=0; i!=5; ++i) {
    int reaction=groups2.selectReaction(randomGenerator);
    ++counts6[reaction+1];//add 1 so index 0 corresponds to returning rxn -1
  }
  std::cout << "selected rxn -1: "<<counts6[0]<<"\n";
  for (std::size_t i=1; i!=counts6.size(); ++i) {
    std::cout << "selected rxn "<<i-1<<": "<<counts6[i]<<"\n";
  }


  std::vector<double> newPropensities2(6,0.0);
  newPropensities2[0]=16;
  newPropensities2[1]=24;
  newPropensities2[2]=1025;
  newPropensities2[3]=800;
  newPropensities2[4]=0;
  newPropensities2[5]=67;

  std::cout << "updating propensities...\n";
  for (std::size_t i=0; i!=newPropensities2.size(); ++i) {
    std::cout << "updating propensity "<<i<<"\n";
    groups2.update(i,propensities2[i],newPropensities2[i]);
  }
  groups2.printGroups();

  std::vector<int> counts62(7,0);
  for (std::size_t i=0; i!=(std::size_t)groups2.getPropensitySum()*10000; ++i) {
    //for (std::size_t i=0; i!=5; ++i) {
    int reaction=groups2.selectReaction(randomGenerator);
    ++counts62[reaction+1];//add 1 so index 0 corresponds to returning rxn -1
  }
  std::cout << "selected rxn -1: "<<counts62[0]<<"\n";
  for (std::size_t i=1; i!=counts62.size(); ++i) {
    std::cout << "selected rxn "<<i-1<<": "<<counts62[i]<<"\n";
  }
  */

  //try running a simulation...

  typedef SSA_ConstantTime<StandardDriverTypes::populationType,
    StandardDriverTypes::stoichiometryType, 
    StandardDriverTypes::propensitiesType,
    StandardDriverTypes::graphType> _solverType;

  typedef _solverType::populationVectorType populationVectorType;
  typedef _solverType::stoichiometryType stoichiometryType;
  typedef _solverType::propensitiesType propensitiesType;
  typedef _solverType::dependencyGraphType dependencyGraphType;
  typedef StandardDriverTypes::outputType outputType;  

  std::string modelFile="models/examples/stochkit_mass_action.xml";
  Input_mass_action<populationVectorType, stoichiometryType, propensitiesType, dependencyGraphType> model(const_cast<char*>(modelFile.c_str()));

  _solverType solver(model.writeInitialPopulation(),
		     model.writeStoichiometry(),
		     model.writePropensities(),
		     model.writeDependencyGraph());

  std::cout << "created solver...\n";

  solver.initialize();

  std::cout << "taking step...\n";
  solver.step();
  std::cout << "after step, t="<<solver.currentTime<<"\n";
  std::cout << "currentPop:\n";
  for (std::size_t i=0; i!=solver.currentPopulation.size(); ++i) {
    std::cout << "["<<i<<"]="<<solver.currentPopulation[i]<<"\n";
  }
  solver.groups.printGroups();

  std::cout << "taking step...\n";
  solver.step();
  std::cout << "after step, t="<<solver.currentTime<<"\n";
  std::cout << "currentPop:\n";
  for (std::size_t i=0; i!=solver.currentPopulation.size(); ++i) {
    std::cout << "["<<i<<"]="<<solver.currentPopulation[i]<<"\n";
  }
  solver.groups.printGroups();

  std::cout << "taking step...\n";
  solver.step();
  std::cout << "after step, t="<<solver.currentTime<<"\n";
  std::cout << "currentPop:\n";
  for (std::size_t i=0; i!=solver.currentPopulation.size(); ++i) {
    std::cout << "["<<i<<"]="<<solver.currentPopulation[i]<<"\n";
  }
  solver.groups.printGroups();

  while (solver.step()) {
    for (std::size_t i=0; i!=solver.currentPopulation.size(); ++i) {
      std::cout << "["<<i<<"]="<<solver.currentPopulation[i]<<" ";
    }
    std::cout<<"\n";
  }
  std::cout << "propensity should be zero, groups should be empty:\n";
  solver.groups.printGroups();

  return 0;
}

