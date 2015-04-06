#ifndef _CONSTANT_TIME_GROUP_COLLECTION_H_
#define _CONSTANT_TIME_GROUP_COLLECTION_H_

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <utility>
#include <limits>
#include <deque>
#include <math.h>
#include "Random.h"
#include "StandardDriverTypes.h"
#include "ConstantTimeGroup.h"

/**
   @file ConstantTimeGroupCollection.h
   
   @brief Header file for ConstantTimeGroupCollection class used in constant time solver

   @warning This file is missing doxygen comments
*/

namespace STOCHKIT
{
  
  /**
     @class ConstantTimeGroupCollection
     
     @brief Helper class for ConstantTimeSSA
  */

  class ConstantTimeGroupCollection
 {
 public:

  ConstantTimeGroupCollection(std::size_t numberOfReactions);

  void build(boost::numeric::ublas::vector<double> &propensities);
  void update(std::size_t reactionIndex, double oldPropensity, double newPropensity);

  int selectReaction(STOCHKIT::RandomGenerator &randomGenerator);

  double getPropensitySum();

  int getGroup(double propensityValue);

  void addGroup(int groupExponent);

  void recalculatePropensitySum();

#ifdef DEBUG
  void printGroups();
#endif

 private:
  ConstantTimeGroupCollection();

  std::deque<ConstantTimeGroup> groups;

  double propensitySum;
  int maxGroupExponent;
  int minGroupExponent;

  std::vector<int> withinGroupIndexes;

#ifdef DEBUG
 public:
#endif
  int selectGroupIndex(STOCHKIT::RandomGenerator &randomGenerator);
 };
}

#endif

