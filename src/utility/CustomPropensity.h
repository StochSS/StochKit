#if !defined(__CustomPropensity_h__)
#define __CustomPropensity_h__

#include <iostream>
#include <vector>
#include "boost/numeric/ublas/vector.hpp"
#include "Parameter.h"

namespace STOCHKIT
{
 template<typename _populationVectorType>
 class CustomPropensity 
 {	
 public:
  
  typedef double (*customPropensityFunction)(_populationVectorType&, ParameterSet&);
  customPropensityFunction propensityFunction;

  CustomPropensity(){};

  CustomPropensity(customPropensityFunction func):
    propensityFunction(func)
    {}

  virtual double operator()(_populationVectorType& x, ParameterSet &ParametersList) {
    return (*propensityFunction)(x, ParametersList);
  }

  virtual ~CustomPropensity(){};

 };
}

#endif
