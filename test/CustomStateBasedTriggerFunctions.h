#ifndef _CUSTOMSTATEBASEDTRIGGERFUNCTIONS_H
#define _CUSTOMSTATEBASEDTRIGGERFUNCTIONS_H

template<typename _populationVectorType>
bool _customStateTrigger0(double t, _populationVectorType& x) {
  return (x[1]>=10.0);
}

template<typename _populationVectorType>
bool _customStateTrigger1(double t, _populationVectorType& x) {
  return (x[0]<=9990.0);
}



#endif
