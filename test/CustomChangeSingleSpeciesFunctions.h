#ifndef _CUSTOM_CHANGE_SINGLE_SPECIES_FUNCTIONS_H
#define _CUSTOM_CHANGE_SINGLE_SPECIES_FUNCTIONS_H

template<typename _populationValueType,
         typename _populationVectorType>
_populationValueType _customChangeSingleSpeciesFunction1(double t, _populationVectorType& x) {
  //std::cout << "original population of x[0] was " << x[0] << std::endl;
  return round(x[0]/2);
}

#endif
