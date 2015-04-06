/******************************************************************************
 *  FILE:    Random.cpp                                                       *
 ******************************************************************************/

#include "Random.h"

namespace STOCHKIT
{
#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
#ifndef BOOST_RANDOM_NO_STREAM_OPERATORS
  // output random generator state to an ostream
  template<class CharT, class Traits>
  std::basic_ostream<CharT,Traits>&
  operator<<(std::basic_ostream<CharT,Traits>& os, const RandomGenerator& randomGenerator)
  {
    os << randomGenerator.generatorUniform;
    return os;
  }

  // input random generator state from an istream
  // WARNING: no check to make sure this works, make sure the istream only contains relevant information, use at your own risk
  template<class CharT, class Traits>
  std::basic_istream<CharT,Traits>&
  operator>>(std::basic_istream<CharT,Traits>& is, RandomGenerator& randomGenerator)
  {
    is >> randomGenerator.generatorUniform;
    return is;
  }
#endif
#endif

}//end namespace
