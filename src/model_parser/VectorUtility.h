#ifndef _VECTOR_UTILITY_H_
#define _VECTOR_UTILITY_H_

#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/type_traits.hpp>

namespace STOCHKIT
{

class VectorUtility
{
private:
   template <class T, class Z, int N>
   class helper
   {
   public:
      inline void func(T &row) const
      {
         assert(false);
         exit(1);
      }
   };

   template <class T, class Z>
   class helper<T, Z, 1>
   {
   public:
      inline void func(T &row) const
      {
         row = Z(row.size() );
      }
   };

   template <class T, class Z>
   class helper<T, Z, 0>
   {
   public:
      inline void func(T &row) const
      {
         for(unsigned int j = 0; j < row.size(); ++j)
         {
            row[j] = 0;
         }
      }
   };

public:
   template <class T>
   static inline void initializeVectorToZero(T &input)
   {
      typedef boost::numeric::ublas::zero_vector<double> zero_double_ublas_vec;
      const int N = boost::is_convertible<zero_double_ublas_vec, typename T::value_type>::value;
      helper<typename T::value_type, zero_double_ublas_vec, N> obj;
      for(unsigned int i=0; i < input.size(); ++i)
      {
         obj.func(input[i]);
      }
   }
};

}
#endif
