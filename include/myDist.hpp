#pragma once

#include "../wrappers/types.hpp"

template<typename T>
void myUniformDist(Base<T>* &eigenv, std::size_t n, Base<T> epsilon, Base<T> dmax){
  for(auto k = 0; k < n; k++){
    eigenv[k] = dmax * (epsilon + (Base<T>)k * (1.0 - epsilon) / (Base<T>)n);
  }
}

template<typename T>
void myGeometricDist(Base<T>* &eigenv, std::size_t n, Base<T> epsilon, Base<T> dmax){
  for(auto k = 0; k < n; k++){
    Base<T> exp = (Base<T>)(n - k) / (Base<T>)n;
    eigenv[k] = dmax * std::pow(epsilon, exp);
  }
}
