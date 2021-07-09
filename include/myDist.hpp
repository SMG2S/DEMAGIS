#pragma once

#include "types.hpp"

template<typename T>
T* myUniformDist(std::size_t n, Base<T> epsilon, Base<T> dmax){
  Base<T>* eigenv = new Base<T>[n];
  for(auto k = 0; k < n; k++){
    eigenv[k] = dmax * (epsilon + (Base<T>)k * (1.0 - epsilon) / (Base<T>)n);
  }
  return eigenv;
}

template<typename T>
T* myGeometricDist(std::size_t n, Base<T> epsilon, Base<T> dmax){
  Base<T>* eigenv = new Base<T>[n];
  for(auto k = 0; k < n; k++){
    Base<T> exp = (Base<T>)(n - k) / (Base<T>)n;
    eigenv[k] = dmax * std::pow(epsilon, exp);
  }
  return eigenv;  
}

