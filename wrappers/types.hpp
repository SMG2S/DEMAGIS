#pragma once

template <class Q>
struct Base_Class {
  typedef Q type;
};

template <class Q>
struct Base_Class<std::complex<Q>> {
  typedef Q type;
};

template <typename Q>
using Base = typename Base_Class<Q>::type;

