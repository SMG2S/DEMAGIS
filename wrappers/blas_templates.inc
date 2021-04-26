#pragma once

#include <assert.h>
#include <memory>

#include "blas_fortran.hpp"
#include "fortran_mangle.h"

template <>
void t_latms(const std::size_t m, const std::size_t n,
             const char dist, std::size_t* iseed, char sym, float *d, const std::size_t mode,
             const float cond, const float dmax, const std::size_t kl, const std::size_t ku,
             char pack, float* a, const std::size_t lda){

  using T = std::remove_reference<decltype((a[0]))>::type;
	     
  BlasInt m_ = m;
  BlasInt n_ = n;
  BlasInt iseed_[4];
  BlasInt mode_ = mode;
  BlasInt kl_ = kl;
  BlasInt ku_ = ku;
  BlasInt lda_ = lda;
  BlasInt info_;

  for(int i = 0; i < 4; i++){
    iseed_[i] = iseed[i];
  }

  T* work;

  int lwork = std::max(1,3*(std::max(n_,m_)));

  auto ptr = std::unique_ptr<T[]> {
    new T[ lwork ]
  };
  work = ptr.get();

  FC_GLOBAL(slatms, SLATMS)
  (&m_, &n_, &dist, iseed_, &sym, d, &mode_, &cond, &dmax, &kl_, &ku_, &pack, a, &lda_, work, &info_);  

}

template <>
void t_latms(const std::size_t m, const std::size_t n,
             const char dist, std::size_t* iseed, char sym, double *d, const std::size_t mode,
             const double cond, const double dmax, const std::size_t kl, const std::size_t ku,
             char pack, double* a, const std::size_t lda){

  using T = std::remove_reference<decltype((a[0]))>::type;

  BlasInt m_ = m;
  BlasInt n_ = n;
  BlasInt iseed_[4];
  BlasInt mode_ = mode;
  BlasInt kl_ = kl;
  BlasInt ku_ = ku;
  BlasInt lda_ = lda;
  BlasInt info_;

  for(int i = 0; i < 4; i++){
    iseed_[i] = iseed[i];
  }

  T* work;

  int lwork = std::max(1,3*(std::max(n_,m_)));

  auto ptr = std::unique_ptr<T[]> {
    new T[ lwork ]
  };
  work = ptr.get();

  FC_GLOBAL(dlatms, DLATMS)
  (&m_, &n_, &dist, iseed_, &sym, d, &mode_, &cond, &dmax, &kl_, &ku_, &pack, a, &lda_, work, &info_);

}


template <>
void t_latms(const std::size_t m, const std::size_t n,
             const char dist, std::size_t* iseed, char sym, float *d, const std::size_t mode,
             const float cond, float dmax, const std::size_t kl, const std::size_t ku,
             char pack, std::complex<float>* a, const std::size_t lda){

  using T = std::remove_reference<decltype((a[0]))>::type;

  BlasInt m_ = m;
  BlasInt n_ = n;
  BlasInt iseed_[4];
  BlasInt mode_ = mode;
  BlasInt kl_ = kl;
  BlasInt ku_ = ku;
  BlasInt lda_ = lda;
  BlasInt info_;

  for(int i = 0; i < 4; i++){
    iseed_[i] = iseed[i];
  }

  T* work;

  int lwork = std::max(1,3*(std::max(n_,m_)));

  auto ptr = std::unique_ptr<T[]> {
    new T[ lwork ]
  };
  work = ptr.get();

  FC_GLOBAL(clatms, CLATMS)
  (&m_, &n_, &dist, iseed_, &sym, d, &mode_, &cond, &dmax, &kl_, &ku_, &pack, a, &lda_, work, &info_);

}

template <>
void t_latms(const std::size_t m, const std::size_t n,
             const char dist, std::size_t* iseed, char sym, double *d, const std::size_t mode,
             const double cond, const double dmax, const std::size_t kl, const std::size_t ku,
             char pack, std::complex<double>* a, const std::size_t lda){

  using T = std::remove_reference<decltype((a[0]))>::type;

  BlasInt m_ = m;
  BlasInt n_ = n;
  BlasInt iseed_[4];
  BlasInt mode_ = mode;
  BlasInt kl_ = kl;
  BlasInt ku_ = ku;
  BlasInt lda_ = lda;
  BlasInt info_;

  for(int i = 0; i < 4; i++){
    iseed_[i] = iseed[i];
  }

  T* work;

  int lwork = std::max(1,3*(std::max(n_,m_)));

  auto ptr = std::unique_ptr<T[]> {
    new T[ lwork ]
  };
  work = ptr.get();

  FC_GLOBAL(zlatms, ZLATMS)
  (&m_, &n_, &dist, iseed_, &sym, d, &mode_, &cond, &dmax, &kl_, &ku_, &pack, a, &lda_, work, &info_);

}

