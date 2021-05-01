#pragma once

#include <assert.h>
#include <memory>

#include "blas_fortran.hpp"
#include "fortran_mangle.h"

template <>
std::size_t t_geqrf(std::size_t m, std::size_t n, double* a,
        std::size_t lda, double* tau)
{
//return LAPACKE_dgeqrf(matrix_layout, m, n, a, lda, tau);
  using T = std::remove_reference<decltype((a[0]))>::type;

  BlasInt m_ = m;
  BlasInt n_ = n;
  BlasInt lda_ = lda;

  T* work;
  T numwork;
  BlasInt lwork, info;

  lwork = -1;
  FC_GLOBAL(dgeqrf,DGEQRF)(&m_, &n_, a, &lda_, tau, &numwork, &lwork, &info);
  assert(info == 0);

  lwork = static_cast<std::size_t>((numwork));
  auto ptr = std::unique_ptr<T[]> {
    new T[ lwork ]
  };
  work = ptr.get();

  FC_GLOBAL(dgeqrf,DGEQRF)(&m_, &n_, a, &lda_, tau, work, &lwork, &info);
  assert(info == 0);
  return info;
}
template <>
std::size_t t_geqrf(std::size_t m, std::size_t n, float* a,
        std::size_t lda, float* tau)
{
  //return LAPACKE_sgeqrf(matrix_layout, m, n, a, lda, tau);
  using T = std::remove_reference<decltype((a[0]))>::type;

  BlasInt m_ = m;
  BlasInt n_ = n;
  BlasInt lda_ = lda;

  T* work;
  T numwork;
  BlasInt lwork, info;

  lwork = -1;
  FC_GLOBAL(sgeqrf,SGEQRF)(&m_, &n_, a, &lda_, tau, &numwork, &lwork, &info);
  assert(info == 0);

  lwork = static_cast<std::size_t>((numwork));
  auto ptr = std::unique_ptr<T[]> {
    new T[ lwork ]
  };
  work = ptr.get();

  FC_GLOBAL(sgeqrf,SGEQRF)(&m_, &n_, a, &lda_, tau, work, &lwork, &info);
  assert(info == 0);
  return info;
}

template <>
std::size_t t_geqrf(std::size_t m, std::size_t n,
                    std::complex<double>* a, std::size_t lda,
                    std::complex<double>* tau) {
  // return LAPACKE_zgeqrf(matrix_layout, m, n, a, lda, tau);

  using T = std::remove_reference<decltype((a[0]))>::type;

  BlasInt m_ = m;
  BlasInt n_ = n;
  BlasInt lda_ = lda;

  T* work;
  T numwork;
  BlasInt lwork, info;

  lwork = -1;
  FC_GLOBAL(zgeqrf, ZGEQRF)(&m_, &n_, a, &lda_, tau, &numwork, &lwork, &info);
  assert(info == 0);

  lwork = static_cast<std::size_t>(real(numwork));
  auto ptr = std::unique_ptr<T[]>{new T[lwork]};
  work = ptr.get();

  FC_GLOBAL(zgeqrf, ZGEQRF)(&m_, &n_, a, &lda_, tau, work, &lwork, &info);
  assert(info == 0);
  return info;
}

template <>
std::size_t t_geqrf(std::size_t m, std::size_t n,
                    std::complex<float>* a, std::size_t lda,
                    std::complex<float>* tau) {
  using T = std::remove_reference<decltype((a[0]))>::type;

  //  return LAPACKE_cgeqrf(matrix_layout, m, n, a, lda, tau);
  BlasInt m_ = m;
  BlasInt n_ = n;
  BlasInt lda_ = lda;

  T* work;
  T numwork;
  BlasInt lwork, info;

  lwork = -1;
  FC_GLOBAL(cgeqrf, CGEQRF)(&m_, &n_, a, &lda_, tau, &numwork, &lwork, &info);
  assert(info == 0);

  lwork = static_cast<std::size_t>(real(numwork));
  auto ptr = std::unique_ptr<T[]>{new T[lwork]};
  work = ptr.get();

  FC_GLOBAL(cgeqrf, CGEQRF)(&m_, &n_, a, &lda_, tau, work, &lwork, &info);
  assert(info == 0);
  return info;
}

template <>
std::size_t t_mqr(char side, char trans, std::size_t m, std::size_t n,  std::size_t k, float* a, std::size_t lda,
                  float *tau, float *c, std::size_t ldc){
  using T = std::remove_reference<decltype((c[0]))>::type;

  BlasInt m_ = m;
  BlasInt n_ = n;
  BlasInt k_ = k;  
  BlasInt lda_ = lda;  
  BlasInt ldc_ = ldc;  

  T* work;
  T numwork;
  BlasInt lwork, info;

  lwork = -1;
  FC_GLOBAL(sormqr, SORMQR)(&side, &trans, &m_, &n_, &k_, a, &lda_, tau, c, &ldc_, &numwork, &lwork, &info);  
  assert(info == 0);

  lwork = static_cast<std::size_t>((numwork));
  auto ptr = std::unique_ptr<T[]>{new T[lwork]};
  work = ptr.get();

  FC_GLOBAL(sormqr, SORMQR)(&side, &trans, &m_, &n_, &k_, a, &lda_, tau, c, &ldc_, work, &lwork, &info);  
  assert(info == 0);
  return info;
}

template <>
std::size_t t_mqr(char side, char trans, std::size_t m, std::size_t n,  std::size_t k, double* a, std::size_t lda,
                  double *tau, double *c, std::size_t ldc){
  using T = std::remove_reference<decltype((c[0]))>::type;

  BlasInt m_ = m;
  BlasInt n_ = n;
  BlasInt k_ = k;
  BlasInt lda_ = lda;
  BlasInt ldc_ = ldc;

  T* work;
  T numwork;
  BlasInt lwork, info;

  lwork = -1;
  FC_GLOBAL(dormqr, DORMQR)(&side, &trans, &m_, &n_, &k_, a, &lda_, tau, c, &ldc_, &numwork, &lwork, &info);
  assert(info == 0);

  lwork = static_cast<std::size_t>((numwork));
  auto ptr = std::unique_ptr<T[]>{new T[lwork]};
  work = ptr.get();

  FC_GLOBAL(dormqr, DORMQR)(&side, &trans, &m_, &n_, &k_, a, &lda_, tau, c, &ldc_, work, &lwork, &info);
  assert(info == 0);
  return info;
}

template <>
std::size_t t_mqr(char side, char trans, std::size_t m, std::size_t n,  std::size_t k, std::complex<float>* a, std::size_t lda,
                  std::complex<float> *tau, std::complex<float> *c, std::size_t ldc){
  using T = std::remove_reference<decltype((c[0]))>::type;

  BlasInt m_ = m;
  BlasInt n_ = n;
  BlasInt k_ = k;
  BlasInt lda_ = lda;
  BlasInt ldc_ = ldc;

  T* work;
  T numwork;
  BlasInt lwork, info;

  lwork = -1;
  FC_GLOBAL(cunmqr, CUNMQR)(&side, &trans, &m_, &n_, &k_, a, &lda_, tau, c, &ldc_, &numwork, &lwork, &info);
  assert(info == 0);

  lwork = static_cast<std::size_t>(real(numwork));
  auto ptr = std::unique_ptr<T[]>{new T[lwork]};
  work = ptr.get();

  FC_GLOBAL(cunmqr, CUNMQR)(&side, &trans, &m_, &n_, &k_, a, &lda_, tau, c, &ldc_, work, &lwork, &info);
  assert(info == 0);
  return info;
}

template <>
std::size_t t_mqr(char side, char trans, std::size_t m, std::size_t n,  std::size_t k, std::complex<double>* a, std::size_t lda,
                  std::complex<double> *tau, std::complex<double> *c, std::size_t ldc){
  using T = std::remove_reference<decltype((c[0]))>::type;

  BlasInt m_ = m;
  BlasInt n_ = n;
  BlasInt k_ = k;
  BlasInt lda_ = lda;
  BlasInt ldc_ = ldc;

  T* work;
  T numwork;
  BlasInt lwork, info;

  lwork = -1;
  FC_GLOBAL(zunmqr, ZUNMQR)(&side, &trans, &m_, &n_, &k_, a, &lda_, tau, c, &ldc_, &numwork, &lwork, &info);
  assert(info == 0);

  lwork = static_cast<std::size_t>(real(numwork));
  auto ptr = std::unique_ptr<T[]>{new T[lwork]};
  work = ptr.get();

  FC_GLOBAL(zunmqr, ZUNMQR)(&side, &trans, &m_, &n_, &k_, a, &lda_, tau, c, &ldc_, work, &lwork, &info);
  assert(info == 0);
  return info;
}


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

////////////
// BLACS //
////////////
void blacs_pinfo(int *mypnum, int *nprocs){
  BlasInt mypnum_, nprocs_;
  FC_GLOBAL(blacs_pinfo, BLACS_PINFO)(&mypnum_, &nprocs_);
  *mypnum = mypnum_;
  *nprocs = nprocs_;
}

void blacs_get(int *icontxt, const int *what, int *val ){
  BlasInt val_;
  BlasInt icontxt_;
  BlasInt what_ = *what;
  FC_GLOBAL(blacs_get, BLACS_GET)(&icontxt_, &what_, &val_);
  *val = val_; 
  *icontxt = icontxt_; 
}

void blacs_gridinit(int *icontxt, const char layout, const int *nprow, const int *npcol){
  BlasInt icontxt_;
  FC_GLOBAL(blacs_gridinit, BLACS_GRIDINIT)(&icontxt_, &layout, nprow, npcol);  
  *icontxt = icontxt_; 
}

void blacs_gridinfo(int *icontxt, int *nprow, int *npcol, int *myprow, int *mypcol){
  BlasInt nprow_, npcol_, myprow_, mypcol_, icontxt_;
  icontxt_ = *icontxt;
  FC_GLOBAL(blacs_gridinfo, BLACS_GRIDINFO)( &icontxt_, &nprow_, &npcol_, &myprow_, &mypcol_);
  *nprow = nprow_;
  *npcol = npcol_;
  *myprow = myprow_;
  *mypcol = mypcol_;
}

std::size_t numroc(std::size_t *n, std::size_t *nb, int *iproc, const int *isrcproc, int *nprocs){
  std::size_t val;
  BlasInt val_, n_, nb_;
  n_ = static_cast<int>(*n);
  nb_ = static_cast<int>(*nb);
  val_ = FC_GLOBAL(numroc, NUMROC)(&n_, &nb_, iproc, isrcproc, nprocs);
  val = static_cast<std::size_t>(val_);
  return val;
}

void descinit(std::size_t *desc, std::size_t *m, std::size_t *n, std::size_t *mb, std::size_t *nb, 
	const int *irsrc, const int *icsrc, int *ictxt, std::size_t *lld, int *info)
{
  BlasInt *desc_ = new BlasInt[9];
  BlasInt m_ = *m;
  BlasInt n_ = *n;
  BlasInt mb_ = *mb;
  BlasInt nb_ = *nb;
  const BlasInt irsrc_ = *irsrc;
  const BlasInt icsrc_ = *icsrc;
  BlasInt ictxt_ = *ictxt;
  BlasInt lld_ = *lld;
  BlasInt info_ = *info;
  FC_GLOBAL(descinit, DESCINIT)(desc_, &m_, &n_, &mb_, &nb_, &irsrc_, &icsrc_, &ictxt_, &lld_, &info_);
  *info = info_;

  for(int i = 0; i < 9; i++){
      desc[i] = static_cast<std::size_t>(desc_[i]);
  }

}



template <>
void t_pgeqrf(const std::size_t m, const std::size_t n, float *a, const std::size_t ia, const std::size_t ja,
              const std::size_t *desc_a, float *tau)
{
  using T = std::remove_reference<decltype((a[0]))>::type;

  BlasInt m_ = m;
  BlasInt n_ = n;  
  BlasInt ia_ = ia;
  BlasInt ja_ = ja;
  BlasInt *desc_a_ = new BlasInt[9];
  BlasInt info_ = 0;
  for(int i = 0; i < 9; i++){
      desc_a_[i] = static_cast<BlasInt>(desc_a[i]);
  }

  T* work; T numwork; BlasInt lwork;
  lwork = -1;
  FC_GLOBAL(psgeqrf, PSGEQRF)(&m_, &n_, a, &ia_, &ja_, desc_a_, tau, &numwork, &lwork, &info_);

  lwork = (BlasInt)numwork;

  auto wptr = std::unique_ptr<T[]> {
      new T[ lwork ]
  };

  work = wptr.get();
  FC_GLOBAL(psgeqrf, PSGEQRF)(&m_, &n_, a, &ia_, &ja_, desc_a_, tau, work, &lwork, &info_);
  assert(info_ == 0);

}

template <>
void t_pgeqrf(const std::size_t m, const std::size_t n, double *a, const std::size_t ia, const std::size_t ja,
              const std::size_t *desc_a, double *tau)
{
  using T = std::remove_reference<decltype((a[0]))>::type;

  BlasInt m_ = m;
  BlasInt n_ = n;
  BlasInt ia_ = ia;
  BlasInt ja_ = ja;
  BlasInt *desc_a_ = new BlasInt[9];
  BlasInt info_ = 0;
  for(int i = 0; i < 9; i++){
      desc_a_[i] = static_cast<BlasInt>(desc_a[i]);
  }

  T* work; T numwork; BlasInt lwork;
  lwork = -1;
  FC_GLOBAL(pdgeqrf, PDGEQRF)(&m_, &n_, a, &ia_, &ja_, desc_a_, tau, &numwork, &lwork, &info_);

  lwork = (BlasInt)numwork;

  auto wptr = std::unique_ptr<T[]> {
      new T[ lwork ]
  };

  work = wptr.get();
  FC_GLOBAL(pdgeqrf, PDGEQRF)(&m_, &n_, a, &ia_, &ja_, desc_a_, tau, work, &lwork, &info_);
  assert(info_ == 0);

}

template <>
void t_pgeqrf(const std::size_t m, const std::size_t n, std::complex<float> *a, const std::size_t ia, const std::size_t ja,
              const std::size_t *desc_a, std::complex<float> *tau)
{
  using T = std::remove_reference<decltype((a[0]))>::type;

  BlasInt m_ = m;
  BlasInt n_ = n;
  BlasInt ia_ = ia;
  BlasInt ja_ = ja;
  BlasInt *desc_a_ = new BlasInt[9];
  BlasInt info_ = 0;
  for(int i = 0; i < 9; i++){
      desc_a_[i] = static_cast<BlasInt>(desc_a[i]);
  }

  T* work; T numwork; BlasInt lwork;
  lwork = -1;
  FC_GLOBAL(pcgeqrf, PCGEQRF)(&m_, &n_, a, &ia_, &ja_, desc_a_, tau, &numwork, &lwork, &info_);

  lwork = static_cast<BlasInt>(real(numwork));

  auto wptr = std::unique_ptr<T[]> {
      new T[ lwork ]
  };

  work = wptr.get();
  FC_GLOBAL(pcgeqrf, PCGEQRF)(&m_, &n_, a, &ia_, &ja_, desc_a_, tau, work, &lwork, &info_);
  assert(info_ == 0);

}

template <>
void t_pgeqrf(const std::size_t m, const std::size_t n, std::complex<double> *a, const std::size_t ia, const std::size_t ja,
              const std::size_t *desc_a, std::complex<double> *tau)
{
  using T = std::remove_reference<decltype((a[0]))>::type;

  BlasInt m_ = m;
  BlasInt n_ = n;
  BlasInt ia_ = ia;
  BlasInt ja_ = ja;
  BlasInt *desc_a_ = new BlasInt[9];
  BlasInt info_ = 0;
  for(int i = 0; i < 9; i++){
      desc_a_[i] = static_cast<BlasInt>(desc_a[i]);
  }

  T* work; T numwork; BlasInt lwork;
  lwork = -1;
  FC_GLOBAL(pzgeqrf, PZGEQRF)(&m_, &n_, a, &ia_, &ja_, desc_a_, tau, &numwork, &lwork, &info_);

  lwork = static_cast<BlasInt>(real(numwork));

  auto wptr = std::unique_ptr<T[]> {
      new T[ lwork ]
  };

  work = wptr.get();
  FC_GLOBAL(pzgeqrf, PZGEQRF)(&m_, &n_, a, &ia_, &ja_, desc_a_, tau, work, &lwork, &info_);
  assert(info_ == 0);

}

template <>
void t_pmqr(const char side, const char trans, const std::size_t m, const std::size_t n, const std::size_t k,
            const float *a, const std::size_t ia, const std::size_t ja, const std::size_t *desc_a, const float *tau,
            float *c, const std::size_t ic, const std::size_t jc, const std::size_t *desc_c)
{
  using T = std::remove_reference<decltype((c[0]))>::type;

  BlasInt m_ = m;
  BlasInt n_ = n;
  BlasInt k_ = k;
  BlasInt ia_ = ia;
  BlasInt ja_ = ja;
  BlasInt ic_ = ic;
  BlasInt jc_ = jc;  

  BlasInt *desc_a_ = new BlasInt[9];
  BlasInt *desc_c_ = new BlasInt[9];

  BlasInt info_ = 0;
  for(int i = 0; i < 9; i++){
      desc_a_[i] = static_cast<BlasInt>(desc_a[i]);
      desc_c_[i] = static_cast<BlasInt>(desc_c[i]);
  }  

  T* work; T numwork; BlasInt lwork;
  lwork = -1;

  FC_GLOBAL(psormqr, PSORMQR)(&side, &trans, &m_, &n_,&k_, a, &ia_, &ja_, desc_a_, tau, c, &ic_, &jc_, desc_c_, &numwork, &lwork, &info_);

  lwork = (BlasInt)numwork;

  auto wptr = std::unique_ptr<T[]> {
      new T[ lwork ]
  };

  work = wptr.get();
  FC_GLOBAL(psormqr, PSORMQR)(&side, &trans, &m_, &n_,&k_, a, &ia_, &ja_, desc_a_, tau, c, &ic_, &jc_, desc_c_, work, &lwork, &info_);
  assert(info_ == 0);

}

template <>
void t_pmqr(const char side, const char trans, const std::size_t m, const std::size_t n, const std::size_t k,
            const double *a, const std::size_t ia, const std::size_t ja, const std::size_t *desc_a, const double *tau,
            double *c, const std::size_t ic, const std::size_t jc, const std::size_t *desc_c)
{
  using T = std::remove_reference<decltype((c[0]))>::type;

  BlasInt m_ = m;
  BlasInt n_ = n;
  BlasInt k_ = k;
  BlasInt ia_ = ia;
  BlasInt ja_ = ja;
  BlasInt ic_ = ic;
  BlasInt jc_ = jc;

  BlasInt *desc_a_ = new BlasInt[9];
  BlasInt *desc_c_ = new BlasInt[9];

  BlasInt info_ = 0;
  for(int i = 0; i < 9; i++){
      desc_a_[i] = static_cast<BlasInt>(desc_a[i]);
      desc_c_[i] = static_cast<BlasInt>(desc_c[i]);
  }

  T* work; T numwork; BlasInt lwork;
  lwork = -1;

  FC_GLOBAL(pdormqr, PDORMQR)(&side, &trans, &m_, &n_,&k_, a, &ia_, &ja_, desc_a_, tau, c, &ic_, &jc_, desc_c_, &numwork, &lwork, &info_);

  lwork = (BlasInt)numwork;

  auto wptr = std::unique_ptr<T[]> {
      new T[ lwork ]
  };

  work = wptr.get();
  FC_GLOBAL(pdormqr, PDORMQR)(&side, &trans, &m_, &n_,&k_, a, &ia_, &ja_, desc_a_, tau, c, &ic_, &jc_, desc_c_, work, &lwork, &info_);
  assert(info_ == 0);

}


template <>
void t_pmqr(const char side, const char trans, const std::size_t m, const std::size_t n, const std::size_t k,
            const std::complex<float> *a, const std::size_t ia, const std::size_t ja, const std::size_t *desc_a, const std::complex<float> *tau,
            std::complex<float> *c, const std::size_t ic, const std::size_t jc, const std::size_t *desc_c)
{
  using T = std::remove_reference<decltype((c[0]))>::type;

  BlasInt m_ = m;
  BlasInt n_ = n;
  BlasInt k_ = k;
  BlasInt ia_ = ia;
  BlasInt ja_ = ja;
  BlasInt ic_ = ic;
  BlasInt jc_ = jc;

  BlasInt *desc_a_ = new BlasInt[9];
  BlasInt *desc_c_ = new BlasInt[9];

  BlasInt info_ = 0;
  for(int i = 0; i < 9; i++){
      desc_a_[i] = static_cast<BlasInt>(desc_a[i]);
      desc_c_[i] = static_cast<BlasInt>(desc_c[i]);
  }

  T* work; T numwork; BlasInt lwork;
  lwork = -1;

  FC_GLOBAL(pcunmqr, PCUNMQR)(&side, &trans, &m_, &n_,&k_, a, &ia_, &ja_, desc_a_, tau, c, &ic_, &jc_, desc_c_, &numwork, &lwork, &info_);

  lwork = static_cast<BlasInt>(real(numwork));

  auto wptr = std::unique_ptr<T[]> {
      new T[ lwork ]
  };

  work = wptr.get();
  FC_GLOBAL(pcunmqr, PCUNMQR)(&side, &trans, &m_, &n_,&k_, a, &ia_, &ja_, desc_a_, tau, c, &ic_, &jc_, desc_c_, work, &lwork, &info_);
  assert(info_ == 0);

}


template <>
void t_pmqr(const char side, const char trans, const std::size_t m, const std::size_t n, const std::size_t k,
            const std::complex<double> *a, const std::size_t ia, const std::size_t ja, const std::size_t *desc_a, const std::complex<double> *tau,
            std::complex<double> *c, const std::size_t ic, const std::size_t jc, const std::size_t *desc_c)
{
  using T = std::remove_reference<decltype((c[0]))>::type;

  BlasInt m_ = m;
  BlasInt n_ = n;
  BlasInt k_ = k;
  BlasInt ia_ = ia;
  BlasInt ja_ = ja;
  BlasInt ic_ = ic;
  BlasInt jc_ = jc;

  BlasInt *desc_a_ = new BlasInt[9];
  BlasInt *desc_c_ = new BlasInt[9];

  BlasInt info_ = 0;
  for(int i = 0; i < 9; i++){
      desc_a_[i] = static_cast<BlasInt>(desc_a[i]);
      desc_c_[i] = static_cast<BlasInt>(desc_c[i]);
  }

  T* work; T numwork; BlasInt lwork;
  lwork = -1;

  FC_GLOBAL(pzunmqr, PZUNMQR)(&side, &trans, &m_, &n_,&k_, a, &ia_, &ja_, desc_a_, tau, c, &ic_, &jc_, desc_c_, &numwork, &lwork, &info_);

  lwork = static_cast<BlasInt>(real(numwork));

  auto wptr = std::unique_ptr<T[]> {
      new T[ lwork ]
  };

  work = wptr.get();
  FC_GLOBAL(pzunmqr, PZUNMQR)(&side, &trans, &m_, &n_,&k_, a, &ia_, &ja_, desc_a_, tau, c, &ic_, &jc_, desc_c_, work, &lwork, &info_);
  assert(info_ == 0);

}

