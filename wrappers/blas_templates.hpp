#pragma once

#include <complex>
#include "./types.hpp"

template <typename T>
std::size_t t_geqrf(std::size_t m, std::size_t n, T* a, std::size_t lda, T* tau);

template <typename T>
std::size_t t_mqr(const char side, const char trans, std::size_t m, std::size_t n,  std::size_t k, T* a, std::size_t lda,
		  T *tau, T *c, std::size_t ldc);

template <typename T>
void t_latms(const std::size_t m, const std::size_t n,
	     const char dist, std::size_t* iseed, char sym, Base<T> *d, const std::size_t mode,
	     const Base<T> cond, const Base<T> dmax, const std::size_t kl, const std::size_t ku,
	     char pack, T* a, const std::size_t lda);

////////////
// BLACS ///
////////////
void blacs_pinfo(int *mypnum, int *nprocs);
void blacs_get(int *icontxt, const int *what, int *val );
void blacs_gridinit(int *icontxt, const char layout, const int *nprow, const int *npcol);
void blacs_gridinfo(int *icontxt, int *nprow, int *npcol, int *myprow, int *mypcol);
std::size_t numroc(std::size_t *n, std::size_t *nb, int *iproc, const int *isrcproc, int *nprocs);
void descinit(std::size_t *desc, std::size_t *m, std::size_t *n, std::size_t *mb, std::size_t *nb,
        const int *irsrc, const int *icsrc, int *ictxt, std::size_t *lld, int *info);

/////////////
//ScaLAPACK//
/////////////

template <typename T>
void t_pgeqrf(const std::size_t m, const std::size_t n, T *a, const std::size_t ia, const std::size_t ja, 
	      const std::size_t *desc_a, T *tau);

template <typename T>
void t_pmqr(const char side, const char trans, const std::size_t m, const std::size_t n, const std::size_t k,
	    const T *a, const std::size_t ia, const std::size_t ja, const std::size_t *desc_a, const T *tau,
	    T *c, const std::size_t ic, const std::size_t jc, const std::size_t *desc_c);

#include "blas_templates.inc"


