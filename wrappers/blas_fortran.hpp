#pragma once

#include <fortran_mangle.h>
#include <complex>

using BlasInt = int;
using dcomplex = std::complex<double>;
using scomplex = std::complex<float>;

extern "C" {

		
void FC_GLOBAL(slatms, SLATMS)(const BlasInt* m,    const BlasInt* n,
			       const char* dist,    BlasInt* iseed,
			       const char* sym,     float* d,
			       const BlasInt* mode, const float* cond,
			       const float* dmax,   const BlasInt* kl,
			       const BlasInt* ku,   const char* pack,
			       float*a,             const BlasInt* lda,  
			       float* work,   	    BlasInt* info);

void FC_GLOBAL(dlatms, DLATMS)(const BlasInt* m,    const BlasInt* n,
                               const char* dist,    BlasInt* iseed,
                               const char* sym,     double* d,
                               const BlasInt* mode, const double* cond,
                               const double* dmax,  const BlasInt* kl,
                               const BlasInt* ku,   const char* pack,
                               double *a,           const BlasInt* lda,
                               double* work,  	    BlasInt* info);

void FC_GLOBAL(clatms, CLATMS)(const BlasInt* m,    const BlasInt* n,
                               const char* dist,    BlasInt* iseed,
                               const char* sym,     float* d,
                               const BlasInt* mode, const float* cond,
                               const float* dmax,   const BlasInt* kl,
                               const BlasInt* ku,   const char* pack,
                               scomplex *a,         const BlasInt* lda,
                               scomplex* work,      BlasInt* info);

void FC_GLOBAL(zlatms, ZLATMS)(const BlasInt* m,    const BlasInt* n,
                               const char* dist,    BlasInt* iseed,
                               const char* sym,     double* d,
                               const BlasInt* mode, const double* cond,
                               const double* dmax,  const BlasInt* kl,
                               const BlasInt* ku,   const char* pack,
                               dcomplex*a,          const BlasInt* lda,
                               dcomplex* work,      BlasInt* info);


}

