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


////////////
// BLACS //
////////////
void FC_GLOBAL(blacs_pinfo, BLACS_PINFO)(BlasInt *mypnum, BlasInt *nprocs);
void FC_GLOBAL(blacs_get, BLACS_GET)(BlasInt *icontxt, const BlasInt *what, BlasInt *val );
void FC_GLOBAL(blacs_gridinit, BLACS_GRIDINIT)(BlasInt *icontxt, const char* layout,
				const BlasInt *nprow, const BlasInt *npcol);
void FC_GLOBAL(blacs_gridinfo, BLACS_GRIDINFO)(BlasInt *icontxt, BlasInt *nprow, BlasInt *npcol,
			        BlasInt *myprow, BlasInt *mypcol);

int FC_GLOBAL(numroc, NUMROC)(BlasInt *n, BlasInt *nb, BlasInt *iproc,
			      const BlasInt *isrcproc, BlasInt *nprocs);

void FC_GLOBAL(descinit, DESCINIT)(BlasInt *desc, BlasInt *m, BlasInt *n, BlasInt *mb, BlasInt *nb,
		const BlasInt *irsrc, const BlasInt *icsrc, BlasInt *ictxt, BlasInt *lld, BlasInt *info);


///////////////
// Scalapack //
//////////////

//QR factorization
void FC_GLOBAL(psgeqrf, PSGEQRF)(const BlasInt *m, const BlasInt *n,
				 float *a, const BlasInt *ia,
			         const BlasInt *ja, const BlasInt *desc_a, 
				 float *tau, float *work, 
				 BlasInt *lwork, BlasInt *info);

void FC_GLOBAL(pdgeqrf, PDGEQRF)(const BlasInt *m, const BlasInt *n,
                                 double *a, const BlasInt *ia,
                                 const BlasInt *ja, const BlasInt *desc_a,
                                 double *tau, double *work,
                                 BlasInt *lwork, BlasInt *info);

void FC_GLOBAL(pcgeqrf, PCGEQRF)(const BlasInt *m, const BlasInt *n,
                                 scomplex *a, const BlasInt *ia,
                                 const BlasInt *ja, const BlasInt *desc_a,
                                 scomplex *tau, scomplex *work,
                                 BlasInt *lwork, BlasInt *info);

void FC_GLOBAL(pzgeqrf, PZGEQRF)(const BlasInt *m, const BlasInt *n,
                                 dcomplex *a, const BlasInt *ia,
                                 const BlasInt *ja, const BlasInt *desc_a,
                                 dcomplex *tau, dcomplex *work,
                                 BlasInt *lwork, BlasInt *info);

//Multiply an orthogonal matrix produced by QR factorization formed by pxgeqrf
void FC_GLOBAL(psormqr, PSORMQR)(const char* side, const char *trans,
				 const BlasInt *m, const BlasInt *n,
                                 const BlasInt *k, const float *a,
				 const BlasInt *ia, const BlasInt *ja,
				 const BlasInt *desc_a, const float *tau,
				 float *c, const BlasInt *ic,
				 const BlasInt *jc, const BlasInt *desc_c,
				 float *work, BlasInt *lwork,
				 BlasInt *info);

void FC_GLOBAL(pdormqr, PDORMQR)(const char* side, const char *trans,
                                 const BlasInt *m, const BlasInt *n,
                                 const BlasInt *k, const double *a,
                                 const BlasInt *ia, const BlasInt *ja,
                                 const BlasInt *desc_a, const double *tau,
                                 double *c, const BlasInt *ic,
                                 const BlasInt *jc, const BlasInt *desc_c,
                                 double *work, BlasInt *lwork,
                                 BlasInt *info);

void FC_GLOBAL(pcunmqr, PCUNMQR)(const char* side, const char *trans,
                                 const BlasInt *m, const BlasInt *n,
                                 const BlasInt *k, const scomplex *a,
                                 const BlasInt *ia, const BlasInt *ja,
                                 const BlasInt *desc_a, const scomplex *tau,
                                 scomplex *c, const BlasInt *ic,
                                 const BlasInt *jc, const BlasInt *desc_c,
                                 scomplex *work, BlasInt *lwork,
                                 BlasInt *info);

void FC_GLOBAL(pzunmqr, PZUNMQR)(const char* side, const char *trans,
                                 const BlasInt *m, const BlasInt *n,
                                 const BlasInt *k, const dcomplex *a,
                                 const BlasInt *ia, const BlasInt *ja,
                                 const BlasInt *desc_a, const dcomplex *tau,
                                 dcomplex *c, const BlasInt *ic,
                                 const BlasInt *jc, const BlasInt *desc_c,
                                 dcomplex *work, BlasInt *lwork,
                                 BlasInt *info);

}

