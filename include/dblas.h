/*
  blas.h
  wrapper for blas
 */

#ifndef _DBLAS_H_
#define _DBLAS_H_

namespace dblas {

  // level 1
  void dcopy(const int& n, const double *x, const int& incx, double *y, const int& incy);
  void dscal(const int& n, const double& alpha, double *x, const int& incx);
  void daxpy(const int& n, const double& alpha, const double *x, const int& incx,
    double *y, const int& incy);

  double ddot(const int& n, const double *x, const int& incx,
    const double *y, const int& incy);
  double dasum(const int& n, const double *x, const int& incx);
  double dnrm2(const int& n, const double *x, const int& incx);
  const double* idamax(const int& n, const double *x, const int& incx);

  // level 2
  void dgemv(const char& trans, const int& m, const int& n, const double& alpha,
    const double *A, const int& lda, const double *x, const int& incx,
    const double& beta, double *y, const int& incy);

  void dger(const int& m, const int& n, const double& alpha,
    const double *x, const int& incx, const double *y, const int& incy,
    double *A, const int& lda);

  // level 3
  void dgemm(const char& transA, const char& transB,
    const int& m, const int& n, const int& k, const double& alpha,
    const double *A, const int& lda, const double *B, const int& ldb,
    const double& beta, double *C, const int& ldc);

}

#endif
