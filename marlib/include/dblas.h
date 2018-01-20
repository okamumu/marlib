/*
  blas.h
  wrapper for blas
 */

namespace dblas {

  // level 1
  void dcopy(const int n, const double *x, const int incx, double *y, const int incy);
  void dscal(const int n, const double alpha, double *x, const int incx);
  void daxpy(const int n, const double alpha, const double *x, const int incx,
    double *y, const int incy);

  double ddot(const int n, const double *x, const int incx,
    const double *y, const int incy);
  double dasum(const int n, const double *x, const int incx);
  double dnrm2(const int n, const double *x, const int incx);
  const double* idamax(const int n, const double *x, const int incx);

  // level 2
  void dgemv(const char& trans, const int m, const int n, const double alpha,
    const double *A, const int lda, const double *x, const int incx,
    const double beta, double *y, const int incy);

  void dger(const int m, const int n, const double alpha,
    const double *x, const int incx, const double *y, const int incy,
    double *A, const int lda);

  // level 3
  void dgemm(const char& transA, const char& transB,
    const int m, const int n, const int k, const double alpha,
    const double *A, const int lda, const double *B, const int ldb,
    const double beta, double *C, const int ldc);

  // use for f77blas
  #ifdef F77BLAS

  extern "C" {
    void dcopy_(const int *n, const double *x, const int *incx, double *y, const int *incy);
    void dscal_(const int *n, const double *alpha, double *x, const int *incx);
    void daxpy_(const int *n, const double *alpha, const double *x, const int *incx, double *y, const int *incy);
    double ddot_(const int *n, const double *x, const int *incx, const double *y, const int *incy);
    double dnrm2_(const int *n, const double *x, const int *incx);
    double dasum_(const int *n, const double *x, const int *incx);
    void dgemv_(const char *trans, const int *m, const int *n, const double *alpha, const double *A, const int *lda,
      const double *x, const int *incx, const double *beta, double *y, const int *incy);
    void dger_(const int *m, const int *n, const double *alpha, const double *x, const int *incx,
      const double *y, const int *incy, double *A, const int *lda);
    void dgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k, const double *alpha,
      const double *A, const int *lda, const double *B, const int *ldb, const double *beta, double *C, const int *ldc);
  }

  #define __DCOPY__ dcopy_
  #define __DSCAL__ dscal_
  #define __DAXPY__ daxpy_
  #define __DDOT__ ddot_
  #define __DNRM2__ dnrm2_
  #define __DASUM__ dasum_
  #define __DGEMV__ dgemv_
  #define __DGER__ dger_
  #define __DGEMM__ dgemm_
  #endif

  #ifdef MKLBLAS

  extern "C" {
    void dcopy_(const int *n, const double *x, const int *incx, double *y, const int *incy);
    void dscal_(const int *n, const double *alpha, double *x, const int *incx);
    void daxpy_(const int *n, const double *alpha, const double *x, const int *incx, double *y, const int *incy);
    double ddot_(const int *n, const double *x, const int *incx, const double *y, const int *incy);
    double dnrm2_(const int *n, const double *x, const int *incx);
    double dasum_(const int *n, const double *x, const int *incx);
    void dgemv_(const char *trans, const int *m, const int *n, const double *alpha, const double *A, const int *lda,
      const double *x, const int *incx, const double *beta, double *y, const int *incy);
    void dger_(const int *m, const int *n, const double *alpha, const double *x, const int *incx,
      const double *y, const int *incy, double *A, const int *lda);
    void dgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k, const double *alpha,
      const double *A, const int *lda, const double *B, const int *ldb, const double *beta, double *C, const int *ldc);
  }

  #define __DCOPY__ dcopy_
  #define __DSCAL__ dscal_
  #define __DAXPY__ daxpy_
  #define __DDOT__ ddot_
  #define __DNRM2__ dnrm2_
  #define __DASUM__ dasum_
  #define __DGEMV__ dgemv_
  #define __DGER__ dger_
  #define __DGEMM__ dgemm_
  #endif

  // level 1

  /*
   * copy from x to y
   * @param n length of vector
   * @param x vector
   * @param incx increment of index of vector x
   * @param y vector
   * @param incy increment of index of vector y
   */

  inline void dcopy(const int n, const double *x, const int incx, double *y, const int incy) {
    __DCOPY__(&n, x, &incx, y, &incy);
  }

  /**
    * x = alpha * x
    * @param n length of vector x
    * @param alpha a value
    * @param x vector
    * @param incx increment of index of vector x
    */

  inline void dscal(const int n, const double alpha, double *x, const int incx) {
    __DSCAL__(&n, &alpha, x, &incx);
  }

  /*
   * y = a * x + y
   * @param n length of vector
   * @param x vector
   * @param incx increment of index of vector x
   * @param y vector
   * @param incy increment of index of vector y
   */

  inline void daxpy(const int n, const double alpha, const double *x, const int incx,
    double *y, const int incy) {
    __DAXPY__(&n, &alpha, x, &incx, y, &incy);
  }

  /*
   * dot product
   * @param n length of vector
   * @param x vector
   * @param incx increment of index of vector x
   * @param y vector
   * @param incy increment of index of vector y
   * @return a dot product of x and y
   */

  inline double ddot(const int n, const double *x, const int incx,
    const double *y, const int incy) {
    return __DDOT__(&n, x, &incx, y, &incy);
  }

  /*
   * l1norm
   * @param n length of vector
   * @param x vector
   * @param incx increment of index of vector x
   * @return a l2norm of x
   */

  inline double dasum(const int n, const double *x, const int incx) {
    return __DASUM__(&n, x, &incx);
  }

  /*
   * l2norm
   * @param n length of vector
   * @param x vector
   * @param incx increment of index of vector x
   * @return a l2norm of x
   */

  inline double dnrm2(const int n, const double *x, const int incx) {
    return __DNRM2__(&n, x, &incx);
  }

  // double blas_dsum(int n, const double *x, int incx) {
  //   int i;
  //   double sum = 0.0;
  //   for (i=0; i<n; i++, x+=incx) {
  //     sum += *x;
  //   }
  //   return sum;
  // }
  //
  // int blas_idmax(int n, const double *x, int incx) {
  //   int i, id;
  //   double max = *x;
  //   id = 0;
  //   for (i=0; i<n; i++, x+=incx) {
  //     if (*x > max) {
  //       max = *x;
  //       id = i;
  //     }
  //   }
  //   return id;
  // }
  //

  inline const double* idamax(const int n, const double *x, const int incx) {
    const double* max = x;
    for (int i=0; i<n; i++, x+=incx) {
      if (fabs(*x) > fabs(*max)) {
        max = x;
      }
    }
    return max;
  }

  // level 2

  /**
    * y = alpha * trans(A) * x + beta * y
    * @param trans a charactor to indicate the transpose of matrix A. If trans is 'T', the matrix A is transposed. Otherwise, if trans is 'N', the matrix A is not indicated.
    * @param m the number of rows of matrix A
    * @param n the number of columns of matrix A
    * @param alpha a value
    * @param A a matrix
    * @param lda length of column data of matrix A
    * @param x a vector
    * @param incx the increment of index of vector x
    * @param beta a value
    * @param y a vector
    * @param incy increment of index of vector y
    */

  inline void dgemv(const char& trans, const int m, const int n, const double alpha,
    const double *A, const int lda, const double *x, const int incx,
    const double beta, double *y, const int incy) {
    __DGEMV__(&trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy);
  }

  /**
    * A = alpha * x * y + A
    * @param m the number of rows of matrix A
    * @param n the number of columns of matrix A
    * @param alpha a value
    * @param x a vector
    * @param incx the increment of index of vector x
    * @param y a vector
    * @param incy increment of index of vector y
    * @param A a matrix
    * @param lda length of column data of matrix A
    */

  inline void dger(const int m, const int n, const double alpha,
    const double *x, const int incx, const double *y, const int incy,
    double *A, const int lda) {
    __DGER__(&m, &n, &alpha, x, &incx, y, &incy, A, &lda);
  }

  // level 3

  /**
    * C = alpha * transa(A) * transb(B) + beta * C
    * @param transa a charactor to indicate the transpose of matrix A. If trans is 'T', the matrix A is transposed. Otherwise, if trans is 'N', the matrix A is not indicated.
    * @param transb a charactor to indicate the transpose of matrix B. If trans is 'T', the matrix B is transposed. Otherwise, if trans is 'N', the matrix B is not indicated.
    * @param m the number of rows of matrix C
    * @param n the number of columns of matrix C
    * @param alpha a value
    * @param A a matrix
    * @param lda length of column data of matrix A
    * @param B a matrix
    * @param ldb length of column data of matrix B
    * @param beta a value
    * @param C a matrix
    * @param ldc length of column data of matrix C
    */

  inline void dgemm(const char& transA, const char& transB,
    const int m, const int n, const int k, const double alpha,
    const double *A, const int lda, const double *B, const int ldb,
    const double beta, double *C, const int ldc) {
    __DGEMM__(&transA, &transB, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
  }
}
