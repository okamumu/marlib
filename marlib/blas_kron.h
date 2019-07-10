// #pragma once
//
// #include <marlib.h>
//
// namespace marlib {
//
//   /*
//   !
//   ! Description:
//   !   This routine computes the following kronecker product:
//   !
//   !     y = alpha (I_{n1} * tras(A) * I_{n2}) x + beta y
//   !
//   !   where A is an m-by-n square matrix,
//   !         I_{n} is an n-by-n identity matrix, and;
//   !         * is kronecker product.
//   !
//   ! Parameters:
//   !   m, n: rows and columns of matrix A
//   !   n1: size of identity matrix located on left
//   !   n2: size of identity matrix located on right
//   !
//   */
//
//   template <typename T1, typename T2, typename T3>
//   void dgemv(NOTRANS, double alpha, const T1& A, const T2& x, double beta, T3& y, KRONMatrixT) {
//     using traits1 = kron_matrix_traits<T1>;
//     using traits2 = vector_traits<T2>;
//     using traits3 = vector_traits<T3>;
//     const int n1 = traits1::left(A);
//     const int n2 = traits1::right(A);
//     const int m = traits1::nrow(A);
//     const int n = traits1::ncol(A);
//     const double* valueA = traits1::value(A);
//     const double* valueX = traits2::value(x);
//     double* valueY = traits3::value(y);
//
//     dscal(beta, y);
//     for (int i=0, j=0; i<n*n1*n2; i+=n*n2, j+=m*n2) {
//       const_dense_matrix X(n, n2, n2, &valueX[i]);
//       dense_matrix Y(m, n2, n2, &valueY[i]);
//       dgemm(TRANS(), NOTRANS(), alpha, traits1::genmat(A), X, 1.0, Y, traits1::MatrixT(), DenseMatrixT());
//     }
//   }
//
// void blas_dkronmv(
//     char trans,
// int m, int n, int n1, int n2,
// double alpha,
// const double *A, int lda,
// const double *x,
// double beta, double *y) {
//
// int i, j;
// if (trans == 'T' || trans == 't') {
// blas_dscal(n*n1*n2, beta, y, 1);
// for (i=0, j=0; i<m*n1*n2; i+=m*n2, j+=n*n2) {
// // blas_dgemm('N', 'N', n2, n, m,
// // 	alpha, &x[i], n2,
// // 	A, lda, 1.0, &y[j], n2);
// blas_dgemm('T', 'T', n, n2, m,
//            alpha, A, lda, &x[i], n2,
// 1.0, &y[j], n2);
// }
// } else if (trans == 'N' || trans == 'n') {
// blas_dscal(m*n1*n2, beta, y, 1);
// for (i=0, j=0; i<n*n1*n2; i+=n*n2, j+=m*n2) {
// // blas_dgemm('N', 'T', n2, m, n,
// // 	alpha, &x[i], n2,
// // 	A, lda, 1.0, &y[j], n2);
// blas_dgemm('T', 'N', m, n2, n,
//            alpha, A, lda, &x[i], n2,
// 1.0, &y[j], n2);
// }
// }
// }
// */
//
// /*
// !
// ! Description:
// !   This routine computes the following aggregation on kronecker product:
// !
// !     B = B + alpha x^T y
// !       and element-wise aggregation of the m-by-n matrix A; I_{n1} * A * I_{n2}
// !
// !   where A is m-by-n matrix
// !         I_{n} is an n-by-n identity matrix, and;
// !         * is kronecker product.
// !
// ! Parameters:
// !   m, n: row and column sizes of square matrix A
// !   n1: size of identity matrix located on left
// !   n2: size of identity matrix located on right
// !
// */
// /*
// void blas_dkronr(int m, int n, int n1, int n2,
//                  double alpha, double *x, double *y,
// double *A, int lda) {
// int i;
// for (i=0; i<n*n1*n2; i+=n*n2) {
// blas_dgemm('T', 'N', m, n, n2,
//            alpha, &x[i], n2,
// &y[i], n2, 1.0, A, lda);
// }
// }
//
// */
//
//   template <typename T1, typename T2, typename T3>
//   void dgemv(NOTRANS, double alpha, const T1& A, const T2& x, double beta, T3& y, CSRMatrixT) {
//     using traits1 = csr_matrix_traits<T1>;
//     using traits2 = vector_traits<T2>;
//     using traits3 = vector_traits<T3>;
//     const int base = traits1::base(A);
//     const int m = traits1::nrow(A);
//     // const int n = traits1::ncol(A);
//     const double* valueA = traits1::value(A);
//     const int* rowptr = traits1::rowptr(A);
//     const int* colind = traits1::colind(A);
//     const double* valueX = traits2::value(x);
//     const int incx = traits2::inc(x);
//     double* valueY = traits3::value(y);
//     const int incy = traits3::inc(y);
//
//     dscal(beta, y);
//     for (int i=0; i<m; i++) {
//       for (int z=rowptr[i]-base; z<rowptr[i+1]-base; z++) {
//         int j = colind[z] - base;
//         valueY[i*incy] += alpha * valueA[z] * valueX[j*incx];
//       }
//     }
//   }
//
//   template <typename T1, typename T2, typename T3>
//   void dgemv(TRANS, double alpha, const T1& A, const T2& x, double beta, T3& y, CSRMatrixT) {
//     using traits1 = csr_matrix_traits<T1>;
//     using traits2 = vector_traits<T2>;
//     using traits3 = vector_traits<T3>;
//     const int base = traits1::base(A);
//     const int m = traits1::nrow(A);
//     // const int n = traits1::ncol(A);
//     const double* valueA = traits1::value(A);
//     const int* rowptr = traits1::rowptr(A);
//     const int* colind = traits1::colind(A);
//     const double* valueX = traits2::value(x);
//     const int incx = traits2::inc(x);
//     double* valueY = traits3::value(y);
//     const int incy = traits3::inc(y);
//
//     dscal(beta, y);
//     for (int i=0; i<m; i++) {
//       for (int z=rowptr[i]-base; z<rowptr[i+1]-base; z++) {
//         int j = colind[z] - base;
//         valueY[j*incy] += alpha * valueA[z] * valueX[i*incx];
//       }
//     }
//   }
//
//   template <typename T1, typename T2, typename T3>
//   void dger(NOTRANS, double alpha, const T1& x, const T2& y, T3& A, CSRMatrixT) {
//     using traits1 = csr_matrix_traits<T3>;
//     using traits2 = vector_traits<T1>;
//     using traits3 = vector_traits<T2>;
//     const int base = traits1::base(A);
//     const int m = traits1::nrow(A);
//     // const int n = traits1::ncol(A);
//     double* valueA = traits1::value(A);
//     const int* rowptr = traits1::rowptr(A);
//     const int* colind = traits1::colind(A);
//     const double* valueX = traits2::value(x);
//     const int incx = traits2::inc(x);
//     const double* valueY = traits3::value(y);
//     const int incy = traits3::inc(y);
//
//     for (int i=0; i<m; i++) {
//       for (int z=rowptr[i]-base; z<rowptr[i+1]-base; z++) {
//         int j = colind[z] - base;
//         valueA[z] += alpha * valueX[i*incx] * valueY[j*incy];
//       }
//     }
//   }
//
//   template <typename T1, typename T2, typename T3>
//   void dger(TRANS, double alpha, const T1& x, const T2& y, T3& A, CSRMatrixT) {
//     dger(NOTRANS(), alpha, y, x, A, CSRMatrixT());
//   }
//
// }
