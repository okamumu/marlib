/*
! Description: matrix power
!
!        ME = MA^m
!
*/


#include <cmath>

#include "marlib.hpp"

namespace marlib {

  /*
  !  SOR step for solving the following linear equation
  !
  !         alpha * trans(A - sigma I) * x = b
  !
  !         The step computes
  !
  !         fwd
  !         x := (D/omega + L)^(-1) (b/alpha - (U - D (1-omega)/omega - sigma I) * x)
  !
  !         bwd
  !         x := (D/omega + U)^(-1) (b/alpha - (L - D (1-omega)/omega - sigma I) * x)
  !
  !         A: square matrix
  !         x: vector (in; initial vector for the step, out; updated vector)
  !         b: constant vector
  */

  template <typename ValueT, typename RangeT>
  void gsstep_impl(const trans_t& trans, const ValueT alpha,
    dense_matrix<ValueT,RangeT> A,
    const ValueT sigma,
    const ValueT omega,
    dense_matrix<ValueT,RangeT> b, dense_matrix<ValueT,RangeT> x) {
    assert(A.nrow() == A.ncol());
    assert(b.ncol() == x.ncol());
    size_type n = A.nrow();
    size_type nrhs = b.ncol();
    A.set_range(range<RangeT>(1,n),range<RangeT>(1,n));
    b.set_range(range<RangeT>(1,n),range<RangeT>(1,nrhs));
    x.set_range(range<RangeT>(1,n),range<RangeT>(1,nrhs));

    vector<ValueT,RangeT> tmpx(nrhs);
    tmpx.set_range(range<RangeT>(1,nrhs));

    switch (trans) {
      case NoTrans:
        for (RangeT i=1; i<=n; i++) {
          ValueT tmpd = 0;
          for (RangeT u=1; u<=nrhs; u++) {
            tmpx(u) = b(i,u) / alpha;
          }
          for (RangeT j=1; j<=n; j++) {
            if (i == j) {
              tmpd = A(i,j);
              for (RangeT u=1; u<=nrhs; u++) {
                tmpx(u) += sigma * x(j,u);
              }
            } else {
              for (RangeT u=1; u<=nrhs; u++) {
                tmpx(u) -= A(i,j) * x(j,u);
              }
            }
          }
          for (RangeT u=1; u<=nrhs; u++) {
            x(i,u) = omega * tmpx(u) / tmpd + (1 - omega) * x(i,u);
          }
        }
        break;
      case Trans:
        for (RangeT j=1; j<=n; j++) {
          ValueT tmpd = 0;
          for (RangeT u=1; u<=nrhs; u++) {
            tmpx(u) = b(j,u) / alpha;
          }
          for (RangeT i=1; i<=n; i++) {
            if (i == j) {
              tmpd = A(i,j);
              for (RangeT u=1; u<=nrhs; u++) {
                tmpx(u) += sigma * x(i,u);
              }
            } else {
              for (RangeT u=1; u<=nrhs; u++) {
                tmpx(u) -= A(i,j) * x(i,u);
              }
            }
          }
          for (RangeT u=1; u<=nrhs; u++) {
            x(j,u) = omega * tmpx(u) / tmpd + (1 - omega) * x(j,u);
          }
        }
        break;
    }
  }

  template <typename ValueT, typename RangeT>
  void gsstep_impl(const trans_t& trans, const ValueT alpha,
    const csr_matrix<ValueT, RangeT>& A,
    const ValueT sigma,
    const ValueT omega,
    dense_matrix<ValueT,RangeT> b, dense_matrix<ValueT,RangeT> x) {
    assert(A.nrow() == A.ncol());
    assert(b.ncol() == x.ncol());
    size_type n = A.nrow();
    size_type nrhs = b.ncol();
    b.set_range(range<RangeT>(0,n-1),range<RangeT>(1,nrhs));
    x.set_range(range<RangeT>(0,n-1),range<RangeT>(1,nrhs));

    vector<ValueT,RangeT> tmpx(nrhs);
    tmpx.set_range(range<RangeT>(1,nrhs));

    switch (trans) {
      case NoTrans:
        for (size_type i=0; i<A.nrow(); i++) {
          ValueT tmpd = 0;
          for (RangeT u=1; u<=nrhs; u++) {
            tmpx(u) = b(i,u) / alpha;
          }
          for (size_type z=A.rowptr(i); z<A.rowptr(i+1); z++) {
            size_type j = A.colind(z);
            if (i == j) {
              tmpd = A.value(z);
              for (RangeT u=1; u<=nrhs; u++) {
                tmpx(u) += sigma * x(i,u);
              }
            } else {
              for (RangeT u=1; u<=nrhs; u++) {
                tmpx(u) -= A.value(z) * x(j,u);
              }
            }
          }
          for (RangeT u=1; u<=nrhs; u++) {
            x(i,u) = omega * tmpx(u) / tmpd + (1 - omega) * x(i,u);
          }
        }
        break;
      case Trans:
        dense_matrix<ValueT,RangeT> tmp(n,nrhs);
        vector<size_type,RangeT> diag(n);
        tmp.set_range(range<RangeT>(0,n-1),range<RangeT>(1,nrhs));
        diag.set_range(range<RangeT>(0,n-1));
        tmp = 0;
        daxpy(1/alpha, b, tmp);
        for (size_type i=0; i<n; i++) {
          for (size_type z=A.rowptr(i); z<A.rowptr(i+1); z++) {
            size_type j = A.colind(z);
            if (i == j) {
              diag(i) = z;
              for (RangeT u=1; u<=nrhs; u++) {
                tmp(i,u) += sigma * x(i,u);
              }
              break;
            } else {
              for (RangeT u=1; u<=nrhs; u++) {
                tmp(j,u) -= A.value(z) * x(i,u);
              }
            }
          }
        }
        for (size_type i=0; i<n; i++) {
          for (RangeT u=1; u<=nrhs; u++) {
            x(i,u) = omega * tmp(i,u) / A.value(diag(i)) + (1 - omega) * x(i,u);
          }
          for (size_type z=diag(i)+1; z<A.rowptr(i+1); z++) {
            size_type j = A.colind(z);
            for (RangeT u=1; u<=nrhs; u++) {
              tmp(j,u) -= A.value(z) * x(i,u);
            }
          }
        }
        break;
    }
  }

  template void gsstep_impl(const trans_t& trans, const double alpha,
    dense_matrix<double,int> A,
    const double sigma, const double omega,
    dense_matrix<double,int> b, dense_matrix<double,int> x);

  template void gsstep_impl(const trans_t& trans, const double alpha,
    const csr_matrix<double,int>& A,
    const double sigma, const double omega,
    dense_matrix<double,int> b, dense_matrix<double,int> x);
}
