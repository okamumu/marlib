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
  void gsstep_impl(const trans_t& trans, const ValueT& alpha,
    dense_matrix<ValueT,RangeT> A,
    const ValueT& sigma,
    const ValueT& omega,
    vector<ValueT,RangeT> b, vector<ValueT,RangeT> x) {
    assert(A.nrow() == A.ncol());
    size_type n = A.nrow();
    A.set_range(range<RangeT>(1,n),range<RangeT>(1,n));
    b.set_range(range<RangeT>(1,n));
    x.set_range(range<RangeT>(1,n));

    switch (trans) {
      case NoTrans:
        for (RangeT i=1; i<=n; i++) {
          ValueT tmpd = 0;
          ValueT tmpx = b(i) / alpha;
          for (RangeT j=1; j<=n; j++) {
            if (i == j) {
              tmpd = A(i,j);
              tmpx += sigma * x(j);
            } else {
              tmpx -= A(i,j) * x(j);
            }
          }
          x(i) = omega * tmpx / tmpd + (1 - omega) * x(i);
        }
        break;
      case Trans:
        for (RangeT j=1; j<=n; j++) {
          ValueT tmpd = 0;
          ValueT tmpx = b(j) / alpha;
          for (RangeT i=1; i<=n; i++) {
            if (i == j) {
              tmpd = A(i,j);
              tmpx += sigma * x(i);
            } else {
              tmpx -= A(i,j) * x(i);
            }
          }
          x(j) = omega * tmpx / tmpd + (1 - omega) * x(j);
        }
        break;
    }
  }

  template <typename ValueT, typename RangeT>
  void gsstep_impl(const trans_t& trans, const ValueT& alpha,
    const csr_matrix<ValueT, RangeT>& A,
    const ValueT& sigma,
    const ValueT& omega,
    vector<ValueT,RangeT> b, vector<ValueT,RangeT> x) {
    assert(A.nrow() == A.ncol());
    size_type n = A.nrow();
    b.set_range(range<RangeT>(0,n-1));
    x.set_range(range<RangeT>(0,n-1));

    switch (trans) {
      case NoTrans:
        for (size_type i=0; i<A.nrow(); i++) {
          ValueT tmpd = 0;
          ValueT tmpx = b(i) / alpha;
          for (size_type z=A.rowptr(i); z<A.rowptr(i+1); z++) {
            size_type j = A.colind(z);
            if (i == j) {
              tmpd = A.value(z);
              tmpx += sigma * x(i);
            } else {
              tmpx -= A.value(z) * x(j);
            }
          }
          x(i) = omega * tmpx / tmpd + (1 - omega) * x(i);
        }
        break;
      case Trans:
        vector<ValueT,RangeT> tmp(n);
        vector<size_type,RangeT> diag(n);
        tmp.set_range(range<RangeT>(0,n-1));
        diag.set_range(range<RangeT>(0,n-1));
        tmp = 0;
        daxpy(1/alpha, b, tmp);
        for (size_type i=0; i<n; i++) {
          for (size_type z=A.rowptr(i); z<A.rowptr(i+1); z++) {
            size_type j = A.colind(z);
            if (i == j) {
              diag(i) = z;
              tmp(i) += sigma * x(i);
              break;
            } else {
              tmp(j) -= A.value(z) * x(i);
            }
          }
        }
        for (size_type i=0; i<n; i++) {
          x(i) = omega * tmp(i) / A.value(diag(i)) + (1 - omega) * x(i);
          for (size_type z=diag(i)+1; z<A.rowptr(i+1); z++) {
            size_type j = A.colind(z);
            tmp(j) -= A.value(z) * x(i);
          }
        }
        break;
    }
  }

  template void gsstep_impl(const trans_t& trans, const double& alpha,
    dense_matrix<double,int> A,
    const double& sigma, const double& omega,
    vector<double,int> b, vector<double,int> x);

  template void gsstep_impl(const trans_t& trans, const double& alpha,
    const csr_matrix<double,int>& A,
    const double& sigma, const double& omega,
    vector<double,int> b, vector<double,int> x);
}
