/*
! Description: GTH
*/


#include <cmath>

#include "marlib.hpp"

namespace marlib {

  template <typename VectorT, typename MatrixT>
  void gth_impl(const MatrixT& Q, VectorT x) {
    assert(Q.nrow() == Q.ncol());
    using ValueT = typename VectorT::ValueType;
    using RangeT = typename VectorT::RangeType;

    size_type n= Q.nrow();
    dense_matrix<ValueT,RangeT> A(n,n);
    A.set_range(range<RangeT>(1,n), range<RangeT>(1,n));
    x.set_range(range<RangeT>(1,n));
    A = Q;

    for (RangeT l=n; l>=2; l--) {
      const range<RangeT> ll(l,l);
      const range<RangeT> rr(1,l-1);
      ValueT tmp = dasum(A(ll,rr));
      for (RangeT j=1; j<=l-1; j++) {
        for (RangeT i=1; i<=l-1; i++) {
          if (i != j) {
            A(i,j) += A(l,j) * A(i,l) / tmp;
          }
        }
      }
      A(rr,ll) /= tmp;
      A(ll,rr) = 0;
      A(l,l) = -1;
    }

    x(1) = 1;
    for (RangeT l=2; l<=n; l++) {
      const range<RangeT> rr(1,l-1);
      x(l) = ddot(x(rr),A(rr,l));
    }
    x /= dasum(x);
  }

  template vector<double,int>& gth(const dense_matrix<double,int>& Q, vector<double,int>& x);
  template vector<double,int>& gth(const csr_matrix<double,int>& Q, vector<double,int>& x);

  template <typename ValueT, typename RangeT, typename MatrixT, typename ConfigT>
  vector<ValueT,RangeT>& stgs(const MatrixT& Q, vector<ValueT,RangeT>& x, ConfigT& conf) {
    size_type n = Q.nrow();
    vector<ValueT,RangeT> prevx(n);
    vector<ValueT,RangeT> b(n);
    b = 0;

    conf.iter = 0;
    while(true) {
      prevx = x;
      for (size_type k=0; k<conf.steps; k++) {
        gsstep(Trans, ValueT(1), Q, ValueT(0), ValueT(1), b, x);
        x /= dasum(x);
      }
      conf.iter += conf.steps;
      prevx -= x;
      conf.aerror = damax(prevx);
      conf.rerror = conf.aerror / damax(x);

      conf.callback(x);

      if (conf.rerror < conf.rtol && conf.aerror < conf.atol) {
        conf.info = 0; // convergence
        break;
      }

      if (conf.iter >= conf.maxiter) {
        conf.info = 1; // maxiter
        break;
      }
    }
    return x;
  }

  template vector<double,int>& stgs(const dense_matrix<double,int>& Q, vector<double,int>& xinit, gsconf<double,vector<double,int>>& conf);
}
