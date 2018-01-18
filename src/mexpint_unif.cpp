/*
  ! Description: integral operation for matrix exp form;
  !
  !                    |t
  !        cME = cME + | exp(Q*s) ds
  !                    |0
  !
  !        ME = exp(Q*t)
  !
  !        Q is uniformized to P and qv
  !        t is involved in the Poisson probability vector.
  !        return value is ME
 */

#include <cmath>

#include "marlib.hpp"

namespace marlib {

  template <typename ValueT, typename RangeT, typename MatrixT, typename VectorT>
  VectorT& mexpint_unif(const trans_t& trans, const MatrixT& P, const ValueT qv,
  const poisson<ValueT,RangeT>& pois, const VectorT& x, VectorT& y, VectorT& cy) {

    const RangeT left = pois.left();
    const RangeT right = pois.right();
    vector<ValueT,RangeT> cpoi(right - left + 1);
    cpoi.set_range(range<RangeT>(left,right));

    cpoi(right) = 0;
    for (RangeT k=right-1; k>=left; k--) {
      cpoi(k) = cpoi(k+1) + pois(k+1);
    }

    VectorT xi = x.clone();
    VectorT tmp = x.clone();

    y = 0;
    cy *= qv*pois.weight();
    daxpy(pois(left), xi, y);
    daxpy(cpoi(left), xi, cy);
    for (RangeT k=left+1; k<=right; k++) {
      tmp = xi;
      dgemm(trans, NoTrans, ValueT(1), P, tmp, ValueT(0), xi);
      daxpy(pois(k), xi, y);
      daxpy(cpoi(k), xi, cy);
    }
    y /= pois.weight();
    cy /= qv*pois.weight();
    return y;
  }

  // dense_matrix

  template dense_matrix<double,int>& mexpint_unif(const trans_t& trans,
    const dense_matrix<double,int>& P, const double qv, const poisson<double,int>& pois,
    const dense_matrix<double,int>& x, dense_matrix<double,int>& y, dense_matrix<double,int>& cy);

  template vector<double,int>& mexpint_unif(const trans_t& trans,
    const dense_matrix<double,int>& P, const double qv, const poisson<double,int>& pois,
    const vector<double,int>& x, vector<double,int>& y, vector<double,int>& cy);

  // csr_matrix

  template dense_matrix<double,int>& mexpint_unif(const trans_t& trans,
    const csr_matrix<double,int>& P, const double qv, const poisson<double,int>& pois,
    const dense_matrix<double,int>& x, dense_matrix<double,int>& y, dense_matrix<double,int>& cy);

  template vector<double,int>& mexpint_unif(const trans_t& trans,
    const csr_matrix<double,int>& P, const double qv, const poisson<double,int>& pois,
    const vector<double,int>& x, vector<double,int>& y, vector<double,int>& cy);

  template <typename ValueT, typename RangeT, typename MatrixT, typename VectorT>
  VectorT& mexpint_unif(const trans_t& trans, const MatrixT& P, const ValueT qv,
  const poisson<ValueT,RangeT>& pois, const VectorT& x, VectorT& y, VectorT& cy, ValueT* p) {

    const RangeT left = pois.left();
    const RangeT right = pois.right();
    vector<ValueT,RangeT> cpoi(right - left + 1);
    cpoi.set_range(range<RangeT>(left,right));

    cpoi(right) = 0;
    for (RangeT k=right-1; k>=left; k--) {
      cpoi(k) = cpoi(k+1) + pois(k+1);
    }

    VectorT xi = x.clone(p);
    VectorT tmp = x.clone(p + x.value().size());

    y = 0;
    cy *= qv*pois.weight();
    daxpy(pois(left), xi, y);
    daxpy(cpoi(left), xi, cy);
    for (RangeT k=left+1; k<=right; k++) {
      tmp = xi;
      dgemm(trans, NoTrans, ValueT(1), P, tmp, ValueT(0), xi);
      daxpy(pois(k), xi, y);
      daxpy(cpoi(k), xi, cy);
    }
    y /= pois.weight();
    cy /= qv*pois.weight();
    return y;
  }

  // dense_matrix

  template dense_matrix<double,int>& mexpint_unif(const trans_t& trans,
    const dense_matrix<double,int>& P, const double qv, const poisson<double,int>& pois,
    const dense_matrix<double,int>& x, dense_matrix<double,int>& y, dense_matrix<double,int>& cy, double* p);

  template vector<double,int>& mexpint_unif(const trans_t& trans,
    const dense_matrix<double,int>& P, const double qv, const poisson<double,int>& pois,
    const vector<double,int>& x, vector<double,int>& y, vector<double,int>& cy, double* p);

  // csr_matrix

  template dense_matrix<double,int>& mexpint_unif(const trans_t& trans,
    const csr_matrix<double,int>& P, const double qv, const poisson<double,int>& pois,
    const dense_matrix<double,int>& x, dense_matrix<double,int>& y, dense_matrix<double,int>& cy, double* p);

  template vector<double,int>& mexpint_unif(const trans_t& trans,
    const csr_matrix<double,int>& P, const double qv, const poisson<double,int>& pois,
    const vector<double,int>& x, vector<double,int>& y, vector<double,int>& cy, double* p);
}
