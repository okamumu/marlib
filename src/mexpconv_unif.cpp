/*
  ! Description: convolution integral operation for matrix exp form;
  !
  !                           |t
  ! transH(MH) = transH(MH) + | exp(transQ(Q)*s) * x * y' * exp(transQ(Q)*(t-s)) ds
  !                           |0
  !
  !        and
  !
  !        z = exp(transQ(Q)*t) * x
  !
  !        t is involved in the Poisson probability vector.
  !        qv is an uniformed parameter
  !        return value is z
 */

#include <cmath>

#include "marlib.hpp"

namespace marlib {

  template <typename ValueT, typename RangeT, typename MatrixT, typename VectorT>
  VectorT& mexpconv_unif(const trans_t& transQ, const trans_t& transH,
  const MatrixT& P, const ValueT qv, const poisson<ValueT,RangeT>& pois,
  const VectorT& x, const VectorT& y, VectorT& z, MatrixT& H) {

    const RangeT left = pois.left();
    const RangeT right = pois.right();

    assert(left == 0 && right >= 1);

    array<VectorT*> vc(right+1);
    size_type s = y.value().size();
    ValueT* vcpool = new ValueT [s * (right - left + 1)];
    ValueT* p = vcpool;
    for (RangeT l=left; l<=right; l++) {
      vc.ptr(l) = new VectorT(y, p);
      p += s;
    }

    vc[right] = 0;
    daxpy(pois(right), y, vc[right]);
    for (RangeT l=right-1; l>=left+1; l--) {
      dgemm(trans_c(transQ), NoTrans, ValueT(1), P, vc[l+1], ValueT(0), vc[l]);
      daxpy(pois(l), y, vc[l]);
    }

    VectorT tmp = x.clone();
    VectorT xi = x.clone();
    z = 0;
    H *= qv*pois.weight();
    daxpy(pois(left), xi, z);
    dger(transH, ValueT(1), xi, vc[left+1], H);

    for (RangeT l=left+1; l<=right-1; l++) {
      tmp = xi;
      dgemm(transQ, NoTrans, ValueT(1), P, tmp, ValueT(0), xi);
      daxpy(pois(l), xi, z);
      dger(transH, ValueT(1), xi, vc[l+1], H);
    }
    z /= pois.weight();
    H /= qv*pois.weight();

    for (RangeT l=left; l<=right; l++) {
      delete vc.ptr(l);
    }
    delete [] vcpool;

    return z;
  }

  template
  vector<double,int>& mexpconv_unif(const trans_t& transQ, const trans_t& transH,
  const dense_matrix<double,int>& P, const double qv, const poisson<double,int>& pois,
  const vector<double,int>& x, const vector<double,int>& y, vector<double,int>& z,
  dense_matrix<double,int>& H);

  template
  vector<double,int>& mexpconv_unif(const trans_t& transQ, const trans_t& transH,
  const csr_matrix<double,int>& P, const double qv, const poisson<double,int>& pois,
  const vector<double,int>& x, const vector<double,int>& y, vector<double,int>& z,
  csr_matrix<double,int>& H);

}
