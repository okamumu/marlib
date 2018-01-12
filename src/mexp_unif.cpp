/*
! Description: matrix exp with uniformization
!
!        y = exp(Q*t)
!
!        Q is uniformized to P and qv
!        t is involved in the Poisson probability vector.
*/

#include "marlib.hpp"

namespace marlib {

  template <typename ValueT, typename RangeT, typename MatrixT, typename VectorT>
  VectorT& mexp_unif(const trans_t& trans, const MatrixT& P, const ValueT& qv,
  const poisson<ValueT,RangeT>& pois, const VectorT& x, VectorT& y,
  const ValueT& atol) {
    VectorT xi = x.clone();
    VectorT tmp = x.clone();

    y = 0;
    daxpy(pois(pois.left()), xi, y);
    for (RangeT k=pois.left()+1; k<=pois.right(); k++) {
      tmp = xi;
      dgemm(trans, NoTrans, ValueT(1), P, tmp, ValueT(0), xi);
      daxpy(pois(k), xi, y);
      if (damax(xi) < atol)
        break;
    }
    y /= pois.weight();
    // dscal(1.0/pois.weight(), y);
    return y;
  }

  // dense

  template dense_matrix<double,int>& mexp_unif(const trans_t& trans,
    const dense_matrix<double,int>& P,
    const double& qv, const poisson<double,int>& pois,
    const dense_matrix<double,int>& x, dense_matrix<double,int>& y, const double& atol);

  template vector<double,int>& mexp_unif(const trans_t& trans,
    const dense_matrix<double,int>& P,
    const double& qv, const poisson<double,int>& pois,
    const vector<double,int>& x, vector<double,int>& y, const double& atol);

  // csr

  template dense_matrix<double,int>& mexp_unif(const trans_t& trans,
    const csr_matrix<double,int>& P,
    const double& qv, const poisson<double,int>& pois,
    const dense_matrix<double,int>& x, dense_matrix<double,int>& y, const double& atol);

  template vector<double,int>& mexp_unif(const trans_t& trans,
    const csr_matrix<double,int>& P,
    const double& qv, const poisson<double,int>& pois,
    const vector<double,int>& x, vector<double,int>& y, const double& atol);
}
