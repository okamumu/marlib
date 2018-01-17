
#include <cmath>

#include "marlib.hpp"

namespace marlib {

  template <typename ValueT, typename RangeT, typename MatrixT>
  dense_matrix<ValueT,RangeT>& mexp_pade(const trans_t& trans,
    const MatrixT& MA, const ValueT& t, dense_matrix<ValueT,RangeT>& ME, const ValueT& eps) {
    size_type n = MA.nrow();

    ME = MA;
    ME *= t;
    ValueT norma = damax(ME);

    RangeT j = static_cast<RangeT>(log(norma) / log(2.0));
    j = (0 < 1+j) ? 1+j : 0;
    ME /= std::ldexp(1.0,j);
    // dscal(1.0/std::ldexp(1.0,j), ME);

    RangeT q = 1;
    ValueT tolerr = 1.0 / 6.0;
    while (tolerr > eps/norma) {
      tolerr /= 16.0 * (3.0 + 4.0 * q * (2.0 + q));
      q++;
    }

    ValueT c = 1.0;
    dense_matrix<ValueT,RangeT> MD = dense_matrix<ValueT,RangeT>::eye(n);
    dense_matrix<ValueT,RangeT> MN = dense_matrix<ValueT,RangeT>::eye(n);
    dense_matrix<ValueT,RangeT> MX = dense_matrix<ValueT,RangeT>::eye(n);
    dense_matrix<ValueT,RangeT> MT = dense_matrix<ValueT,RangeT>::zeros(n,n);

    for (RangeT k=1; k<=q; k++) {
      c *= (q - k + 1.0) / ((2.0 * q - k + 1.0) * k);
      dgemm(trans, NoTrans, ValueT(1), ME, MX, ValueT(0), MT);
      MX = MT;
      daxpy(c, MX, MN);
      if (k % 2 == 0) {
        daxpy(c, MX, MD);
      } else {
        daxpy(-c, MX, MD);
      }
    }
    ME = MN;
    dgesv(MD, ME);
    for (RangeT k=1; k<=j; k++) {
      dgemm(NoTrans, NoTrans, ValueT(1), ME, ME, ValueT(0), MT);
      ME = MT;
    }
    return ME;
  }

  template dense_matrix<double,int>& mexp_pade(const trans_t& trans,
    const dense_matrix<double,int>& MA, const double& t, dense_matrix<double,int>& ME, const double& eps);
  template dense_matrix<double,int>& mexp_pade(const trans_t& trans,
    const csr_matrix<double,int>& MA, const double& t, dense_matrix<double,int>& ME, const double& eps);

}
