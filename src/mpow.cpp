
#include <cmath>

#include "marlib.hpp"

namespace marlib {

  template <typename MatrixT, typename IntegerT>
  dense_matrix<typename MatrixT::ValueType, typename MatrixT::RangeType>& mpow(const trans_t& trans,
    const MatrixT& MA, IntegerT m, dense_matrix<typename MatrixT::ValueType, typename MatrixT::RangeType>& ME) {

    using ValueT = typename MatrixT::ValueType;
    using RangeT = typename MatrixT::RangeType;

    assert(MA.nrow() == MA.ncol());
    size_type n = MA.nrow();
    dense_matrix<ValueT,RangeT> MX = dense_matrix<ValueT,RangeT>::eye(n);
    dense_matrix<ValueT,RangeT> MT(n,n);

    int info = 0;
    if (m < 0) {
      MT = MA;
      if (trans == Trans) {
        MT = MT.tr();
      }
      dgesv(MT, MX);
      m = -m;
    } else {
      MX = MA;
      if (trans == Trans) {
        MX = MX.tr();
      }
    }

    ME.eye();
    while (m != 0) {
      if (m % 2 == 1) {
        MT = ME;
        dgemm(NoTrans, NoTrans, ValueT(1), MX, MT, ValueT(0), ME);
      }
      m /= 2;
      dgemm(NoTrans, NoTrans, ValueT(1), MX, MX, ValueT(0), MT);
      MX = MT;
    }
    return ME;
  }

  template dense_matrix<double,int>& mpow(const trans_t& trans,
    const dense_matrix<double,int>& MA, int m, dense_matrix<double,int>& ME);
  template dense_matrix<double,int>& mpow(const trans_t& trans,
    const csr_matrix<double,int>& MA, int m, dense_matrix<double,int>& ME);

}
