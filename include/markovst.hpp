#pragma once

namespace marlib {

  template <typename VectorT, typename MatrixT>
  void gth_impl(const MatrixT& Q, VectorT x);

  template <typename VectorT, typename MatrixT>
  VectorT& gth(const MatrixT& Q, VectorT& x) {
    gth_impl(Q, x);
    return x;
  }

}
